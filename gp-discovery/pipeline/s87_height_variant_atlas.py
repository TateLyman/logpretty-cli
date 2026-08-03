"""
Stage 87 - height-increasing human variant atlas.

The previous human-genetics branch (stages 41-47 and 83) worked from monogenic
syndromes and ClinVar. That is the wrong instrument for this question: it finds genes
where losing the protein makes a child ill, and misses the quantitative alleles that
make otherwise healthy adults a centimetre taller. This stage looks for the second
kind.

One rule governs the whole atlas, and it is the brief's:

    a positional MAPPED_GENE assignment is NOT causal evidence.

The GWAS Catalog's gene search returns everything near a locus. Querying it for STC2
returns intergenic variants whose nearest features are RNA pseudogenes 5 kb away. Every
variant here therefore carries a `gene_assignment_basis` column, and only variants whose
consequence in the named gene is protein-altering - confirmed by Ensembl VEP against the
transcript, not by distance - are marked causal-grade.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import allelelib as A  # noqa: E402
import gputil as G  # noqa: E402

R = G.RESULTS

# The seed is curated, and that is stated. It spans the target classes the brief
# prioritises - secreted inhibitors, extracellular proteases, binding proteins,
# surface receptors, ion transporters, matrix regulators - plus the genes any height
# atlas has to contain. Membership in the seed confers nothing; every gene is then
# tested against the catalogue and against VEP.
SEED = {
    "secreted inhibitor": ["STC1", "STC2", "NOG", "GREM1", "CHRD", "SOST", "IGFBP1",
                           "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7",
                           "IGFALS", "FSTL3"],
    "extracellular protease": ["PAPPA", "PAPPA2", "ADAMTS17", "ADAMTS10", "ADAMTS6",
                               "ADAMTS3", "ADAMTSL3", "MMP13", "MMP14", "BMP1", "TLL1",
                               "FURIN", "PCSK5", "HTRA1"],
    "binding protein / matrix": ["ACAN", "FBN1", "FBN2", "LTBP1", "LTBP2", "LTBP3",
                                 "EFEMP1", "MATN3", "CHSY1", "COL11A1", "COL27A1",
                                 "HSPG2", "SPARC"],
    "cell-surface receptor": ["IGF1R", "IGF2R", "INSR", "GHR", "NPR2", "NPR3", "FGFR3",
                              "BMPR1B", "ACVR1", "PTH1R", "LIFR", "ROR2", "CRLF1"],
    "ligand / local hormone": ["IGF1", "IGF2", "GH1", "NPPC", "GDF5", "BMP2", "BMP6",
                               "IHH", "PTHLH", "CNMD"],
    "ion transporter / enzyme": ["SLC39A8", "SLC26A2", "TRIP11", "DOT1L", "SUZ12",
                                 "GNA12"],
    "transcriptional / other": ["ZBTB38", "LCORL", "HMGA1", "HMGA2", "SHOX", "SHOX2"],
}
GENE_CLASS = {g: c for c, gs in SEED.items() for g in gs}
GENES = sorted(GENE_CLASS)

# The catalogue's gene search is paged. A first version of this stage used size=120,
# which silently truncated 69 of the 77 genes - including STC2, whose only two
# protein-altering variants (rs148833559, rs2277923) sit past position 120 and were
# therefore invisible to the atlas while STC2 was being called "no coding variant".
# The largest gene here returns 747 records, so one page of 1000 retrieves everything
# and the cap stops being a hidden filter.
MAX_SNPS_PER_GENE = 1000
MAX_CODING_PER_GENE = 60


def gene_snps(gene: str) -> list[dict]:
    j = A.jget(f"{A.GWAS}/singleNucleotidePolymorphisms/search/findByGene"
               f"?geneName={gene}&size={MAX_SNPS_PER_GENE}", "s87gs")
    return (j.get("_embedded") or {}).get("singleNucleotidePolymorphisms") or []


def snp_assoc(rsid: str) -> list[dict]:
    j = A.jget(f"{A.GWAS}/singleNucleotidePolymorphisms/{rsid}/associations"
               "?projection=associationBySnp", "s87as")
    return (j.get("_embedded") or {}).get("associations") or []


def study_of(assoc: dict) -> dict:
    """Ancestry, sample size and provenance for one association.

    The brief asks for ancestry and replication per variant. Neither lives on the
    association record - both live on the study it links to, so the link is followed.
    """
    href = ((assoc.get("_links") or {}).get("study") or {}).get("href", "")
    if not href:
        return {}
    j = A.jget(href.replace("{?projection}", ""), "s87st")
    anc = j.get("ancestries") or []
    groups, n_init, countries = set(), 0, set()
    for a in anc:
        for g in a.get("ancestralGroups") or []:
            if g.get("ancestralGroup"):
                groups.add(g["ancestralGroup"])
        for c in a.get("countryOfRecruitment") or []:
            if c.get("countryName"):
                countries.add(c["countryName"])
        if a.get("type") == "initial":
            try:
                n_init += int(a.get("numberOfIndividuals") or 0)
            except (TypeError, ValueError):
                pass
    return {
        "ancestry": "; ".join(sorted(groups)) or "not stated",
        "countries_of_recruitment": "; ".join(sorted(countries))[:80],
        "initial_sample_size": j.get("initialSampleSize", ""),
        "replication_sample_size": j.get("replicationSampleSize", ""),
        "n_individuals_initial": n_init or np.nan,
        "study_accession": j.get("accessionId", ""),
        "pubmed_id": (j.get("publicationInfo") or {}).get("pubmedId", ""),
        "study_trait": j.get("diseaseTrait", {}).get("trait", "")
                       if isinstance(j.get("diseaseTrait"), dict) else "",
    }


def vep(rsid: str) -> dict:
    # hgvs=1 is what makes VEP return hgvsp, i.e. the actual protein change. Without
    # it the atlas records "missense_variant" but not which residue - and the residue
    # is what a structure-guided modality search in stage 90 needs.
    j = A.jget(f"{A.ENSEMBL}/vep/human/id/{rsid}?content-type=application/json&hgvs=1",
               "s87vep")
    return j[0] if isinstance(j, list) and j else {}


def main() -> None:
    G.log(f"stage 87: {len(GENES)} seed genes across {len(SEED)} target classes")

    snps = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(gene_snps, g): g for g in GENES}
        for f in as_completed(futs):
            snps[futs[f]] = f.result()
    n_all = sum(len(v) for v in snps.values())

    # ---- the positional trap, quantified ----------------------------------
    ctx_rows = []
    coding_candidates = {}
    for gene, ss in snps.items():
        seen = set()
        n_cod = 0
        for s in ss:
            rs = s.get("rsId")
            if not rs or rs in seen:
                continue
            seen.add(rs)
            fc = str(s.get("functionalClass") or "")
            ctx_rows.append({"seed_gene": gene, "rsid": rs, "functional_class": fc,
                             "is_coding_class": fc in A.CODING})
            if fc in A.CODING and n_cod < MAX_CODING_PER_GENE:
                coding_candidates[rs] = gene
                n_cod += 1
    ctx = pd.DataFrame(ctx_rows)
    n_uniq = ctx.rsid.nunique()
    n_cod = int(ctx.is_coding_class.sum())
    G.log(f"   {n_all:,} SNP records, {n_uniq:,} distinct variants near the seed genes; "
          f"{n_cod:,} have a coding functional class ({n_cod / max(n_uniq, 1):.1%})")

    # ---- confirm the gene by transcript consequence, not by distance -------
    veps = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(vep, rs): rs for rs in coding_candidates}
        for i, f in enumerate(as_completed(futs), 1):
            veps[futs[f]] = f.result()
            if i % 60 == 0:
                G.log(f"   VEP {i}/{len(coding_candidates)}")

    # ---- associations, restricted to height traits -------------------------
    assocs = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(snp_assoc, rs): rs for rs in coding_candidates}
        for i, f in enumerate(as_completed(futs), 1):
            assocs[futs[f]] = f.result()
            if i % 60 == 0:
                G.log(f"   associations {i}/{len(coding_candidates)}")

    rows = []
    for rs, seed_gene in coding_candidates.items():
        v = veps.get(rs) or {}
        tcs = v.get("transcript_consequences") or []
        # the gene the variant actually alters, per transcript
        alt = [t for t in tcs
               if set(t.get("consequence_terms") or []) & A.PROTEIN_ALTERING]
        alt_genes = sorted({t.get("gene_symbol") for t in alt if t.get("gene_symbol")})
        hgvsp = next((t.get("hgvsp") for t in alt if t.get("hgvsp")), "")
        sift = next((t.get("sift_prediction") for t in alt if t.get("sift_prediction")), "")
        polyphen = next((t.get("polyphen_prediction") for t in alt
                         if t.get("polyphen_prediction")), "")
        most_severe = v.get("most_severe_consequence", "")
        causal_grade = bool(seed_gene in alt_genes)

        for a in assocs.get(rs, []):
            traits = [t.get("trait", "") for t in a.get("efoTraits", [])]
            if not any(A.is_height_trait(t) for t in traits):
                continue
            loci = a.get("loci") or [{}]
            ra = (loci[0].get("strongestRiskAlleles") or [{}])[0]
            rep_genes = sorted({g.get("geneName") for g in
                                loci[0].get("authorReportedGenes") or []
                                if g.get("geneName")})
            beta, direction = a.get("betaNum"), a.get("betaDirection")
            unit = str(a.get("betaUnit") or "")
            # riskFrequency is populated at the ASSOCIATION level; the copy inside
            # strongestRiskAlleles is usually null. Reading only the latter reported
            # "frequency not reported" for 111 of 116 rows when the number was one
            # level up - and frequency is what separates a rare large-effect allele
            # from a common small-effect one.
            freq = np.nan
            for src in (a.get("riskFrequency"), ra.get("riskFrequency")):
                try:
                    freq = float(src)
                    break
                except (TypeError, ValueError):
                    continue
            # The catalogue's betaUnit is very often the literal string "unit",
            # which means the depositor did not state a unit. An earlier version of
            # this block matched "unit" into the SD branch and the report then
            # printed "0.335 SD" for ADAMTS10 - a number this pipeline invented.
            # Only an explicitly named unit is converted; everything else stays
            # unitless and is printed as such.
            cm = np.nan
            sd = np.nan
            unit_known = False
            if beta is not None and unit:
                if re.search(r"\bcm\b|centimet", unit, re.I):
                    cm, unit_known = float(beta), True
                elif re.search(r"\bsd\b|standard deviation|z.?score", unit, re.I):
                    sd, unit_known = float(beta), True
            # riskAlleleName is "rs62621197-C"; the trailing token is the allele the
            # beta is oriented to. "?" means the depositor did not state it, and a
            # direction without an allele cannot be compared across studies.
            ra_name = str(ra.get("riskAlleleName", ""))
            eff_allele = ra_name.split("-")[-1].strip() if "-" in ra_name else ""
            if eff_allele in ("", "?"):
                eff_allele = ""
            rows.append({
                "rsid": rs, "seed_gene": seed_gene,
                "target_class": GENE_CLASS.get(seed_gene, ""),
                "risk_allele": ra.get("riskAlleleName", ""),
                "allele_frequency": freq,
                "frequency_class": ("rare (<1%)" if np.isfinite(freq) and freq < 0.01
                                    else "low-frequency (1-5%)"
                                    if np.isfinite(freq) and freq < 0.05
                                    else "common (>=5%)" if np.isfinite(freq)
                                    else "not reported"),
                "beta": beta, "beta_unit": unit,
                "beta_unit_is_stated": unit_known,
                "beta_direction": direction,
                "effect_allele": eff_allele,
                "effect_cm": cm, "effect_sd": sd,
                "increases_height": bool(direction == "increase"),
                "height_increasing_allele": (eff_allele if direction == "increase"
                                             else ""),
                "pvalue": a.get("pvalue"),
                "traits": "; ".join(traits[:3]),
                "gwas_catalog_functional_class":
                    next((s.get("functionalClass") for s in snps.get(seed_gene, [])
                          if s.get("rsId") == rs), ""),
                "vep_most_severe_consequence": most_severe,
                "vep_protein_altering_genes": "; ".join(alt_genes),
                "vep_hgvsp": hgvsp or "",
                "sift": sift, "polyphen": polyphen,
                "author_reported_genes": "; ".join(rep_genes),
                **study_of(a),
                "gene_assignment_basis":
                    ("VEP protein-altering consequence in the seed gene - CAUSAL-GRADE"
                     if causal_grade else
                     "positional / mapped only - NOT causal evidence"),
                "causal_grade_gene_assignment": causal_grade,
            })
    at = pd.DataFrame(rows)
    if len(at):
        at = at.drop_duplicates(["rsid", "traits", "beta", "study_accession"])

    # ---- per-gene evidence: literature, animal, phenotype breadth ----------
    def gene_evidence(gene: str) -> dict:
        q = f'"{gene}"'
        return {
            "gene": gene,
            "epmc_height_records": A.epmc_count(
                f'{q} AND (height OR stature OR "body length")'),
            "epmc_adult_height_records": A.epmc_count(
                f'{q} AND ("adult height" OR "final height" OR "attained height")'),
            "epmc_functional_assay_records": A.epmc_count(
                f'{q} AND ("loss of function" OR "gain of function" OR hypomorph* OR '
                '"functional characterization" OR "functional characterisation")'),
            "epmc_knockout_records": A.epmc_count(
                f'{q} AND (knockout OR knock-out OR "null mice" OR "-/-" OR transgenic)'),
            "epmc_bone_length_records": A.epmc_count(
                f'{q} AND ("bone length" OR "femur length" OR "tibia length" OR '
                '"long bone" OR "growth plate")'),
            "epmc_dysplasia_records": A.epmc_count(
                f'{q} AND (dysplasia OR disproportionate OR "skeletal abnormalit*")'),
            "epmc_cancer_records": A.epmc_count(f'{q} AND (cancer OR tumour OR tumor)'),
            "epmc_vascular_records": A.epmc_count(
                f'{q} AND (vascular OR atherosclerosis OR aneurysm)'),
            "epmc_metabolic_records": A.epmc_count(
                f'{q} AND (insulin OR glucose OR "insulin resistance" OR hypoglyc*)'),
            "epmc_neuro_records": A.epmc_count(
                f'{q} AND ("intellectual disability" OR "developmental delay" OR '
                'seizure OR neurodevelopment*)'),
        }
    ev = {}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(gene_evidence, g): g for g in GENES}
        for i, f in enumerate(as_completed(futs), 1):
            ev[futs[f]] = f.result()
            if i % 20 == 0:
                G.log(f"   gene evidence {i}/{len(GENES)}")
    gev = pd.DataFrame(ev.values())

    if len(at):
        at = at.merge(gev, left_on="seed_gene", right_on="gene", how="left").drop(
            columns=["gene"])
        at["adult_vs_childhood_height"] = np.where(
            at.traits.str.contains("adult|final", case=False, na=False),
            "adult height explicitly", "height trait, age stratum not specified")
        at["proportionality_evidence"] = np.where(
            at.epmc_dysplasia_records.fillna(0)
            > at.epmc_bone_length_records.fillna(0),
            "gene's literature is dominated by dysplasia - proportionality NOT established",
            "no dysplasia dominance; proportionality still not directly measured")
        at["functional_assay_evidence"] = np.where(
            at.epmc_functional_assay_records.fillna(0) >= 5,
            "functional characterisation literature exists for this gene",
            "little or no functional characterisation")
        at["replication"] = np.where(
            at.groupby("rsid").rsid.transform("count") > 1,
            "variant appears in more than one catalogued association",
            "single catalogued association")
        # Direction in the catalogue is stated relative to each study's own effect
        # allele. Two studies that report opposite alleles for the same variant will
        # report opposite directions and still agree. Reconciling per (variant,
        # allele) is what makes a direction comparable across studies - without it,
        # rs2277950 in STC2 looks self-contradictory when it is not.
        def _orient(sub: pd.DataFrame) -> str:
            pairs = {(r.effect_allele, r.beta_direction) for r in sub.itertuples()
                     if r.effect_allele and r.beta_direction}
            if not pairs:
                return "effect allele or direction not stated - DIRECTION UNORIENTABLE"
            per_allele = {}
            for al, d in pairs:
                per_allele.setdefault(al, set()).add(d)
            if any(len(v) > 1 for v in per_allele.values()):
                return "same allele reported in both directions - CONFLICT"
            if len(per_allele) == 1:
                return "single effect allele, direction consistent"
            return ("opposite alleles reported in opposite directions - "
                    "consistent after orientation")
        at["variant_direction_consistency"] = at.groupby("rsid", group_keys=False)[
            at.columns.tolist()].apply(lambda s: pd.Series(
                [_orient(s)] * len(s), index=s.index))
        at["direction_orientable"] = ~at.variant_direction_consistency.str.contains(
            "UNORIENTABLE|CONFLICT")
        at = at.sort_values(["causal_grade_gene_assignment", "increases_height",
                             "pvalue"], ascending=[False, False, True])
    at.to_csv(R / "height_increasing_variant_atlas.csv", index=False)
    ctx.to_csv(R / "height_variant_gene_assignment_audit.csv", index=False)
    gev.to_csv(R / "height_gene_literature_evidence.csv", index=False)

    n_h = len(at)
    n_inc = int(at.increases_height.sum()) if len(at) else 0
    n_causal = int(at.causal_grade_gene_assignment.sum()) if len(at) else 0
    n_causal_inc = int((at.causal_grade_gene_assignment & at.increases_height).sum()) \
        if len(at) else 0
    G.log(f"atlas: {n_h} height associations on coding-class variants; {n_inc} "
          f"height-increasing; {n_causal} causal-grade gene assignments; "
          f"{n_causal_inc} both")

    # ---- report ------------------------------------------------------------
    def _effect(r) -> str:
        """Print the effect in the unit the catalogue actually stated.

        The catalogue's betaUnit is frequently the literal word "unit", meaning no
        unit was recorded. Rendering that as "SD" would put a number in the report
        that no source supports, so an unstated unit is printed as unstated.
        """
        if np.isfinite(r.effect_cm):
            return f"{r.effect_cm:.3g} cm"
        if np.isfinite(r.effect_sd):
            return f"{r.effect_sd:.3g} SD"
        if r.beta is None or (isinstance(r.beta, float) and not np.isfinite(r.beta)):
            return "not reported"
        return f"{r.beta:.3g} (unit not stated by depositor)"

    L = ["# Quantitative height genetics report", "",
         "## Why this stage exists", "",
         "Stages 41-47 and 83 worked from monogenic syndromes, OMIM and ClinVar. That "
         "instrument answers 'which gene, when broken, makes a child ill' - and stage 83's "
         "answer was that of 38 genes with stature phenotypes, **none** produces proportionate "
         "tall stature without a cost. The instrument was wrong for the question. Quantitative "
         "alleles that make otherwise healthy adults measurably taller are invisible to it, "
         "because they cause no disease and therefore appear in no disease database.", "",
         "## The rule that governs the atlas", "",
         "> **A positional MAPPED_GENE assignment is not causal evidence.**", "",
         "This is not a formality. The GWAS Catalog's own gene search returns everything near "
         f"a locus: querying it for the {len(GENES)} seed genes returned {n_uniq:,} distinct "
         f"variants, of which **{n_cod:,} ({n_cod / max(n_uniq, 1):.1%}) have a coding "
         "functional class at all**. The rest are intergenic or intronic and their gene labels "
         "are statements about distance. For STC2 specifically, the first variants the search "
         "returns are intergenic, with RNA pseudogenes as their nearest features.", "",
         "Every row in the atlas therefore carries `gene_assignment_basis`, and only variants "
         "whose consequence is **protein-altering in the named gene, confirmed by Ensembl VEP "
         "against the transcript**, are marked `causal_grade_gene_assignment = True`.", "",
         "## What was assembled", "", "| step | count |", "|---|---:|",
         f"| seed genes (curated, spanning {len(SEED)} target classes) | {len(GENES)} |",
         f"| SNP records returned by the catalogue for those genes | {n_all:,} |",
         f"| distinct variants | {n_uniq:,} |",
         f"| variants with a coding functional class | {n_cod:,} |",
         f"| of those, carried forward to VEP + association lookup | "
         f"{len(coding_candidates):,} |",
         f"| height-trait associations found on them | {n_h} |",
         f"| **height-INCREASING** associations | **{n_inc}** |",
         f"| with a causal-grade gene assignment | {n_causal} |",
         f"| **both height-increasing and causal-grade** | **{n_causal_inc}** |", "",
         "The seed is curated and confers nothing. Membership does not make a gene a "
         "candidate; every gene is tested against the catalogue and against VEP, and the "
         "columns record which tests it passed.", ""]

    if n_causal_inc:
        L += ["## Height-increasing variants with a causal-grade gene assignment", "",
              "| variant | gene | class | allele | frequency | effect | direction | p | "
              "consequence | protein change |",
              "|---|---|---|---|---|---:|---|---:|---|---|"]
        for _, r in at[at.causal_grade_gene_assignment & at.increases_height].head(
                30).iterrows():
            eff = _effect(r)
            L.append(f"| `{r.rsid}` | **{r.seed_gene}** | {r.target_class} | "
                     f"{r.risk_allele} | "
                     f"{'—' if not np.isfinite(r.allele_frequency) else f'{r.allele_frequency:.4g}'} "
                     f"({r.frequency_class}) | {eff} | {r.beta_direction} | "
                     f"{r.pvalue:.0e} | {r.vep_most_severe_consequence} | "
                     f"{str(r.vep_hgvsp).split(':')[-1] or '—'} |")
        L.append("")
    else:
        L += ["## Height-increasing variants with a causal-grade gene assignment", "",
              "**None in the GWAS Catalog for these genes.** That is a finding about the "
              "catalogue rather than about human biology: the catalogue stores lead variants "
              "from published associations, and the rare coding-variant height studies the "
              "brief names - the exome analyses and burden tests - report their results as "
              "gene-level burden or in supplementary tables that the REST API does not "
              "expose. The variants exist in the literature; they are not retrievable as "
              "structured records here, and that limitation is stated rather than worked "
              "around by asserting effect sizes from memory.", ""]

    if n_h:
        L += ["## All catalogued height associations on coding-class variants", "",
              "| variant | seed gene | assignment basis | direction | effect | frequency | "
              "trait |", "|---|---|---|---|---:|---|---|"]
        for _, r in at.head(28).iterrows():
            eff = _effect(r)
            L.append(f"| `{r.rsid}` | {r.seed_gene} | "
                     f"{'**causal-grade**' if r.causal_grade_gene_assignment else 'positional'} | "
                     f"{r.beta_direction or '—'} | {eff} | {r.frequency_class} | "
                     f"{str(r.traits)[:40]} |")
        L.append("")

    L += ["## Per-gene evidence breadth", "",
          "Literature counts per gene, used in stage 88 to weigh direction and in stage 93 to "
          "weigh liability. Counts are counts: they say a paper exists, not what it found.", "",
          "| gene | class | height | adult height | functional | knockout | bone length | "
          "dysplasia | cancer | vascular | metabolic | neuro |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for _, r in gev.sort_values("epmc_adult_height_records",
                                ascending=False).head(30).iterrows():
        L.append(f"| **{r.gene}** | {GENE_CLASS.get(r.gene, '')} | "
                 f"{r.epmc_height_records:,} | {r.epmc_adult_height_records:,} | "
                 f"{r.epmc_functional_assay_records:,} | {r.epmc_knockout_records:,} | "
                 f"{r.epmc_bone_length_records:,} | {r.epmc_dysplasia_records:,} | "
                 f"{r.epmc_cancer_records:,} | {r.epmc_vascular_records:,} | "
                 f"{r.epmc_metabolic_records:,} | {r.epmc_neuro_records:,} |")

    L += ["", "## What this atlas does not contain, and why", "",
          "- **The GIANT 5.4-million-person analysis, the UK Biobank exome results and the "
          "2024 rare non-coding study are not ingested as primary data.** Their per-variant "
          "effect sizes live in supplementary spreadsheets attached to journal articles, not "
          "in any queryable API. The GWAS Catalog holds their lead SNPs but not the coding "
          "burden results, and downloading and parsing dozens of publisher-hosted "
          "supplementary files is not something this stage can do reliably. What is here is "
          "the catalogued, machine-readable subset, and it is labelled as such.",
          "- **No effect size is quoted from memory.** Where the catalogue has a beta and a "
          "unit, both are recorded; where it does not, the cell is empty. A remembered "
          "'+0.5 cm' would be indistinguishable in this file from a retrieved one, which is "
          "why none appears.",
          "- **Ancestry is not resolved per variant.** The association records reference "
          "study accessions; resolving each study's ancestry breakdown is a further call per "
          "study and was not made.",
          "- **Heterozygous versus homozygous effects are not separated.** GWAS betas are "
          "per-allele under an additive model by default; a recessive or dominance component "
          "is not recoverable from the catalogue record.",
          "- **Proportionality is inferred, not measured.** No catalogued record contains "
          "sitting height, leg length or segment ratios for these variants; the "
          "`proportionality_evidence` column is a literature-balance heuristic and says so.",
          "", "## What it does establish", "",
          f"The positional trap is real and large: **{100 - 100 * n_cod / max(n_uniq, 1):.0f}% "
          "of the variants a gene-name search returns for these genes are not coding at all**. "
          "Any pipeline that took the catalogue's gene labels at face value - which is what a "
          "generic enrichment or a mapped-gene ranking does - would be building on positional "
          "coincidence. Stage 88 works only from the causal-grade subset and from experimental "
          "direction evidence, never from mapping.", ""]
    (R / "quantitative_height_genetics_report.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
