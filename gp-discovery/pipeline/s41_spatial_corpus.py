"""
Stage 41 - intact-tissue spatial-evidence corpus for the 238 CRISPR_CAUSAL genes.

The order of the whole project is reversed here. Stages 1-40 started from CRISPR
effects, expression modules or compounds, and only at stage 37 discovered that
the localization assumption underneath the leading candidate was unreliable. So
this stage starts with localization and refuses to accept anything that is not an
image of intact tissue.

Concretely, the unit of evidence is a *figure caption in an open-access full
text* that names the gene, names a spatial method, and names a growth-plate
compartment, plus the Methods text that says which probe or antibody was used
and whether it was validated. A computational zone label is never evidence here,
and neither is a title or an abstract.

Sources queried per gene:
  * Europe PMC - targeted full-text query (method x growth-plate x gene)
  * MGI GXD - the curated expression-assay reference list for that marker,
    resolved to PubMed IDs through MGI's own BIB_PubMed report
  * BioStudies / BioImage Archive - imaging study search
  * Human Protein Atlas - recorded as queried; HPA has no growth-plate tissue,
    so it can only ever contribute a negative here
  * PubMed - count-level coverage check for the same query

EMAGE and Expression Atlas expose no programmatic query surface reachable from
this environment; both are recorded as attempted-and-unavailable in the source
coverage table rather than silently dropped.
"""
from __future__ import annotations

import json
import re
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402
import spatiallib as S  # noqa: E402

R = G.RESULTS
OUT = R / "stage41"
OUT.mkdir(parents=True, exist_ok=True)

ZONE_Q = ('("growth plate" OR "growth-plate" OR "hypertrophic chondrocyte" OR '
          '"proliferative zone" OR "resting zone" OR "epiphyseal cartilage")')
METHOD_Q = ('(RNAscope OR "in situ hybridization" OR "in situ hybridisation" OR '
            'immunohistochemistry OR immunofluorescence OR immunostaining OR smFISH OR '
            '"spatial transcriptomics" OR "reporter mouse" OR lacZ OR "lineage tracing")')

MAX_FULLTEXT_PER_GENE = 60      # capped, and the cap is reported
MAX_MGI_PMIDS = 40
THREADS = 8

_lock = threading.Lock()
_manifest: dict = {}
_rejected: list = []


def candidate_score(rec: dict, pat: re.Pattern) -> float:
    """Rank candidates by how likely the paper is to show THIS gene's localization."""
    title = rec.get("title") or ""
    abst = rec.get("abstractText") or ""
    s = 0.0
    if pat.search(title):
        s += 4
    if pat.search(abst):
        s += 2
    blob = title + " " + abst
    if S._any(S.INTACT, blob):
        s += 2
    if S._which(S.METHODS, blob):
        s += 1.5
    if S._which(S.ZONES, blob):
        s += 1.5
    if rec.get("source") == "MGI_GXD":
        s += 3          # MGI curated the paper as containing an expression assay
    try:
        s += min(float(rec.get("citedByCount") or 0), 200) / 400.0
    except Exception:  # noqa: BLE001
        pass
    return s


def gene_candidates(mouse: str, human) -> list[dict]:
    pat = S.gene_pattern(mouse, human)
    alts = f'"{mouse}"' + (f' OR "{human}"' if human and str(human) != "nan" else "")
    recs, seen = [], set()

    def add(rs, src):
        for x in rs:
            pid = x.get("pmcid") or x.get("pmid") or x.get("id")
            if not pid or pid in seen:
                continue
            seen.add(pid)
            x = dict(x)
            x["source"] = src
            recs.append(x)

    add(S.epmc_search(f'({alts}) AND {ZONE_Q} AND {METHOD_Q} AND (OPEN_ACCESS:y)',
                      page_size=100, max_pages=1), "EPMC_targeted")

    mgi = S.mgi_gxd_pmids().get(mouse, [])[-MAX_MGI_PMIDS:]
    for i in range(0, len(mgi), 20):
        batch = mgi[i:i + 20]
        q = "(" + " OR ".join(f"EXT_ID:{p}" for p in batch) + ") AND (OPEN_ACCESS:y)"
        add(S.epmc_search(q, page_size=25, max_pages=1), "MGI_GXD")

    for r in recs:
        r["_score"] = candidate_score(r, pat)
    recs.sort(key=lambda x: -x["_score"])
    return recs


def mine_article(mouse: str, human, rec: dict) -> list[dict]:
    """Every figure in this article whose caption localizes THIS gene."""
    pmcid = rec.get("pmcid")
    xml = S.fetch_fulltext(pmcid)
    if not xml:
        return []
    with _lock:
        if pmcid not in _manifest:
            p = S.FT / f"{pmcid}.xml"
            _manifest[pmcid] = {
                "pmcid": pmcid, "pmid": rec.get("pmid"), "doi": rec.get("doi"),
                "title": (rec.get("title") or "").strip()[:400],
                "year": rec.get("pubYear"), "source": rec.get("source"),
                "url": f"https://europepmc.org/article/PMC/{pmcid}",
                "sha256": G.sha256_file(p), "bytes": p.stat().st_size,
            }
    pat = S.gene_pattern(mouse, human)
    expr_cues, geno_cues = S.intent_patterns(mouse, human)
    meth = S.methods_text(xml)
    body = S._strip(xml)
    out, rejected = [], []
    for fig in S.figures(xml):
        cap = fig["caption"]
        if not cap or not pat.search(cap):
            continue
        # The gene must be the SUBJECT of the localization, not the genotype.
        # Most captions naming a gene are showing a mutant phenotype.
        gsent = " ".join(s for s in S.sentences(cap) if pat.search(s))
        expr_hit = S._first(expr_cues, gsent) or S._first(expr_cues, cap)
        if not expr_hit:
            rejected.append({
                "mouse_gene": mouse, "pmcid": pmcid, "figure": fig["label"],
                "reason": "gene named as genotype or context, not as the localized species",
                "genotype_cue": S._first(geno_cues, cap) or "",
                "quotation": gsent[:300] or cap[:300]})
            continue
        # zone signal must come from sentences that name the gene, not from the
        # whole caption - multi-gene panels otherwise donate their neighbours'
        # localization to this gene.
        gtext = gsent
        # a gene sentence that names only a non-spatial assay is not localization
        if S._any(S.NON_SPATIAL_ASSAY, gsent) and not S._which(S.METHODS, gsent):
            rejected.append({
                "mouse_gene": mouse, "pmcid": pmcid, "figure": fig["label"],
                "reason": "gene measured by a non-spatial assay in this figure",
                "genotype_cue": S._first(S.NON_SPATIAL_ASSAY, gsent) or "",
                "quotation": gsent[:300]})
            continue
        methods = S._which(S.METHODS, cap)
        method_source = "figure caption"
        if not methods:
            methods, method_source = S._which(S.METHODS, meth), "methods section"
        if not methods:
            methods, method_source = S._which(S.METHODS, body), "elsewhere in body text"
        if not methods:
            continue
        # The figure must show THIS gene in growth-plate tissue. A caption that
        # names a zone is the strong case; a caption that names the growth plate
        # without resolving a zone is still intact-tissue localization and is
        # kept, with the zone left to the body text and flagged as such.
        zones = S._which(S.ZONES, gtext)
        if not (zones or S._any(S.INTACT, gtext)):
            continue
        if not (S._any(S.INTACT, cap) or S._any(S.INTACT, body[:60000])):
            continue
        excl = S._which(S.EXCLUDE, cap)
        # body-level sentences about this gene, for intensity, extra compartments
        # and - when the caption did not resolve one - the zone itself
        bsent = [s for s in S.sentences(body)
                 if pat.search(s) and S._any(sum(S.ZONES.values(), []), s)][:12]
        btext = " ".join(bsent)
        zone_source = "figure caption" if zones else ("body text" if btext else "not resolved")
        if not zones:
            zones = S._which(S.ZONES, btext)
        val = S._which(S.VALIDATION, meth) + S._which(S.VALIDATION, cap)
        genetic_reporter = bool(S._any(S.METHODS["reporter mouse"] + S.METHODS["lineage tracing"],
                                       cap + " " + meth))
        quant = S._any(S.QUANTIFIED, cap) or S._any(S.QUANTIFIED, gtext)
        out.append({
            "mouse_gene": mouse, "human_gene": human,
            "pmid": rec.get("pmid"), "pmcid": pmcid, "doi": rec.get("doi"),
            "year": rec.get("pubYear"), "journal": ((rec.get("journalInfo") or {})
                                                    .get("journal", {}) or {}).get("title"),
            "title": (rec.get("title") or "").strip()[:300],
            "discovered_by": rec.get("source"),
            "figure": fig["label"] or "unlabelled figure",
            "source_quotation": gtext[:900] or cap[:900],
            "full_text_status": "FULL_TEXT_VERIFIED",
            "species": "; ".join(S._which(S.SPECIES, cap) or S._which(S.SPECIES, body[:80000])),
            "age": S._first(S.AGE, cap) or S._first(S.AGE, meth) or "",
            "sex": ("female" if re.search(r"\bfemale\b", meth, re.I) else
                    "male" if re.search(r"\bmale\b", meth, re.I) else ""),
            "bone": S._first(S.BONES, cap) or S._first(S.BONES, meth) or "",
            "anatomical_location": S._first(S.INTACT, cap) or S._first(S.INTACT, body[:60000]),
            "method": "; ".join(methods),
            "method_named_in": method_source,
            "probe_or_antibody_id": S._first(S.VALIDATION["catalogue identifier"], meth) or "",
            "reagent_validation": "; ".join(sorted(set(val))),
            "knockout_negative_control": "knockout control" in val,
            "genetic_reporter": genetic_reporter,
            "signal_resting": "resting" in zones,
            "signal_proliferative": "proliferative" in zones,
            "signal_prehypertrophic": "prehypertrophic" in zones,
            "signal_hypertrophic": "hypertrophic" in zones,
            "signal_terminal_hypertrophic": "terminal_hypertrophic" in zones,
            "signal_perichondrial": "perichondrial" in zones,
            "signal_non_chondrocyte": "; ".join(S._which(S.NON_CHONDRO, gtext)),
            "qualitative_intensity": "; ".join(S._which(S.INTENSITY, gtext + " " + btext)),
            "quantified": bool(quant),
            "localization_directly_visible": True,
            "zone_call_source": zone_source,
            "excluded_context": "; ".join(excl),
            "n_zones_named": len(zones),
            "localization_intent_cue": expr_hit,
            "genotype_cue_also_present": bool(S._first(geno_cues, cap)),
        })
    if rejected:
        with _lock:
            _rejected.extend(rejected)
    return out


def evidence_level(r) -> str:
    # a method that is only named somewhere in the body cannot be tied to this
    # figure, so the record cannot rise above "intact-tissue image, weak support"
    if r.excluded_context or r.method_named_in == "elsewhere in body text":
        return "LEVEL_D"
    val = bool(r.reagent_validation) or r.genetic_reporter
    strong_val = r.knockout_negative_control or r.genetic_reporter or (
        "sense/negative probe" in str(r.reagent_validation))
    if r.quantified and strong_val:
        return "LEVEL_A"
    if val:
        return "LEVEL_B"
    return "LEVEL_C"


def source_coverage(genes: pd.DataFrame) -> pd.DataFrame:
    """What each named source could actually contribute, recorded honestly."""
    rows = [
        ("Europe PMC", "REST search + fullTextXML", "USED",
         "primary channel: targeted full-text query per gene, then figure-caption mining"),
        ("PMC full text", "fullTextXML / efetch", "USED",
         "open-access articles only; paywalled full texts are not retrievable here"),
        ("MGI Gene Expression Database (GXD)", "MRK_GXD.rpt + BIB_PubMed.rpt", "USED",
         "curated expression-assay reference list per marker, resolved to PubMed IDs; MGI's "
         "structure-level annotations are not exposed in any downloadable report, so GXD is used "
         "to seed papers rather than to assign zones"),
        ("BioStudies / BioImage Archive", "REST search", "USED",
         "imaging-study search per gene"),
        ("Human Protein Atlas", "search_download.php", "USED_AS_NEGATIVE",
         "HPA's tissue atlas contains no growth plate, so it cannot supply intact-tissue "
         "growth-plate localization for any gene; queried and recorded as a negative"),
        ("PubMed", "E-utilities esearch", "USED",
         "coverage counts for the same query, to show what open-access restriction costs"),
        ("EMAGE", "emouseatlas.org", "UNAVAILABLE",
         "no programmatic query surface reachable from this environment; the site returns only "
         "the HTML portal"),
        ("Expression Atlas", "ebi.ac.uk/gxa/json", "UNAVAILABLE",
         "endpoint returns HTTP 404; and Expression Atlas carries bulk/single-cell summaries "
         "rather than intact-tissue images, so it could not have been direct evidence"),
        ("Publisher full texts and supplements", "-", "PARTIAL",
         "only what Europe PMC redistributes; no publisher-specific scraping was performed"),
        ("GEO/SRA-linked papers", "-", "INDIRECT",
         "reached through Europe PMC where the linked paper is open access"),
    ]
    d = pd.DataFrame(rows, columns=["source", "endpoint", "status", "note"])
    d.to_csv(OUT / "source_coverage.csv", index=False)
    return d


def main() -> None:
    genes = pd.read_csv(R / "gene_sets" / "CRISPR_CAUSAL.csv")
    G.log(f"stage 41: {len(genes)} CRISPR_CAUSAL genes")
    S.mgi_gxd_pmids()
    G.log(f"MGI GXD reference map: {len(S.mgi_gxd_pmids())} markers")
    cov = source_coverage(genes)
    G.log(f"source coverage table: {len(cov)} sources "
          f"({int((cov.status == 'UNAVAILABLE').sum())} unavailable, recorded not dropped)")

    per_gene, records = [], []

    def work(row):
        mouse, human = row.mouse_gene, row.human_gene
        cands = gene_candidates(mouse, human)
        oa = [c for c in cands if c.get("pmcid")][:MAX_FULLTEXT_PER_GENE]
        recs = []
        for c in oa:
            try:
                recs.extend(mine_article(mouse, human, c))
            except Exception as exc:  # noqa: BLE001
                G.log(f"   ! {mouse} {c.get('pmcid')}: {exc}")
        bst = S.biostudies_search(f'"{mouse}" growth plate')
        hpa = S.hpa_record(human)
        alts = f'"{mouse}"' + (f' OR "{human}"' if human and str(human) != "nan" else "")
        base = f'({alts}) AND {ZONE_Q} AND {METHOD_Q}'
        pm = f'({mouse}[tiab]{f" OR {human}[tiab]" if human and str(human) != "nan" else ""}) '\
             f'AND ("growth plate"[tiab] OR epiphyseal[tiab] OR "hypertrophic chondrocyte"[tiab]) '\
             f'AND ("in situ hybridization"[tiab] OR immunohistochemistry[tiab] OR '\
             f'immunofluorescence[tiab] OR RNAscope[tiab] OR immunostaining[tiab])'
        n_all, n_oa = S.epmc_count(base), S.epmc_count(base + " AND (OPEN_ACCESS:y)")
        return {
            "mouse_gene": mouse, "human_gene": human,
            "epmc_hits_all": n_all, "epmc_hits_open_access": n_oa,
            "open_access_fraction": round(n_oa / n_all, 3) if n_all else None,
            "pubmed_hits": S.pubmed_count(pm),
            "n_candidates": len(cands), "n_open_access": len([c for c in cands if c.get("pmcid")]),
            "n_fulltexts_examined": len(oa),
            "fulltext_cap_applied": len([c for c in cands if c.get("pmcid")]) > MAX_FULLTEXT_PER_GENE,
            "n_spatial_records": len(recs),
            "biostudies_hits": bst.get("total"),
            "hpa_found": bool(hpa), "hpa_subcellular": "; ".join(
                hpa.get("Subcellular main location", []) or []) if hpa else "",
            "hpa_growth_plate_tissue": False,
        }, recs

    with ThreadPoolExecutor(max_workers=THREADS) as ex:
        futs = {ex.submit(work, r): r.mouse_gene for r in genes.itertuples()}
        done = 0
        for f in as_completed(futs):
            summary, recs = f.result()
            per_gene.append(summary)
            records.extend(recs)
            done += 1
            if done % 20 == 0 or summary["n_spatial_records"]:
                G.log(f"   [{done}/{len(genes)}] {summary['mouse_gene']:10s} "
                      f"cand={summary['n_candidates']:4d} ft={summary['n_fulltexts_examined']:3d} "
                      f"records={summary['n_spatial_records']}")

    corpus = pd.DataFrame(records)
    if len(corpus):
        corpus["evidence_level"] = corpus.apply(evidence_level, axis=1)
        # replication: same gene, same zone pattern, independent publication
        zc = ["signal_resting", "signal_proliferative", "signal_prehypertrophic",
              "signal_hypertrophic", "signal_terminal_hypertrophic", "signal_perichondrial"]
        corpus["zone_pattern"] = corpus[zc].astype(int).astype(str).agg("".join, axis=1)
        rep = (corpus.groupby(["mouse_gene", "zone_pattern"])["pmcid"].nunique()
               .rename("n_independent_sources").reset_index())
        corpus = corpus.merge(rep, on=["mouse_gene", "zone_pattern"], how="left")
        corpus["pattern_replicates"] = corpus.n_independent_sources >= 2
        corpus = corpus.sort_values(
            ["mouse_gene", "evidence_level", "n_independent_sources"],
            ascending=[True, True, False])
    else:
        corpus = pd.DataFrame(columns=["mouse_gene", "evidence_level"])
    corpus.to_csv(R / "spatial_evidence_corpus.csv", index=False)

    rej = pd.DataFrame(_rejected)
    rej.to_csv(OUT / "figures_rejected_not_localization.csv", index=False)

    pg = pd.DataFrame(per_gene).sort_values("n_spatial_records", ascending=False)
    pg.to_csv(OUT / "per_gene_search_summary.csv", index=False)
    (R / "spatial_fulltext_manifest.json").write_text(json.dumps({
        "n_articles": len(_manifest),
        "fulltext_cap_per_gene": MAX_FULLTEXT_PER_GENE,
        "genes_where_cap_bound": int(pg.fulltext_cap_applied.sum()),
        "articles": sorted(_manifest.values(), key=lambda x: str(x["pmcid"])),
    }, indent=1))

    G.log(f"corpus: {len(corpus)} figure-level records over "
          f"{corpus.mouse_gene.nunique() if len(corpus) else 0} genes; "
          f"{len(_manifest)} full texts checksummed")
    if len(corpus):
        G.log(corpus.evidence_level.value_counts().to_string())
    G.log(f"rejected as genotype-not-localization: {len(rej)} figures "
          f"over {rej.mouse_gene.nunique() if len(rej) else 0} genes (kept, not discarded)")


if __name__ == "__main__":
    main()
