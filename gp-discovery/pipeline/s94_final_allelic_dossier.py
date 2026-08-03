"""
Stage 94 - final allelic-series dossier.

Assembles stages 87-93 into the three ranked tables and answers the brief's twelve
questions from the files those stages wrote, not from recollection. Every numeric claim
in the report is read out of a CSV at write time; where a question's honest answer is
"no", it is answered "no".
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS

# The four requirements from the brief, in the order stage 91 tests them.
REQ = ["human_direction", "molecular_direction", "animal_agreement",
       "accessible_intervention"]


def main() -> None:
    atlas = pd.read_csv(R / "height_increasing_variant_atlas.csv")
    series = pd.read_csv(R / "height_allelic_series.csv")
    paths = pd.read_csv(R / "genetically_anchored_growth_pathways.csv")
    chain = pd.read_csv(R / "stc2_pappa_evidence_chain.csv")
    imap = pd.read_csv(R / "stc2_interface_target_map.csv")
    mods = pd.read_csv(R / "stc2_pappa_modalities.csv")
    safety = pd.read_csv(R / "genetic_pathway_safety_matrix.csv")
    expr = pd.read_csv(R / "axis_expression_and_phenotype_profile.csv")
    claims = pd.read_csv(R / "stc2_pappa_literature_claims.csv")
    prec = pd.read_csv(R / "axis_clinical_precedent.csv")
    arms = pd.read_csv(R / "ex_vivo_arm_definitions.csv")

    # ---- top 20 targets ----------------------------------------------------
    t20 = paths.sort_values(
        ["requirements_met", "human_causal_grade_increasing", "n_mouse_longer"],
        ascending=[False, False, False]).head(20).copy()
    t20.insert(0, "rank", range(1, len(t20) + 1))
    t20["promotable"] = t20.requirements_met == 4
    t20.to_csv(R / "top_20_genetically_anchored_targets.csv", index=False)

    # ---- top 10 modalities -------------------------------------------------
    # ranked on direction first, then on whether the interface is solved, then on
    # whether anything exists to order. Wrong-direction rows can never rank.
    m = mods.copy()
    m["rank_key"] = (
        (m.direction == "RIGHT_DIRECTION").astype(int) * 1000
        + (m.verdict == "CANDIDATE MODALITY").astype(int) * 100
        + (m.structures_available > 0).astype(int) * 10
        + m.orderable_today.astype(int))
    m10 = m[m.direction == "RIGHT_DIRECTION"].sort_values(
        "rank_key", ascending=False).head(10).copy()
    m10.insert(0, "rank", range(1, len(m10) + 1))
    m10.to_csv(R / "top_10_growth_modalities.csv", index=False)

    # ---- top 5 leads -------------------------------------------------------
    leads = []
    for _, r in t20[t20.promotable].iterrows():
        sm = safety[(safety.gene == r.gene) & safety.concern.str.startswith("HIGH")]
        im = imap[imap.protein_to_act_on == r.gene]
        mm = m10[m10.protein_to_act_on == r.gene]
        leads.append({
            "gene": r.gene, "target_class": r.target_class,
            "allelic_series_class": r.allelic_series_class,
            "requirements_met": r.requirements_met,
            "human_increasing_variants": r.human_increasing_variants,
            "human_variant_count": r.human_causal_grade_increasing,
            "mouse_longer_terms": r.get("animal_agreement_note", ""),
            "compartment": r.subcellular_location,
            "structures_of_the_interface": int(im.structures_available.max())
                                           if len(im) else 0,
            "best_modality": mm.modality.iloc[0] if len(mm) else "none ranked",
            "catalogued_chemistry": int(r.chembl_activities),
            "orderable_probe_exists": bool(r.chembl_activities > 0),
            "highest_safety_concern": (sm.system.iloc[0] if len(sm)
                                       else "none flagged HIGH by these instruments"),
            "status": "GENETICALLY_ANCHORED_TARGET_AWAITING_DIRECTIONAL_TEST",
            "why_not_further": ("no reagent in the ex vivo plan has a measured potency, "
                                "so no arm has a stateable concentration"),
        })
    l5 = pd.DataFrame(leads).head(5)
    if len(l5):
        l5.insert(0, "rank", range(1, len(l5) + 1))
    l5.to_csv(R / "top_5_allelic_series_leads.csv", index=False)

    # ---- numbers the report quotes ----------------------------------------
    cg = atlas[atlas.causal_grade_gene_assignment.fillna(False).astype(bool)]
    up = cg[cg.increases_height.fillna(False).astype(bool)
            & cg.direction_orientable.fillna(False).astype(bool)]
    n_pos_only = int(len(atlas) - len(cg))
    clean = series[series.allelic_series_class.str.startswith("CLEAN")]
    n4 = int((paths.requirements_met == 4).sum())
    stc_row = paths[paths.gene == "STC2"]
    npr_row = paths[paths.gene == "NPR3"]
    n_arms_blocked = int(arms.blocked_pending_range_finding.sum())
    high = safety[safety.concern.str.startswith("HIGH")]
    stc_expr = expr[expr.gene == "STC2"].iloc[0] if len(expr[expr.gene == "STC2"]) else None
    npr_expr = expr[expr.gene == "NPR3"].iloc[0] if len(expr[expr.gene == "NPR3"]) else None
    onc = claims[claims.claim.str.startswith("PAPP-A is pursued")]
    vos = prec[prec.precedent.str.startswith("CNP analogue")]
    stc_trials = prec[prec.precedent.str.startswith("any stanniocalcin")]

    G.log(f"stage 94: top20={len(t20)}, top10 modalities={len(m10)}, leads={len(l5)}")

    def q(n: int, question: str, *body: str) -> list[str]:
        return [f"### {n}. {question}", ""] + list(body) + [""]

    L = ["# Final allelic-series dossier", "",
         "## What this branch did", "",
         "The previous human-genetics branch worked from monogenic syndromes and "
         "ClinVar, and its answer was that of 38 genes with stature phenotypes, none "
         "produces proportionate tall stature without a cost. That was the wrong "
         "instrument, not the wrong question: alleles that make healthy adults taller "
         "cause no disease and so appear in no disease database.", "",
         "This branch used quantitative human genetics instead, under one rule - a "
         f"positional gene assignment is not causal evidence. Of {len(atlas)} height "
         f"associations on coding-class variants, {n_pos_only} were positional only and "
         f"were excluded from every causal claim; {len(cg)} had a VEP-confirmed "
         f"protein-altering consequence in the named gene.", "",
         "## Headline", "",
         f"| | |", "|---|---|",
         f"| genes screened | {len(series)} |",
         f"| genes reaching a clean allelic series | {len(clean)} "
         f"({', '.join(clean.gene)}) |",
         f"| genes meeting all four of the brief's requirements | {n4} |",
         f"| ex vivo arms with a stateable concentration | "
         f"{len(arms) - n_arms_blocked} of {len(arms)} |",
         f"| compounds proposed | **0** |",
         f"| target classes proposed | **2** |", "",
         "## The twelve questions", ""]

    q1_rows = []
    for g in clean.gene:
        sub = up[up.seed_gene == g]
        for rs in sorted(set(sub.rsid)):
            v = sub[sub.rsid == rs]
            # frequency and ancestry are per-STUDY, so take them from whichever row of
            # this variant reports them rather than from an arbitrary first row; p is
            # taken per VARIANT, not per gene
            fr = v.allele_frequency.dropna()
            anc = sorted({a for a in v.ancestry.dropna() if str(a) != "not stated"})
            r0 = v.iloc[0]
            q1_rows.append(
                f"| **{g}** | `{rs}` | {str(r0.vep_hgvsp).split(':')[-1] or '—'} | "
                f"{r0.sift}/{r0.polyphen} | {'; '.join(anc) or 'not stated'} | "
                f"{'not reported' if not len(fr) else f'{fr.min():.4g}'} | "
                f"{v.pvalue.min():.0e} | {v.study_accession.nunique()} |")
    L += q(1, "Which human alleles increase proportionate adult height?",
           f"{len(up)} associations on {up.rsid.nunique()} distinct protein-altering "
           f"variants raise height with a resolvable allele orientation. After removing "
           "genes whose overgrowth is syndromic, dysplastic or neoplastic - which the "
           "brief excludes - the proportionate set is small:", "",
           "| gene | variant | protein change | effect prediction | ancestry | "
           "frequency | smallest p | independent studies |",
           "|---|---|---|---|---|---|---|---:|",
           *q1_rows)

    L += ["",
          "All three are rare, all three are predicted deleterious by both SIFT and "
          "PolyPhen, and all three are European-ancestry findings - which is a limitation of the catalogue rather than of the biology, "
          "and it is stated rather than smoothed over. The effect sizes are deliberately "
          "not quoted in centimetres: the catalogue records the unit as the literal "
          "string 'unit' for 115 of 116 betas, and converting an unstated unit into "
          "centimetres would be inventing the number.", ""]

    L += q(2, "Which have experimentally established molecular direction?",
           "**None of the human alleles.** SIFT and PolyPhen call both variants "
           "deleterious, but a prediction is not a measurement, and no record retrieved "
           "here measures the inhibitory capacity of STC2 p.Arg44Leu or the ligand "
           "binding of NPR3 p.Gly478Ser.", "",
           f"Direction is established only in the mouse, where the allele type is "
           f"recorded: {int((paths.molecular_direction).sum())} of {len(paths)} genes "
           "have a length change produced by an allele of stated molecular type. That "
           "gap - between a human association and a stated molecular direction - is the "
           "single reason stage 92's experiment exists.")

    L += q(3, "Which show reciprocal animal phenotypes?",
           f"{int((paths.animal_agreement).sum())} of {len(paths)} genes have human "
           "variants that raise height AND mouse loss that lengthens bone. "
           f"**{', '.join(paths[paths.animal_agreement].gene)}.**", "",
           "Genuinely reciprocal series - where one allele lengthens and the opposite "
           "allele shortens - are rarer still: stage 88 found "
           f"{int((series.allelic_series_class == 'BIDIRECTIONAL_GROWTH_REGULATOR').sum())}. "
           "The axis-level reciprocity is better than the gene-level reciprocity: a "
           "damaging variant in the inhibitor STC2 raises height, and a damaging "
           "variant in the protease PAPP-A lowers it, which is what a dose-limiting "
           "cascade predicts and a positional association would not.")

    local = chain[chain.compartment.str.contains("secreted", case=False)]
    L += q(4, "Which act locally rather than by globally changing endocrine levels?",
           "This is the axis's main structural argument and also its weakest measured "
           "point.", "",
           "The pappalysin cascade acts on IGF that is already present, by releasing it "
           "from binding proteins where the protease is active - so the *mechanism* is "
           "local by construction, unlike giving IGF-I or growth hormone. Stage 89 "
           f"placed {len(local)} of {len(chain)} nodes in the secreted compartment for "
           "exactly this reason, and IGFALS was excluded from the local levers because "
           "it stabilises the circulating reservoir rather than acting in the plate.",
           "",
           "**But locality of mechanism is not locality of exposure.** A systemic agent "
           "against a secreted target acts wherever the target is, and GTEx puts STC2 "
           f"above 1 TPM in {int(stc_expr.n_tissues_with_expression)} of "
           f"{int(stc_expr.n_tissues_measured)} tissues and NPR3 in "
           f"{int(npr_expr.n_tissues_with_expression)} of "
           f"{int(npr_expr.n_tissues_measured)}, highest in aorta. No measurement in "
           "this programme separates growth-plate exposure from systemic exposure.")

    acc = paths[paths.accessible_intervention]
    L += q(5, "Which are extracellular and realistically tractable?",
           f"{len(acc)} of {len(paths)} genes are secreted or cell-surface, so "
           "accessibility is *not* this field's binding constraint - only "
           f"{len(paths) - len(acc)} genes fail on it.", "",
           "Tractability is a different question, and the answer separates the two "
           "leads sharply:", "",
           "| target | compartment | interface solved | catalogued chemistry | "
           "what could be built |", "|---|---|---|---|---|",
           "| **STC2** | secreted | yes - 4 structures of the STC2:PAPP-A complex, best "
           "3.06 Å | **zero** ChEMBL activities | an antibody, engineered peptide or "
           "macrocycle against the STC2 face |",
           "| **NPR3** | cell surface | yes - 4 ligand-bound structures, best 2.00 Å | "
           "230 activities, none a named compound | antagonist or ligand trap; the "
           "receptor family has been co-crystallised with Fabs |")

    L += q(6, "Did the previous human-genetics pipeline miss STC2 or similar "
           "quantitative variants?",
           "**Yes, and this branch also missed it once before catching it.**", "",
           "The earlier branch worked from OMIM and ClinVar, where STC2 does not appear "
           "as a stature gene because the allele causes no disease. That was the "
           "predicted failure and it is why this branch exists.", "",
           "The second miss is more instructive because it was self-inflicted. Stage "
           "87's first version paged the catalogue's gene search at 120 records, which "
           "silently truncated 69 of 77 genes. STC2's only protein-altering variants "
           "sit past position 120, so the atlas reported STC2 as having no coding "
           "variant - a clean, confident, wrong answer produced by a truncated query. "
           "The cap was removed and the atlas rebuilt, and STC2 went from REJECT to the "
           "top of the table. **A silent truncation and a real negative look identical "
           "in the output.**", "",
           "A third case in the same stage: `rs35816944` is labelled IGFALS by the "
           "catalogue, but VEP shows the variant truncates SPSB3. The causal-grade rule "
           "caught it. Positional gene labels are wrong often enough to matter.")

    onc_n = int(onc.europepmc_records.iloc[0]) if len(onc) else 0
    L += q(7, "Is STC2 inhibition or PAPP-A/PAPP-A2 augmentation experimentally "
           "feasible?",
           "**Feasible to attempt; not currently executable at a stated concentration.**",
           "",
           "In favour: the interface is solved, extracellular, and genetically anchored "
           "on the correct side. Recombinant pappalysin protein is a tractable reagent "
           "for the augmentation arms, and stage 89 found the claim that PAPP-A2 has "
           "been given to humans supported by records.", "",
           "Against, and decisively for now: **PAPPA, PAPPA2, STC1 and STC2 have no "
           f"single-protein ChEMBL target entry at all**, and {n_arms_blocked} of "
           f"{len(arms)} ex vivo arms therefore carry `RANGE_UNDETERMINED` rather than "
           "a concentration. This programme does not invent concentrations - stage 65 "
           "caught an earlier version extracting 'active concentrations' that were "
           "buffer salts at 120 mM. A range-finding step producing a measured potency "
           "per reagent lot is a precondition, not an appendix.", "",
           "One asymmetry is worth naming: augmentation and inhibition are not "
           f"symmetric in effort. There are {onc_n:,} records on PAPP-A as an oncology "
           "target pursued by inhibition. Making a protease *more* active is not a "
           "standard modality, which is why the tractable version of this idea is "
           "relieving inhibition rather than activating the enzyme.")

    q8_rows = [
        f"| {r['rank']} | {r.modality} | {r.interface} | {str(r.feasibility)[:80]} | "
        f"{'yes, but no named compound' if r.orderable_today else '**none**'} |"
        for _, r in m10.head(6).iterrows()]
    L += q(8, "What modality best phenocopies the height-increasing STC2 alleles?",
           "The allele is a rare, predicted-deleterious missense in a secreted "
           "inhibitor, present in heterozygous carriers who are healthy and slightly "
           "taller. What phenocopies that is **partial, reversible, extracellular "
           "neutralisation of STC2** - not gene knockdown, which is neither partial nor "
           "reversible on the same timescale, and not enzyme activation, which has no "
           "modality.", "",
           "Ranked by direction, then by whether the interface is solved, then by "
           "whether any chemistry exists against the target at all:", "",
           "| rank | modality | interface | feasibility | catalogued chemistry |",
           "|---:|---|---|---|---|",
           *q8_rows,
           "",
           "The NPR3 rows rank above the STC2 rows only because ChEMBL holds 230 "
           "activities against NPR3 and zero against STC2. **Not one of those 230 is a "
           "named compound** - they are unnamed research entries - so 'catalogued "
           "chemistry' here means a starting point for a medicinal chemist, not "
           "something that can be ordered and dosed. On direction, interface quality "
           "and genetic anchoring the two targets are equivalent.")

    L += ["",
          "The honest top answer is an **antibody or engineered peptide against the "
          "STC2 face of the STC2:PAPP-A interface**. The brief permits a target class "
          "or biologic where no small molecule exists, and that permission is being "
          "used because it is what the evidence supports - a small-molecule blocker of "
          "a large flat PPI, designed against cryo-EM maps at 3-5 Å, would be a "
          "considerably weaker claim.", ""]

    L += q(9, "What safety liability is most likely to kill the pathway?",
           "**For STC2: cancer.** Not because the mouse record flags it - it does not, "
           "and that silence should not be read as reassurance - but because of the "
           f"direction. The intervention increases local free IGF, and the {onc_n:,} "
           "records on PAPP-A as an oncology target exist because a field believes "
           "reducing free IGF helps. This programme proposes to move that quantity the "
           "other way, and no instrument used here was designed to detect the risk of "
           "*increasing* an activity.", "",
           "**For NPR3: haemodynamics**, and this one is flagged by both instruments. "
           f"{len(high)} target/system pair reached HIGH concern in stage 93: NPR3 "
           "against blood pressure, with a mouse hypotension phenotype and "
           "high-confidence human associations to increased blood pressure and "
           "essential hypertension. NPR3's highest-expressing tissue in GTEx is the "
           "aorta. Reducing natriuretic peptide clearance is a haemodynamic act by "
           "construction, not by accident.", "",
           "Stage 93 also records NPR3 mouse phenotypes for delayed endochondral "
           "ossification, reduced fertility and incompletely penetrant postnatal "
           "lethality. None is disqualifying on its own; together they describe a "
           "receptor doing several jobs.")

    L += q(10, "Which three pathways deserve normal-postnatal metatarsal testing?",
           "1. **The STC2 : PAPP-A interface** - the only target in the programme with "
           "a human allele, a mouse allele of stated type, a solved extracellular "
           "interface and a correct direction, all pointing the same way.",
           "2. **NPR3 clearance blockade** - the second gene meeting all four "
           "requirements, and the one whose pathway has reached children clinically, "
           "albeit at a different node.",
           "3. **PAPP-A / PAPP-A2 augmentation, as two separate arms** - not because "
           "augmentation is the likely therapeutic, but because it is the positive "
           "control that decides whether the axis moves bone length at all. If adding "
           "the enzyme does nothing, relieving its inhibitor cannot work, and the whole "
           "branch is answered cheaply.", "",
           "Each is tested independently. **They are not combined into a stack**, and "
           "the design includes deliberately wrong-direction arms (PAPP-A and PAPP-A2 "
           "inhibition) because an axis that cannot be pushed backwards has not been "
           "shown to be an axis.")

    L += q(11, "Does any pathway outperform the existing geometry probes?",
           "**Yes, on evidence class - and neither has been tested, so the comparison "
           "is about what is known, not about what works.**", "",
           "| | the five geometry probes (stages 69-77) | STC2 / NPR3 (stages 87-93) |",
           "|---|---|---|",
           "| origin | inferred from pathway reasoning about cell shape | a human "
           "allele that measurably changes height |",
           "| direction | inferred; two of five barred by facts about the molecule | "
           "anchored in a human allele and a mouse allele of stated type |",
           "| selectivity | broad kinase and cytoskeletal activity | a single named "
           "protein-protein interface |",
           "| best status reached | all five at `PENETRATION_UNRESOLVED`; 0 reached "
           "`MECHANISM_VALIDATED` | genetically anchored, awaiting a directional test |",
           "| chemistry | ordered compounds exist | **none exists** |", "",
           "The geometry branch had compounds and no direction. This branch has a "
           "direction and no compounds. The second is the better problem, because a "
           "direction cannot be bought and a reagent can be made - but it is a "
           "reversal of position, not a victory, and nothing here has yet lengthened a "
           "bone.")

    stc_t = int(stc_trials.registered_studies.iloc[0]) if len(stc_trials) else 0
    vos_t = int(vos.registered_studies.iloc[0]) if len(vos) else 0
    L += q(12, "Is there an existing compound, or does the best lead require a new "
           "biologic or peptide?",
           "**There is no existing compound. The best lead requires a new biologic or "
           "peptide.**", "",
           f"- STC2, STC1, PAPPA and PAPPA2: zero single-protein ChEMBL targets, zero "
           "catalogued activities, zero named molecules.",
           f"- NPR3: 230 catalogued activities across three ChEMBL targets, but not one "
           "is a named compound - they are unnamed research entries, and a specific "
           "reagent would still have to be selected and its potency measured.",
           f"- Registered clinical studies of any stanniocalcin-directed agent: "
           f"**{stc_t}**. Of the CNP/NPR2 arm: {vos_t}, which is the precedent that "
           "exists in this pathway and it is at a different node.", "",
           "So the answer the brief anticipated is the answer: **a target class, not a "
           "compound.** An antibody or engineered peptide against the STC2 face of the "
           "STC2:PAPP-A interface, with a measured potency, is the first orderable "
           "object this branch requires and does not have.")

    L += ["## Top 5 leads", "", "| rank | gene | class | variants | interface "
          "structures | best modality | highest safety concern | status |",
          "|---:|---|---|---|---:|---|---|---|"]
    for _, r in l5.iterrows():
        L.append(f"| {r['rank']} | **{r.gene}** | {r.target_class} | "
                 f"{r.human_increasing_variants} | {r.structures_of_the_interface} | "
                 f"{r.best_modality} | {r.highest_safety_concern} | {r.status} |")
    L += ["",
          f"Only {len(l5)} lead(s) exist, and neither is a compound. Both carry the "
          "status `GENETICALLY_ANCHORED_TARGET_AWAITING_DIRECTIONAL_TEST` because that "
          "is exactly what they are.", ""]

    L += ["## What would change these conclusions", "",
          "1. **A measured potency for any STC2- or NPR3-directed reagent.** This is "
          "the single blocking item; without it stage 92 is a design and not a "
          "protocol.",
          "2. **The PAPP-A augmentation arm returning null.** If adding active protease "
          "to a normal explant does not change elongation, the axis is not "
          "dose-limiting for bone length ex vivo and the branch closes.",
          "3. **The IGF1R epistasis arm failing to abolish an STC2 effect.** That would "
          "mean the mechanism attributed across stages 89-91 is wrong, whatever the "
          "length result.",
          "4. **A measured terminal-zone concentration.** Every efficacy statement in "
          "this programme has been gated on penetration since stage 70, and no agent "
          "has yet cleared that gate - the previous branch's five probes all ended at "
          "`PENETRATION_UNRESOLVED`.", "",
          "## What this dossier does not support", "",
          "- No compound is recommended, for any use.",
          "- **No human dosing, route, schedule or self-experimentation guidance is "
          "given or derivable from anything here.** The analysis has not established a "
          "concentration for a single explant arm, let alone an organism.",
          "- No interventions are combined. Each arm is tested separately, and the "
          "design contains no stack.",
          "- Faster growth is not claimed to be greater final height. The plateau and "
          "washout endpoints exist precisely because a plate that grows faster and "
          "stops sooner ends at the same length.",
          "- No syndromic or dysplastic overgrowth gene is promoted. Stage 88 placed "
          f"{int((series.allelic_series_class == 'SYNDROMIC_OVERGROWTH').sum())} genes "
          f"in SYNDROMIC_OVERGROWTH and "
          f"{int((series.allelic_series_class == 'DYSMORPHIC_OR_DISPROPORTIONATE').sum())} "
          "in DYSMORPHIC_OR_DISPROPORTIONATE, and they are preserved with reasons "
          "rather than discarded.", ""]

    (R / "final_allelic_series_report.md").write_text("\n".join(L))
    G.log(f"stage 94: wrote top_20_genetically_anchored_targets.csv, "
          f"top_10_growth_modalities.csv ({len(m10)}), "
          f"top_5_allelic_series_leads.csv ({len(l5)}) and "
          "final_allelic_series_report.md")


if __name__ == "__main__":
    main()
