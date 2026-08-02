"""Shared vocabulary and helpers for the geometry-first stages (61-68)."""
from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import spatiallib as S  # noqa: E402

# ---------------------------------------------------------------------------
# the distinction the whole strategy rests on
# ---------------------------------------------------------------------------
AXIAL_GEOMETRY = [
    r"cell height", r"height of (?:the )?(?:terminal |hypertrophic )?chondrocyte",
    r"axial (?:height|dimension|length) of", r"height[- ]to[- ]width", r"aspect ratio",
    r"cell (?:length|elongation) along", r"longitudinal (?:cell )?(?:height|diameter)",
    r"cell shape index", r"anisotrop\w+", r"oblate", r"prolate",
    r"long[- ]axis orientation", r"angular deviation", r"cell height/width",
]
VOLUME_ONLY = [r"cell volume", r"cell size", r"cell enlarge\w+", r"hypertroph\w+ volume",
               r"cell diameter", r"cell area", r"cross-?sectional area", r"2D area"]
SWELLING = [r"swell\w+", r"osmotic", r"hypotonic", r"hypertonic", r"regulatory volume",
            r"water influx", r"turgor", r"round(?:ed|ing) up", r"cell rounding"]
APPOSITIONAL = [r"appositional", r"bone width", r"diaphyseal (?:width|diameter)",
                r"radial growth", r"cortical thickness", r"transverse growth"]
COLUMN = [r"column\w* (?:organi[sz]ation|alignment|formation|orientation)",
          r"columnar (?:organi[sz]ation|arrangement|stack)", r"cell stack\w*",
          r"rotation", r"planar cell polarity", r"\bPCP\b", r"clonal column"]

METHODS_3D = [r"confocal", r"z-?stack", r"two-?photon", r"light[- ]sheet",
              r"3D (?:reconstruction|segmentation|morphometr)", r"stereolog\w+",
              r"optical section"]

# organ-culture / tissue context
POSTNATAL = [r"\bP\d{1,3}\b", r"postnatal", r"\d{1,2}[- ]day[- ]old", r"\d{1,2}[- ]week[- ]old",
             r"juvenile", r"weanling"]
EMBRYONIC = [r"\bE\d{1,2}\.?\d?\b", r"embryonic day", r"\bdpc\b", r"fetal", r"foetal"]

# ---------------------------------------------------------------------------
# stage-62 target families, verbatim from the brief
# ---------------------------------------------------------------------------
TARGET_FAMILIES = {
    "A cortical tension": ["ROCK1", "ROCK2", "MYLK", "MYH9", "MYH10", "LIMK1", "LIMK2",
                           "CFL1", "RHOA", "RAC1", "CDC42"],
    "B adhesion and matrix coupling": ["PTK2", "PXN", "TLN1", "VCL", "ITGB1", "ITGA10",
                                       "ITGA5", "CDH2", "ILK", "SRC", "PTK2B"],
    "C planar polarity and column orientation": ["VANGL1", "VANGL2", "CELSR1", "CELSR2",
                                                 "CELSR3", "FZD3", "FZD6", "DVL1", "DVL2",
                                                 "DVL3", "PRICKLE1", "PRICKLE2", "DAAM1",
                                                 "DAAM2", "IFT88", "KIF3A", "PKD1"],
    "D microtubule and centrosome organization": ["TUBB", "TUBA1A", "MAPT", "MAP4", "CLASP1",
                                                  "CLASP2", "CAMSAP3", "PCNT", "CEP164",
                                                  "PARD3", "PARD6A", "PRKCI"],
    "E ion and water mechanics": ["SLC12A2", "SLC9A1", "TRPV4", "PIEZO1", "PIEZO2", "AQP1",
                                  "AQP3", "AQP4", "AQP9", "KCNK2", "ANO1", "CLCN3",
                                  "SLC6A6", "SLC5A3"],
    "F lipid and nuclear signaling": ["RORA", "HMGCR", "SREBF2", "NPC1", "SOAT1", "DHCR7",
                                      "LSS", "SCAP", "INSIG1"],
}
FAMILY_OF = {g: f for f, gs in TARGET_FAMILIES.items() for g in gs}

# compound families to search for in stage 63
COMPOUND_QUERIES = {
    "ROCK1/2 inhibitor": ["ROCK1", "ROCK2"],
    "LIMK inhibitor": ["LIMK1", "LIMK2"],
    "myosin-II modulator": ["MYH9", "MYH10", "MYLK"],
    "Rho/Rac/Cdc42 modulator": ["RHOA", "RAC1", "CDC42"],
    "FAK / adhesion turnover": ["PTK2", "PTK2B", "SRC", "ILK"],
    "integrin-directed": ["ITGB1", "ITGA5", "ITGA10", "ITGAV"],
    "cadherin / junction": ["CDH2", "CDH1"],
    "polarity / cilia": ["PRKCI", "PKD1", "IFT88", "SMO"],
    "microtubule regulator (non-mitotic)": ["MAPT", "CLASP1", "CAMSAP3"],
    "ion / volume regulation": ["SLC12A2", "SLC9A1", "TRPV4", "PIEZO1", "PIEZO2", "AQP1",
                                "AQP4", "ANO1", "KCNK2"],
    "RORalpha / lipid pathway": ["RORA", "HMGCR", "NPC1", "SOAT1", "SREBF2"],
}

# hard rejects as intervention candidates - the brief's list
BROAD_POISONS = {
    "cytochalasin-family actin depolymerizer": [
        r"cytochalasin", r"latrunculin", r"mycalolide", r"swinholide"],
    "jasplakinolide-family actin stabilizer": [r"jasplakinolide", r"phalloidin",
                                               r"chondramide", r"dolastatin 11"],
    "microtubule poison": [r"colchicine", r"nocodazole", r"vinblastine", r"vincristine",
                           r"vinorelbine", r"paclitaxel", r"docetaxel", r"taxane",
                           r"combretastatin", r"podophyllotox", r"epothilone",
                           r"tubulin polymerization"],
    "broad myosin poison": [r"blebbistatin", r"\bBDM\b", r"butanedione monoxime"],
}


# ---------------------------------------------------------------------------
# Named reference compounds the ChEMBL target route cannot reach.
#
# The stage-63 universe is built from ChEMBL activity records against the
# stage-62 target symbols. Compounds whose target is the polymer itself (actin,
# tubulin) or a machine with no gene symbol in the family list (the V-ATPase) have
# no such record, so they are absent from the universe even though the brief names
# them. Injecting them explicitly keeps them visible and rejectable rather than
# silently missing. Every one carries the reason it is here.
# ---------------------------------------------------------------------------
REFERENCE_COMPOUNDS = [
    # name, compound_family, role, why it must appear
    ("cytochalasin D", "actin depolymeriser", "DISORGANIZATION_CONTROL",
     "largest length gain in the anchor paper and the clearest appositional widening; the "
     "brief hard-rejects it as an intervention candidate"),
    ("latrunculin B", "actin depolymeriser", "DISORGANIZATION_CONTROL",
     "second actin depolymeriser with a different binding site, so the control is not "
     "single-compound"),
    ("jasplakinolide", "actin stabiliser", "DISORGANIZATION_CONTROL",
     "opposite-direction actin perturbation with the same phenotype; hard-rejected as a "
     "candidate"),
    ("nocodazole", "microtubule poison", "REJECT",
     "microtubule depolymeriser; blocks mitosis, so it cannot preserve proliferation"),
    ("blebbistatin", "broad myosin poison", "REJECT",
     "non-selective myosin-II ATPase inhibitor; the brief's named broad poison"),
    ("bafilomycin A1", "V-ATPase inhibitor", "MECHANISTIC_PROBE",
     "the brief requires it in the panel; stage 29-35 established it is NOT an established "
     "growth intervention, so it enters as a probe of a prior claim, not a candidate"),
    ("IGF1", "growth factor", "POSITIVE_GEOMETRY_CONTROL",
     "the brief's positive control for longitudinal growth; acts on proliferation and "
     "hypertrophy, not specifically on axial shape"),
    ("mannitol", "osmotic agent", "SWELLING_CONTROL",
     "hyperosmotic agent; the negative direction of the swelling axis"),
    ("hypotonic medium", "osmotic agent", "SWELLING_CONTROL",
     "the brief's osmotic swelling control; changes volume without changing shape "
     "programme, which is exactly the confound to be excluded"),
]


def any_(pats, text) -> bool:
    return S._any(pats, str(text or ""))


def which(d, text) -> list:
    return S._which(d, str(text or ""))


def evidence_class(text: str) -> str:
    """The three classes the brief demands be kept separate."""
    t = str(text or "")
    if any_(AXIAL_GEOMETRY, t):
        return "1 direct measured axial geometry"
    if any_(VOLUME_ONLY + SWELLING + COLUMN, t) or re.search(r"histolog|morpholog|section",
                                                             t, re.I):
        return "2 general morphology without axial measurement"
    return "3 inferred mechanics without morphology data"
