"""Frozen configuration: gene sets and sample table.

Gene sets are fixed in docs/PREREGISTRATION.md and must not be edited after
the pre-registration commit. Any change belongs in docs/DEVIATIONS.md.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
for _p in (PROC, TAB, FIG):
    _p.mkdir(parents=True, exist_ok=True)

SEED = 0

# --- Gene sets (pre-registered section 1) -----------------------------------

ZONE_MARKERS = {
    "resting": ["SFRP5", "APOE", "CLU", "PTCH1", "CYTL1", "GAS1", "PTHLH", "RAMP3"],
    "proliferative": ["CCND1", "MKI67", "E2F1"],
    "prehypertrophic": ["IHH", "ALPL"],
    "hypertrophic": ["COL10A1", "MMP13"],
}

# Full protocol repression-target set.
TARGETS_FULL = ["COL10A1", "MMP13", "IBSP", "SPP1", "VEGFA", "ALPL", "PANX3", "SP7"]

# Genes that are simultaneously zone markers and repression targets. Per the
# pre-registration these are used for zone assignment, NOT for the primary score.
TARGET_ZONE_OVERLAP = ["COL10A1", "MMP13", "ALPL"]

# PRIMARY score set: disjoint from the zone markers.
TARGETS_PRIMARY = [g for g in TARGETS_FULL if g not in TARGET_ZONE_OVERLAP]

PATHWAY_COMPONENTS = [
    "HDAC4", "HDAC5", "HDAC7", "SIK1", "SIK2", "SIK3", "CAMK4", "CAMK2D",
    "PPP2CA", "PPP2R1A", "PTH1R", "PTHLH", "MEF2C", "RUNX2", "SMAD2", "SMAD3",
    "HIF1A", "EPAS1",
]

HOUSEKEEPING = ["ACTB", "GAPDH", "RPL13A", "TBP", "PPIA", "B2M", "HPRT1", "SDHA"]

# Dissociation / stress signature regressed out during QC.
STRESS_GENES = [
    "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "HSPA1A", "HSPA1B",
    "HSPB1", "HSP90AA1", "DNAJB1", "SOCS3", "ZFP36", "IER2", "DUSP1",
]

ZONE_ORDER = ["resting", "proliferative", "prehypertrophic", "hypertrophic"]

# --- Samples ----------------------------------------------------------------
# Age columns are deliberately None: GEO deposits no per-sample age for any of
# these libraries. See docs/DATA_AUDIT.md. Do not fill these in by guessing.

HUMAN_SAMPLES = [
    # gsm, file, donor, arm
    ("GSM9328218", "GSM9328218_P30453_1001.h5", "P30453", "primary"),
    ("GSM9328219", "GSM9328219_P30453_1002.h5", "P30453", "vehicle"),
    ("GSM9328220", "GSM9328220_P30453_1003.h5", "P30453", "gh"),
    ("GSM9328221", "GSM9328221_P31011_1001.h5", "P31011", "primary"),
    ("GSM9328222", "GSM9328222_P31011_1002.h5", "P31011", "vehicle"),
    ("GSM9328223", "GSM9328223_P31011_1003.h5", "P31011", "gh"),
    ("GSM9328224", "GSM9328224_P25452_001.h5", "P25452", "primary"),
    ("GSM9328225", "GSM9328225_P25452_004.h5", "P25452", "vehicle"),
    ("GSM9328226", "GSM9328226_P25452_005.h5", "P25452", "vehicle"),
    ("GSM9328227", "GSM9328227_P25452_007.h5", "P25452", "gh"),
    ("GSM9328228", "GSM9328228_P25452_008.h5", "P25452", "gh"),
    ("GSM9328229", "GSM9328229_P22202_1015.h5", "P22202", "primary"),
]

MOUSE_SAMPLES = [
    ("GSM9328230", "GSE288028/GSM9328230_P27153_1001.h5", "GSE288028", "M_P27153_r1"),
    ("GSM9328231", "GSE288028/GSM9328231_P27153_1002.h5", "GSE288028", "M_P27153_r2"),
    ("GSM8769462", "GSE288529/GSM8769462_3402_filtered_feature_bc_matrix.h5",
     "GSE288529", "M_4wk_female"),
]

# QC thresholds
QC = dict(min_genes=500, max_genes=7500, max_pct_mt=10.0, min_cells_per_gene=3)
MIN_CELLS_PER_ZONE = 100


def to_mouse(genes):
    """Human symbol -> mouse symbol (Title case). Valid for these orthologue-1:1 sets."""
    return [g.capitalize() for g in genes]
