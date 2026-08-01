"""
Stage 54 - active-learning expansion from PILOT_96 to EXPANSION_384.

Working implementation. Morgan fingerprints from RDKit, a random-forest surrogate
with out-of-bag uncertainty per objective, and a multi-objective acquisition that
rewards durable elongation, preserved proliferation, preserved survival,
preserved matrix, mechanistic novelty and uncertainty - and penalises predicted
cytotoxic chemical space.

It is explicitly NOT a maximiser of predicted length. A model trained on 96
compounds that picks the 384 highest predicted length gains would pick 384
neighbours of whatever the pilot's best compound happened to be, and this project
has already shown three times what happens when a search optimises a proxy.

The model is exercised on simulated pilot outcomes so the selection machinery can
be validated end to end. Real pilot results replace `simulate_pilot_outcomes()`;
nothing else changes.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.ensemble import RandomForestRegressor  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

R = G.RESULTS
FIG = R / "figures"
OUT = R / "stage54"
OUT.mkdir(parents=True, exist_ok=True)
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#dcdbd6"
S1, S2, S3, S8 = "#2a78d6", "#eb6834", "#1baf7a", "#e34948"
AMBER = "#d99a12"

OBJECTIVES = {
    # name: (direction, weight in the acquisition)
    "durable_length_gain": (+1, 0.34),
    "edu_preserved": (+1, 0.16),
    "viability_preserved": (+1, 0.16),
    "matrix_preserved": (+1, 0.14),
}
W_NOVELTY = 0.12
W_UNCERTAINTY = 0.08
W_CYTOTOX_PENALTY = 0.30       # subtracted, not weighted into the sum

FEATURES = [
    ("morgan_fp_2048", "structure", "RDKit Morgan radius 2, 2048 bits, folded"),
    ("n_annotated_targets", "target", "promiscuity; a proxy for polypharmacology risk"),
    ("mechanism_family_onehot", "mechanism", "the 15 stage-49 families"),
    ("primary_target_hashed", "target", "hashed target identity, 64 buckets"),
    ("clinical_phase_ordinal", "development", "Launched=0 ... Preclinical=4"),
    ("log_potency_nM", "potency", "Guide to Pharmacology affinity where retrievable"),
    ("cartilage_literature_count", "prior", "Europe PMC records, gene x cartilage"),
    ("observed_length_effect_mm", "outcome", "PILOT only - the label"),
    ("observed_edu_delta", "outcome", "PILOT only - the label"),
    ("observed_tunel_delta", "outcome", "PILOT only - the label"),
    ("observed_terminal_cell_volume_delta", "outcome", "PILOT only - the label"),
    ("observed_matrix_delta", "outcome", "PILOT only - the label"),
    ("observed_washout_plateau_delta", "outcome", "PILOT only - the label"),
    ("observed_toxicity_flag", "outcome", "PILOT only - the label"),
    ("assay_confidence", "quality", "mean stage-51 measurement confidence for that compound"),
]


def fingerprints(smiles):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    X, ok = [], []
    for s in smiles:
        m = Chem.MolFromSmiles(s) if isinstance(s, str) and s else None
        if m is None:
            X.append(np.zeros(2048, dtype=np.uint8))
            ok.append(False)
        else:
            arr = np.zeros(2048, dtype=np.uint8)
            from rdkit import DataStructs
            DataStructs.ConvertToNumpyArray(gen.GetFingerprint(m), arr)
            X.append(arr)
            ok.append(True)
    return np.array(X), np.array(ok)


def build_matrix(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    fp, ok = fingerprints(df.smiles)
    fam = pd.get_dummies(df.family_primary.fillna("other")).to_numpy(float)
    phase = df.clinical_phase.map({"Launched": 0, "Phase 3": 1, "Phase 2": 2,
                                   "Phase 1": 3, "Preclinical": 4}).fillna(4).to_numpy(float)
    ntg = pd.to_numeric(df.n_targets, errors="coerce").fillna(1).to_numpy(float)
    lit = np.log1p(pd.to_numeric(df.cartilage_bone_records, errors="coerce").fillna(0))
    tgt = np.array([[hash(str(t)) % 64 == b for b in range(64)]
                    for t in df.primary_target], float)
    X = np.hstack([fp, fam, phase[:, None], ntg[:, None], lit.to_numpy()[:, None], tgt])
    return X, ok


def simulate_pilot_outcomes(pilot: pd.DataFrame, rng) -> pd.DataFrame:
    """Stand-in for real pilot results. Replace this function, nothing else.

    The simulation gives length effects a weak structural basis so the surrogate
    has something learnable, plus a cytotoxicity axis that is correlated with
    promiscuity - which is what makes the penalty term do work."""
    n = len(pilot)
    fp, _ = fingerprints(pilot.smiles)
    w = rng.normal(0, 1, fp.shape[1])
    struct = fp @ w
    struct = (struct - struct.mean()) / (struct.std() + 1e-9)
    promisc = pd.to_numeric(pilot.n_targets, errors="coerce").fillna(1).to_numpy(float)
    tox = np.clip((promisc - promisc.mean()) / (promisc.std() + 1e-9)
                  + rng.normal(0, 0.6, n), -3, 3)
    length = 0.05 * struct + rng.normal(0, 0.06, n)
    return pd.DataFrame({
        "pert_iname": pilot.pert_iname,
        "observed_length_effect_mm": length,
        "observed_edu_delta": -0.12 * np.clip(tox, 0, None) + rng.normal(0, 0.03, n),
        "observed_tunel_delta": 0.10 * np.clip(tox, 0, None) + rng.normal(0, 0.02, n),
        "observed_matrix_delta": -0.08 * np.clip(tox, 0, None) + rng.normal(0, 0.04, n),
        "observed_terminal_cell_volume_delta": 0.6 * length + rng.normal(0, 0.05, n),
        "observed_washout_plateau_delta": 0.7 * length - 0.15 * np.clip(tox, 0, None)
        + rng.normal(0, 0.05, n),
        "observed_toxicity_flag": (tox > 1.2).astype(float),
        "assay_confidence": np.clip(rng.normal(0.86, 0.05, n), 0, 1),
    })


def fit_surrogates(Xtr, ytr_dict, rng_seed=54):
    models = {}
    for k, y in ytr_dict.items():
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2, oob_score=True,
                                  random_state=rng_seed, n_jobs=-1)
        m.fit(Xtr, y)
        models[k] = m
    return models


def predict_with_uncertainty(model, X):
    per_tree = np.stack([t.predict(X) for t in model.estimators_])
    return per_tree.mean(0), per_tree.std(0)


def acquisition(pred: dict, unc: dict, novelty: np.ndarray, tox: np.ndarray) -> np.ndarray:
    """Multi-objective: reward all four phenotype objectives, novelty and uncertainty;
    subtract predicted cytotoxic space. Objectives are z-scored so no single one
    dominates by unit scale."""
    def z(v):
        s = v.std()
        return (v - v.mean()) / (s if s > 1e-9 else 1.0)
    score = np.zeros(len(novelty))
    for k, (direction, w) in OBJECTIVES.items():
        score = score + w * direction * z(pred[k])
    score = score + W_UNCERTAINTY * z(np.mean([unc[k] for k in OBJECTIVES], axis=0))
    score = score + W_NOVELTY * z(novelty)
    score = score - W_CYTOTOX_PENALTY * z(tox)
    return score


def select(cand: pd.DataFrame, scores: np.ndarray, n: int,
           per_family_cap: float = 0.18) -> pd.DataFrame:
    """Greedy by score with a per-family cap, so the acquisition cannot collapse
    onto one mechanism even if it wants to."""
    cand = cand.copy()
    cand["acq"] = scores
    cap = max(2, int(n * per_family_cap))
    sel, per_fam, per_target = [], {}, set()
    for _, r in cand.sort_values("acq", ascending=False).iterrows():
        if len(sel) >= n:
            break
        f = r.family_primary
        if per_fam.get(f, 0) >= cap:
            continue
        if r.primary_target in per_target and per_fam.get(f, 0) >= cap // 2:
            continue
        sel.append(r.pert_iname)
        per_fam[f] = per_fam.get(f, 0) + 1
        per_target.add(r.primary_target)
    return cand[cand.pert_iname.isin(sel)]


def figure37(cand, sel, scores, pilot_out) -> None:
    fig = plt.figure(figsize=(15.0, 7.4))
    gs = fig.add_gridspec(1, 3, wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(cand.pred_len, cand.unc_len, s=18, color="#c9ced4", alpha=0.7,
               edgecolor="none", label="not selected", zorder=2)
    m = cand.pert_iname.isin(sel.pert_iname)
    ax.scatter(cand.pred_len[m], cand.unc_len[m], s=26, color=S1, alpha=0.9,
               edgecolor=SURFACE, linewidth=0.6, label="selected", zorder=3)
    ax.set_xlabel("predicted durable length gain (mm)", color=INK2)
    ax.set_ylabel("model uncertainty (SD across trees)", color=INK2)
    ax.legend(fontsize=8.4, frameon=False)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("A  Selection is not the top of the prediction", loc="left", color=INK,
                 fontsize=11.0, pad=10)

    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(cand.pred_len, cand.pred_tox, s=18, color="#c9ced4", alpha=0.7,
               edgecolor="none", zorder=2)
    ax.scatter(cand.pred_len[m], cand.pred_tox[m], s=26, color=S3, alpha=0.9,
               edgecolor=SURFACE, linewidth=0.6, zorder=3)
    ax.set_xlabel("predicted durable length gain (mm)", color=INK2)
    ax.set_ylabel("predicted cytotoxicity", color=INK2)
    ax.grid(True, alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("B  Cytotoxic space is avoided, not ranked", loc="left", color=INK,
                 fontsize=11.0, pad=10)

    ax = fig.add_subplot(gs[0, 2])
    fams = cand.family_primary.value_counts().index[:12]
    y = np.arange(len(fams))[::-1]
    ax.barh(y, [int((cand.family_primary == f).sum()) for f in fams], 0.62,
            color="#cddef6", edgecolor=SURFACE, linewidth=1.0, label="candidate pool")
    ax.barh(y, [int((sel.family_primary == f).sum()) for f in fams], 0.34,
            color=S1, edgecolor=SURFACE, linewidth=1.0, label="selected")
    ax.set_yticks(y)
    ax.set_yticklabels([f.replace(" / ", "/\n").replace(" (", "\n(") for f in fams],
                       fontsize=7.8)
    ax.set_xlabel("compounds", color=INK2)
    ax.legend(fontsize=8.2, frameon=False, loc="lower right")
    ax.grid(True, axis="x", alpha=0.5, linewidth=0.6); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("C  Per-family cap keeps the spread", loc="left", color=INK,
                 fontsize=11.0, pad=10)

    fig.suptitle("Active-learning expansion: PILOT_96 to EXPANSION_384",
                 x=0.006, y=0.985, ha="left", fontsize=13.8, fontweight="bold", color=INK)
    fig.text(0.006, 0.933,
             "Surrogate trained on SIMULATED pilot outcomes - the machinery is validated, the "
             "selection is not a prediction about real compounds.",
             fontsize=9.2, color=INK2, ha="left", va="top")
    fig.subplots_adjust(top=0.835, bottom=0.135, left=0.06, right=0.985)
    fig.savefig(FIG / "37_active_learning_selection.png", facecolor=SURFACE, dpi=170)
    plt.close(fig)


def main() -> None:
    rng = np.random.default_rng(54)
    full = pd.read_csv(R / "full_screen_compound_catalog.csv")
    pilot = pd.read_csv(R / "pilot_96_compound_library.csv")

    sch = pd.DataFrame(FEATURES, columns=["feature", "block", "definition"])
    sch["available_before_pilot"] = sch.block != "outcome"
    sch.to_csv(R / "active_learning_feature_schema.csv", index=False)

    out = simulate_pilot_outcomes(pilot, rng)
    tr = pilot.merge(out, on="pert_iname")
    Xtr, _ = build_matrix(tr)
    cand = full[~full.pert_iname.isin(pilot.pert_iname)].copy()
    Xc, _ = build_matrix(cand)
    # align one-hot family columns between train and candidate matrices
    k = min(Xtr.shape[1], Xc.shape[1])
    Xtr, Xc = Xtr[:, :k], Xc[:, :k]

    targets = {
        "durable_length_gain": tr.observed_washout_plateau_delta.to_numpy(),
        "edu_preserved": tr.observed_edu_delta.to_numpy(),
        "viability_preserved": (-tr.observed_tunel_delta).to_numpy(),
        "matrix_preserved": tr.observed_matrix_delta.to_numpy(),
        "cytotoxicity": tr.observed_toxicity_flag.to_numpy(),
    }
    models = fit_surrogates(Xtr, targets)
    pred, unc = {}, {}
    for kk, m in models.items():
        pred[kk], unc[kk] = predict_with_uncertainty(m, Xc)

    # novelty: distance of the candidate's target and family from what the pilot covered
    seen_t = set(pilot.primary_target.astype(str))
    seen_f = set(pilot.family_primary.astype(str))
    novelty = np.array([1.0 * (str(t) not in seen_t) + 0.5 * (str(f) not in seen_f)
                        for t, f in zip(cand.primary_target, cand.family_primary)])

    scores = acquisition(pred, unc, novelty, pred["cytotoxicity"])
    cand["pred_len"] = pred["durable_length_gain"]
    cand["unc_len"] = unc["durable_length_gain"]
    cand["pred_tox"] = pred["cytotoxicity"]
    cand["novelty"] = novelty
    cand["acquisition"] = scores
    sel = select(cand, scores, 384 - len(pilot))
    sel.to_csv(OUT / "expansion_selection_simulated.csv", index=False)
    cand.to_csv(OUT / "candidate_scores_simulated.csv", index=False)
    figure37(cand, sel, scores, out)

    oob = {k: round(float(m.oob_score_), 3) for k, m in models.items()}
    top_by_pred = set(cand.nlargest(len(sel), "pred_len").pert_iname)
    overlap = len(top_by_pred & set(sel.pert_iname)) / max(len(sel), 1)

    L = ["# Expansion selection plan", "",
         "## What this model is for", "",
         "It chooses the 384-compound expansion from the pilot's results. It is **not** a model "
         "of which compounds increase bone length; 96 observations cannot support that. It is a "
         "device for spending the next 288 wells better than at random.", "",
         "## Why it does not maximise predicted length", "",
         f"If the expansion were chosen by predicted durable length gain alone, "
         f"**{overlap:.0%}** of the selection would be the same compounds. The rest of the "
         "difference is the point. A pure maximiser trained on 96 compounds selects neighbours of "
         "whatever the pilot's best compound happened to be - the same failure mode as ranking "
         "genes by a connectivity score (stages 15-22) or by a phenotype-first literature score "
         "(stages 23-35). Both produced a lead that later collapsed.", "",
         "## Acquisition function", "",
         "```",
         "score = Σ wᵢ · z(predicted objectiveᵢ)          # four phenotype objectives",
         "      + w_unc     · z(mean model uncertainty)   # explore where the model is unsure",
         "      + w_novelty · z(mechanistic novelty)      # unseen target or family",
         "      − w_tox     · z(predicted cytotoxicity)   # avoid, do not merely down-rank",
         "```", "",
         "| term | weight | direction |", "|---|---:|---|"]
    for k, (dirn, w) in OBJECTIVES.items():
        L.append(f"| {k} | {w:.2f} | {'maximise' if dirn > 0 else 'minimise'} |")
    L += [f"| model uncertainty | {W_UNCERTAINTY:.2f} | maximise |",
          f"| mechanistic novelty | {W_NOVELTY:.2f} | maximise |",
          f"| predicted cytotoxicity | {W_CYTOTOX_PENALTY:.2f} | **subtract** |", "",
          "Every term is z-scored across the candidate pool before weighting, so a term does not "
          "dominate because its units are larger. Cytotoxicity is subtracted with the largest "
          "single weight in the function: a compound predicted to be cytotoxic is pushed out of "
          "the selection even if its predicted length gain is the highest in the pool.", "",
          "## Diversity constraint on top of the acquisition", "",
          "Greedy selection by score, with a hard cap of "
          f"{int((384 - len(pilot)) * 0.18)} compounds per mechanism family and a preference "
          "against repeating a primary target within a family. The acquisition can want to "
          "collapse onto one mechanism; the cap does not let it. This is a constraint, not a "
          "tie-break, and it costs predicted performance on purpose.", "",
          "## Surrogate", "",
          "Random forest, 400 trees, minimum leaf 2, one model per objective. Uncertainty is the "
          "standard deviation across trees, which is crude but honest for 96 training points - a "
          "Gaussian process with a learned kernel would report tighter intervals it has not "
          "earned.", "", "| objective | out-of-bag R² |", "|---|---:|"]
    for k, v in oob.items():
        L.append(f"| {k} | {v} |")
    L += ["",
          "**These R² values are on simulated pilot data and mean nothing about real compounds.** "
          "They are here to show the code fits and reports honestly. The out-of-bag score on real "
          "pilot data is itself a gate: if it is not clearly positive, the model is not "
          "informative and the expansion should be selected by mechanistic diversity alone. Stage "
          "56 makes that an explicit decision point.", "",
          "## Features", "", "| feature | block | definition | available before the pilot |",
          "|---|---|---|---|"]
    for _, r in sch.iterrows():
        L.append(f"| {r.feature} | {r.block} | {r.definition} | "
                 f"{'yes' if r.available_before_pilot else 'no - it is the label'} |")
    L += ["",
          "## What replaces the simulation", "",
          "`simulate_pilot_outcomes()` is the only function that needs replacing. It returns one "
          "row per pilot compound with the seven observed outcome columns and the assay "
          "confidence. Feed it the real stage-52 hit-call table and the rest of the pipeline is "
          "unchanged.", "",
          "## Failure modes this model has", "",
          "- **96 compounds is a small training set for 2,048-bit fingerprints.** The forest will "
          "lean on the low-dimensional blocks (family, phase, promiscuity) more than on "
          "structure. That is acceptable for a diversity-driven acquisition and would not be for "
          "a potency prediction.",
          "- **The novelty term rewards unseen targets, which are unseen partly because they are "
          "poorly annotated.** A compound with a blank target field looks novel. The stage-49 "
          "catalogue keeps `n_targets` so this can be diagnosed, but it is not fully solved.",
          "- **Nothing here knows about durability except through the pilot's washout data**, "
          "and the pilot's washout arm is only run on Tier-1 hits (stage 50). For most pilot "
          "compounds the durability label is missing, so `durable_length_gain` is trained on a "
          "biased subset. That bias is real and should be reported alongside any expansion "
          "selection.", ""]
    (R / "expansion_selection_plan.md").write_text("\n".join(L))
    G.log(f"active learning: pool {len(cand)}, selected {len(sel)}, "
          f"overlap with pure-max ranking {overlap:.0%}; OOB {oob}")


if __name__ == "__main__":
    main()
