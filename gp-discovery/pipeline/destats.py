"""
Differential-expression statistics.

No R is available in this environment, so limma's moderated t-test is
reimplemented here (empirical-Bayes variance shrinkage against a fitted prior),
which is the correct tool for the log-CPM / normalised-intensity matrices that
several of these series ship. Raw-count series are handled with PyDESeq2
elsewhere. Both routes report effect sizes, not just p-values (requirement E).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import special, stats


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def _trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """Solve trigamma(y) = x for y (Newton iteration, as in limma)."""
    x = np.asarray(x, dtype=float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.max(np.abs(dif / np.maximum(y, 1e-12))) < 1e-8:
            break
    return y


def fit_prior(s2: np.ndarray, df: float):
    """Estimate prior variance s0^2 and prior df (limma's fitFDist)."""
    s2 = s2[np.isfinite(s2) & (s2 > 0)]
    if s2.size < 10:
        return np.nan, 0.0
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = e.var(ddof=1) - special.polygamma(1, df / 2)
    if evar <= 0:  # no excess variance -> infinite prior df (fully shrunk)
        return float(np.exp(emean)), np.inf
    df_prior = 2 * _trigamma_inverse(np.array([evar]))[0]
    s2_prior = float(np.exp(emean + special.digamma(df_prior / 2) - np.log(df_prior / 2)))
    return s2_prior, float(df_prior)


def moderated_ttest(mat: pd.DataFrame, group: pd.Series, ref: str, alt: str) -> pd.DataFrame:
    """
    Two-group moderated t-test on a log-scale matrix (genes x samples).

    Returns log2FC (alt - ref), moderated t, p, BH-FDR, group means, and
    Cohen's d, preserving biological replicates as columns.
    """
    g = pd.Series(group).reindex(mat.columns)
    a = mat.loc[:, (g == alt).values].astype(float)
    b = mat.loc[:, (g == ref).values].astype(float)
    na, nb = a.shape[1], b.shape[1]
    if na < 2 or nb < 2:
        raise ValueError(f"need >=2 replicates per group, got {alt}={na} {ref}={nb}")

    ma, mb = a.mean(axis=1), b.mean(axis=1)
    lfc = ma - mb
    df_res = na + nb - 2
    ss = ((a.sub(ma, axis=0) ** 2).sum(axis=1) + (b.sub(mb, axis=0) ** 2).sum(axis=1))
    s2 = ss / df_res

    s2_prior, df_prior = fit_prior(s2.values, df_res)
    if not np.isfinite(s2_prior):
        s2_post, df_total = s2.values, np.full(s2.size, df_res)
    elif np.isinf(df_prior):
        s2_post = np.full(s2.size, s2_prior)
        df_total = np.full(s2.size, np.inf)
    else:
        s2_post = (df_res * s2.values + df_prior * s2_prior) / (df_res + df_prior)
        df_total = np.full(s2.size, df_res + df_prior)

    se = np.sqrt(s2_post * (1.0 / na + 1.0 / nb))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = lfc.values / se
        p = 2 * stats.t.sf(np.abs(t), df_total)
        d = lfc.values / np.sqrt(np.where(s2_post > 0, s2_post, np.nan))

    return pd.DataFrame(
        {
            "log2FC": lfc.values,
            "mean_" + alt: ma.values,
            "mean_" + ref: mb.values,
            "t_mod": t,
            "pvalue": p,
            "FDR": bh_fdr(p),
            "cohens_d": d,
            "n_alt": na,
            "n_ref": nb,
        },
        index=mat.index,
    )


def cpm_to_log(df: pd.DataFrame, prior: float = 1.0) -> pd.DataFrame:
    return np.log2(df.astype(float) + prior)


def collapse_duplicate_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate gene rows by max mean expression (array probes)."""
    if df.index.duplicated().any():
        df = df.assign(_m=df.mean(axis=1)).sort_values("_m", ascending=False)
        df = df[~df.index.duplicated(keep="first")].drop(columns="_m")
    return df
