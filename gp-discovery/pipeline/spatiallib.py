"""
Shared machinery for the spatial-first stages (41-47).

The rule that shapes everything here: a computational zone label is not spatial
evidence. Only an image of intact tissue containing recognisable growth-plate
architecture counts, and the pipeline has to be able to point at the figure.

So the unit of evidence is a *figure caption in an open-access full text* that
names the gene, names a spatial method, and names a growth-plate compartment -
plus the surrounding Methods text that says which probe or antibody was used and
whether it was validated.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import gputil as G  # noqa: E402

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1"
HPA = "https://www.proteinatlas.org/api/search_download.php"
MGI = "https://www.informatics.jax.org/downloads/reports"

CDIR = G.CACHE / "spatial"
FT = G.CACHE / "spatial_fulltext"
for _d in (CDIR, FT):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# vocabulary
# ---------------------------------------------------------------------------
ZONES = {
    "resting": [r"resting zone", r"reserve zone", r"resting chondrocyte",
                r"round(?:ed)? chondrocyte", r"reserve chondrocyte"],
    "proliferative": [r"proliferat\w* zone", r"proliferat\w* chondrocyte",
                      r"columnar (?:zone|chondrocyte)", r"flattened chondrocyte"],
    "prehypertrophic": [r"pre-?hypertrophic"],
    "hypertrophic": [r"(?<!pre-)(?<!pre)hypertroph\w* zone", r"(?<!pre-)(?<!pre)hypertroph\w* chondrocyte",
                     r"zone of hypertrophy"],
    "terminal_hypertrophic": [r"terminal(?:ly)? hypertroph\w*", r"late hypertroph\w*",
                              r"chondro-?osseous junction"],
    "perichondrial": [r"perichondri\w+", r"groove of ranvier", r"periosteu\w*", r"periosteal"],
}
NON_CHONDRO = {
    "osteoblast": [r"osteoblast\w*", r"primary spongiosa", r"trabecular bone"],
    "osteoclast": [r"osteoclast\w*", r"TRAP-?positive"],
    "vascular": [r"vascular invasion", r"blood vessel", r"endothelial", r"CD31", r"vasculature"],
    "marrow": [r"bone marrow", r"marrow cavity", r"haematopoietic", r"hematopoietic"],
}
METHODS = {
    "RNAscope": [r"RNAscope"],
    "smFISH": [r"smFISH", r"single-?molecule FISH", r"single-?molecule fluorescen\w+ in situ"],
    "in situ hybridization": [r"in situ hybridi[sz]ation", r"\bISH\b", r"riboprobe",
                              r"digoxigenin-?labell?ed probe"],
    "immunohistochemistry": [r"immunohistochemi\w+", r"\bIHC\b", r"DAB stain"],
    "immunofluorescence": [r"immunofluorescen\w+", r"\bIF\b staining", r"immunostain\w+"],
    "spatial transcriptomics": [r"spatial transcriptom\w+", r"Visium", r"MERFISH", r"Slide-?seq",
                                r"Stereo-?seq", r"GeoMx"],
    "reporter mouse": [r"reporter (?:mouse|mice|allele|line)", r"lacZ", r"LacZ",
                       r"knock-?in reporter", r"-?GFP knock-?in", r"tdTomato reporter"],
    "lineage tracing": [r"lineage[- ]trac\w+", r"CreER", r"confetti"],
    "laser capture": [r"laser[- ]capture", r"laser microdissect\w+"],
}
# tissue that must be present for the record to be intact-tissue evidence
INTACT = [r"growth plate", r"growth-plate", r"epiphys\w+", r"physis", r"physeal",
          r"tibial section", r"femoral section", r"metatarsal", r"long bone section",
          r"proximal tibia", r"distal femur", r"cartilage section"]
# hard exclusions - these make a record NOT direct spatial proof
EXCLUDE = {
    "cultured cells": [r"ATDC5", r"chondrogenic cell line", r"micromass culture",
                       r"cultured chondrocytes", r"primary chondrocyte culture",
                       r"monolayer culture", r"C3H10T1/2", r"HEK ?293", r"passaged chondrocytes"],
    "dissociated": [r"single-?cell RNA-?seq", r"scRNA-?seq", r"dissociat\w+ into a single-?cell",
                    r"10x Genomics", r"flow cytometr\w+", r"\bFACS\b", r"fluorescence-activated",
                    r"violin plot", r"\bUMAP\b", r"t-?SNE", r"feature plot", r"dot plot",
                    r"cell clusters?", r"pseudotime", r"gating strateg\w+", r"\bsorted\b"],
    "marker panel / sorting definition": [
        r"(?:CD\d+|TER-?119|Thy1?|PDPN|AlphaV|6C3)\s?[\u2212+-]\s?/",
        r"(?:are|were) defined as .{0,40}(?:SSC|stem cells?|progenitors?)"],
    "bulk without zones": [r"whole cartilage lysate", r"bulk RNA-?seq of cartilage"],
}
# assays that measure the gene but carry no spatial information; if the gene
# sentence names one of these and no spatial method, the figure is not evidence
NON_SPATIAL_ASSAY = [r"immunoblot\w*", r"western blot", r"\bqPCR\b", r"qRT-?PCR",
                     r"RT-?PCR", r"ELISA", r"heatmap", r"bulk RNA-?seq", r"microarray",
                     r"luciferase", r"flow cytometr\w+"]

VALIDATION = {
    "knockout control": [r"knock-?out (?:tissue|control|section)", r"\bKO\b control",
                         r"null (?:tissue|control)", r"-?/-? (?:control )?section",
                         r"specificity was confirmed in", r"absent in .{0,30}knock-?out"],
    "catalogue identifier": [r"\bcat(?:alog(?:ue)?)?\.? ?(?:no\.?|#|number)\s*[:#]?\s*[A-Z0-9-]{4,}",
                             r"RRID\s*:\s*AB_\d+", r"\bab\d{4,6}\b", r"\bsc-\d{4,6}\b",
                             r"ACD ?\d{5,6}", r"Advanced Cell Diagnostics"],
    "sense/negative probe": [r"sense probe", r"negative control probe", r"dapB",
                             r"scrambled probe", r"isotype control", r"secondary[- ]only"],
    "peptide block": [r"peptide (?:block|competition|pre-?absorption)", r"pre-?absorb\w+"],
}
SPECIES = {"mouse": [r"\bmouse\b", r"\bmice\b", r"\bmurine\b", r"C57BL"],
           "rat": [r"\brat\b", r"\brats\b", r"Sprague", r"Wistar"],
           "human": [r"\bhuman\b", r"\bpatient\b", r"fetal human", r"foetal human"],
           "chick": [r"\bchick\b", r"\bchicken\b"],
           "zebrafish": [r"zebrafish", r"Danio"]}
AGE = [r"E\d{1,2}\.?\d?\b", r"P\d{1,3}\b", r"\b\d{1,2}[- ]week[- ]old", r"\b\d{1,2}[- ]day[- ]old",
       r"\b\d{1,2}[- ]month[- ]old", r"embryonic day \d{1,2}", r"postnatal day \d{1,3}",
       r"\b\d{1,2} (?:weeks?|days?|months?) of age", r"gestational week \d{1,2}"]
BONES = [r"tibia\w*", r"femur", r"femoral", r"metatarsal", r"humerus", r"radius", r"ulna",
         r"vertebra\w*", r"rib\b", r"iliac crest", r"fibula"]
INTENSITY = {"strong": [r"\bstrong(?:ly)?\b", r"\bintense(?:ly)?\b", r"\brobust(?:ly)?\b",
                        r"\babundant(?:ly)?\b", r"\bhigh(?:ly)? express\w+"],
             "moderate": [r"\bmoderate(?:ly)?\b", r"\bmodest\b", r"\bdetectable\b"],
             "weak": [r"\bweak(?:ly)?\b", r"\bfaint(?:ly)?\b", r"\blow[- ]level\b",
                      r"\bsparse(?:ly)?\b"],
             "absent": [r"\bnot detect\w+", r"\bno (?:signal|staining|expression)\b",
                        r"\babsent\b", r"\bundetectable\b", r"\bnegative for\b"]}
QUANTIFIED = [r"quantif\w+", r"\bmean (?:intensity|fluorescence)", r"H-?score",
              r"percent\w* (?:of )?(?:positive|cells)", r"puncta per cell", r"dots per cell",
              r"corrected total cell fluorescence", r"\bImageJ\b", r"\bQuPath\b",
              r"integrated density", r"\bn ?= ?\d+ (?:mice|animals|sections|bones)"]


# ---------------------------------------------------------------------------
# localization intent
#
# The gene appearing in a figure caption is not enough. In most captions the
# gene name is a *genotype* - "Sufu f/f", "Itgb1 iDEC", "Gnas R201H" - and the
# figure shows a mutant phenotype, not where the gene is expressed. Those
# figures say nothing about localization and must not be counted as spatial
# evidence. {G} below is substituted with the gene alternation.
# ---------------------------------------------------------------------------
EXPRESSION_CUES = [
    r"expression (?:of|for) {G}", r"{G} (?:mRNA )?expression", r"{G} (?:mRNA|transcripts?)",
    r"{G} protein", r"{G} immunoreactivit\w+", r"{G} (?:immuno)?stain\w+", r"{G} signal",
    r"{G}[- ]positive",
    r"{G}\s?\(?\+\)?\s?(?:cells?|chondrocytes?|nuclei|population|signal|staining|area|region)",
    r"anti-?{G}", r"stained for {G}",
    r"labell?ed (?:for|with) {G}", r"probes? (?:for|against) {G}",
    r"{G} (?:was|were|is|are)\s+(?:strongly |weakly |also |not |only |mainly |predominantly |"
    r"specifically |exclusively |highly |broadly )*(?:detect\w+|express\w+|localis\w+|localiz\w+|"
    r"present|restricted|confined|enriched|observed|seen|found|absent|visible)",
    r"localis\w+ of {G}", r"localiz\w+ of {G}", r"distribution of {G}",
    r"{G}[- ](?:reporter|lacZ|LacZ|GFP|tdTomato|mCherry|EGFP|Venus)",
    r"{G}[- ]?Cre(?:ER(?:T2)?)?", r"{G} in situ", r"in situ[^.]{0,60}{G}",
    r"RNAscope[^.]{0,60}{G}", r"{G}[^.]{0,40}RNAscope", r"{G} riboprobe",
    r"immunostain\w+ for {G}", r"immunohistochemistry for {G}",
    r"{G} was (?:up|down)-?regulated", r"{G} levels",
]
GENOTYPE_CUES = [
    r"{G} ?(?:fl/fl|f/f|flox\w*|fl/\+)", r"{G}[- ]?c?KO\b", r"{G} ?(?:-/-|\+/-|Δ/Δ|Δ\w*)",
    r"{G}[- ](?:mutant|knock-?out|null|deficient|deletion|deleted|haploinsufficient|transgenic)",
    r"(?:deletion|ablation|loss|knock-?down|knock-?out|overexpression) of {G}",
    r"{G}[A-Z]\d{1,4}[A-Z]\b", r"{G}[- ]?i?Δ[A-Za-z]{1,4}", r"{G} ?siRNA", r"{G} ?shRNA",
]


def intent_patterns(mouse: str, human) -> tuple[list, list]:
    alts = {mouse, str(mouse).upper()}
    if human and str(human) != "nan":
        alts.add(str(human))
    G_ = "(?:" + "|".join(re.escape(a) for a in sorted(alts, key=len, reverse=True) if a) + ")"
    return ([p.replace("{G}", G_) for p in EXPRESSION_CUES],
            [p.replace("{G}", G_) for p in GENOTYPE_CUES])


def _any(pats, text) -> bool:
    return any(re.search(p, text, re.I) for p in pats)


def _which(d: dict, text: str) -> list[str]:
    return [k for k, pats in d.items() if _any(pats, text)]


def _first(pats, text):
    for p in pats:
        m = re.search(p, text, re.I)
        if m:
            return m.group(0)
    return None


# ---------------------------------------------------------------------------
# caching
# ---------------------------------------------------------------------------
def cached(key: str, fn):
    f = CDIR / f"{key}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            pass
    v = fn()
    f.write_text(json.dumps(v))
    return v


def _k(prefix: str, s: str) -> str:
    return f"{prefix}_{hashlib.sha1(s.encode()).hexdigest()[:16]}"


# ---------------------------------------------------------------------------
# sources
# ---------------------------------------------------------------------------
def epmc_search(query: str, page_size: int = 100, max_pages: int = 2) -> list[dict]:
    def go():
        out, cursor = [], "*"
        for _ in range(max_pages):
            u = (f"{EPMC}/search?query={urllib.parse.quote(query)}"
                 f"&format=json&pageSize={page_size}&resultType=core&cursorMark={cursor}")
            j = G.get(u, timeout=180).json()
            res = j.get("resultList", {}).get("result", [])
            out.extend(res)
            nxt = j.get("nextCursorMark")
            if not nxt or nxt == cursor or len(res) < page_size:
                break
            cursor = nxt
            time.sleep(0.15)
        return out
    return cached(_k("epmc", query), go)


def pubmed_search(term: str, retmax: int = 40) -> list[str]:
    def go():
        r = G.get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax={retmax}"
                  f"&term={urllib.parse.quote_plus(term)}", timeout=120)
        time.sleep(0.34)
        return r.json().get("esearchresult", {}).get("idlist", []) or []
    try:
        return cached(_k("pm", f"{term}|{retmax}"), go)
    except Exception:  # noqa: BLE001
        return []


def pubmed_count(term: str) -> int:
    def go():
        r = G.get(f"{EUTILS}/esearch.fcgi?db=pubmed&retmode=json&retmax=0"
                  f"&term={urllib.parse.quote_plus(term)}", timeout=120)
        time.sleep(0.34)
        return int(r.json()["esearchresult"]["count"])
    try:
        return cached(_k("pmc", term), go)
    except Exception:  # noqa: BLE001
        return -1


def epmc_count(query: str) -> int:
    def go():
        u = (f"{EPMC}/search?query={urllib.parse.quote(query)}"
             f"&format=json&pageSize=1&resultType=lite")
        return int(G.get(u, timeout=120).json().get("hitCount", 0))
    try:
        return cached(_k("epmcn", query), go)
    except Exception:  # noqa: BLE001
        return -1


def biostudies_search(query: str, page_size: int = 20) -> dict:
    def go():
        u = f"{BIOSTUDIES}/search?query={urllib.parse.quote(query)}&pageSize={page_size}"
        j = G.get(u, timeout=180).json()
        return {"total": j.get("totalHits", 0),
                "hits": [{"accession": h.get("accession"), "title": h.get("title"),
                          "release_date": h.get("release_date")}
                         for h in j.get("hits", [])]}
    try:
        return cached(_k("bst", query), go)
    except Exception:  # noqa: BLE001
        return {"total": -1, "hits": []}


def hpa_record(human_gene: str) -> dict:
    """Human Protein Atlas. HPA has no growth-plate tissue, so this can never be
    direct evidence here - it is recorded as queried and as a negative."""
    def go():
        cols = "g,gs,eg,scml,rnatsm,rnatd,ab"
        u = (f"{HPA}?search={urllib.parse.quote(human_gene)}&format=json"
             f"&columns={cols}&compress=no")
        j = G.get(u, timeout=120).json()
        for rec in j:
            if str(rec.get("Gene", "")).upper() == human_gene.upper():
                return rec
        return {}
    if not human_gene or str(human_gene) == "nan":
        return {}
    try:
        return cached(_k("hpa", human_gene), go)
    except Exception:  # noqa: BLE001
        return {}


# ---------------------------------------------------------------------------
# MGI GXD: curated expression-assay references per gene
# ---------------------------------------------------------------------------
_mgi_cache: dict = {}


def mgi_gxd_pmids() -> dict:
    """{mouse_symbol: [pmid, ...]} from MGI's curated GXD reference report."""
    if _mgi_cache:
        return _mgi_cache
    f = CDIR / "mgi_gxd_pmids.json"
    if f.exists():
        _mgi_cache.update(json.loads(f.read_text()))
        return _mgi_cache
    j2p = {}
    for line in G.get(f"{MGI}/BIB_PubMed.rpt", timeout=600).text.splitlines():
        p = line.split("\t")
        if len(p) >= 3 and p[1].strip().isdigit():
            j2p[p[2].strip()] = p[1].strip()
    out = {}
    for line in G.get(f"{MGI}/MRK_GXD.rpt", timeout=600).text.splitlines():
        p = line.split("\t")
        if len(p) < 3:
            continue
        sym, jnums = p[1].strip(), p[2].strip()
        pmids = [j2p[j] for j in jnums.split(",") if j in j2p]
        if sym and pmids:
            out[sym] = pmids
    f.write_text(json.dumps(out))
    _mgi_cache.update(out)
    return _mgi_cache


# ---------------------------------------------------------------------------
# full text
# ---------------------------------------------------------------------------
def fetch_fulltext(pmcid: str) -> str | None:
    if not pmcid or str(pmcid) == "nan":
        return None
    dest = FT / f"{pmcid}.xml"
    if dest.exists():
        t = dest.read_text(encoding="utf-8", errors="replace")
        return t if len(t) > 2000 else None
    txt = None
    try:
        r = G.get(f"{EPMC}/{pmcid}/fullTextXML", timeout=240, tries=2)
        if r.status_code == 200 and len(r.text) > 2000:
            txt = r.text
    except Exception:  # noqa: BLE001
        pass
    if txt is None:
        dest.write_text("", encoding="utf-8")
        return None
    dest.write_text(txt, encoding="utf-8", errors="replace")
    return txt


_TAG = re.compile(r"<[^>]+>")


def _strip(x: str) -> str:
    return re.sub(r"\s+", " ", _TAG.sub(" ", x)).strip()


def figures(xml: str) -> list[dict]:
    """Every <fig> in the article, as {label, caption}."""
    out = []
    for m in re.finditer(r"<fig\b.*?</fig>", xml, re.S):
        blk = m.group(0)
        lab = re.search(r"<label>(.*?)</label>", blk, re.S)
        cap = re.search(r"<caption>(.*?)</caption>", blk, re.S)
        out.append({"label": _strip(lab.group(1)) if lab else "",
                    "caption": _strip(cap.group(1)) if cap else ""})
    return out


def methods_text(xml: str) -> str:
    """Methods-like sections, where reagent identity and validation live."""
    chunks = []
    for m in re.finditer(r"<sec\b[^>]*>(.*?)</sec>", xml, re.S):
        blk = m.group(1)
        t = re.search(r"<title>(.*?)</title>", blk, re.S)
        title = _strip(t.group(1)) if t else ""
        if re.search(r"method|material|experimental|procedure", title, re.I):
            chunks.append(_strip(blk))
    return " ".join(chunks)[:200000]


def gene_pattern(mouse: str, human) -> re.Pattern:
    alts = {mouse, str(mouse).upper()}
    if human and str(human) != "nan":
        alts.add(str(human))
    alts = [re.escape(a) for a in sorted(alts, key=len, reverse=True) if a]
    return re.compile(r"(?<![A-Za-z0-9])(" + "|".join(alts) + r")(?![A-Za-z0-9])")


def sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.;!?])\s+", text) if s.strip()]


# ---------------------------------------------------------------------------
# MGI mouse phenotypes (MP terms) - real knockout skeletal phenotypes per gene
# ---------------------------------------------------------------------------
_pheno_cache: dict = {}

SKELETAL_MP = [
    r"body length", r"body size", r"\bstature\b", r"\bdwarf", r"\bdwarfis",
    r"long bone", r"growth plate", r"chondrocyte", r"cartilage", r"ossification",
    r"\blimb\b", r"limbs\b", r"\bfemur\b", r"\btibia\b", r"\bhumerus\b",
    r"skeleton", r"skeletal", r"bone mineral", r"endochondral", r"chondrodysplasia",
    r"achondroplasia", r"osteopenia", r"osteopetro", r"epiphys", r"physeal",
    r"\bvertebra", r"\bsnout\b", r"craniofacial",
]
LENGTH_DIRECTION = {
    "shorter": [r"decreased body length", r"decreased body size", r"short(?:ened)? limbs?",
                r"decreased length of long bones", r"dwarf", r"chondrodysplasia",
                r"short(?:ened)? (?:femur|tibia|humerus)", r"decreased (?:femur|tibia) length",
                r"disproportionate dwarf", r"micromelia", r"achondroplasia"],
    "longer": [r"increased body length", r"increased body size",
               r"increased length of long bones", r"long(?:er)? limbs?",
               r"increased (?:femur|tibia) length"],
    "disorganized": [r"abnormal growth plate", r"disorganized (?:growth plate|chondrocyte)",
                     r"abnormal chondrocyte (?:morphology|differentiation|proliferation)",
                     r"abnormal endochondral", r"abnormal ossification",
                     r"abnormal cartilage", r"abnormal skeleton morphology"],
    "lethal": [r"(?:embryonic|perinatal|postnatal|prenatal|neonatal) lethality",
               r"complete lethality", r"lethality.{0,20}(?:complete|incomplete)"],
}


def mgi_phenotypes() -> dict:
    """{mouse_symbol: [{'mp': ..., 'term': ..., 'pmid': ...}, ...]} from MGI reports."""
    if _pheno_cache:
        return _pheno_cache
    f = CDIR / "mgi_phenotypes.json"
    if f.exists():
        _pheno_cache.update(json.loads(f.read_text()))
        return _pheno_cache
    mp2term = {}
    for line in G.get(f"{MGI}/VOC_MammalianPhenotype.rpt", timeout=600).text.splitlines():
        p = line.split("\t")
        if len(p) >= 2 and p[0].startswith("MP:"):
            mp2term[p[0]] = p[1].strip()
    mgi2sym = {}
    for line in G.get(f"{MGI}/MRK_List2.rpt", timeout=600).text.splitlines():
        p = line.split("\t")
        if len(p) >= 7 and p[0].startswith("MGI:"):
            mgi2sym[p[0]] = p[6].strip()
    out: dict = {}
    for line in G.get(f"{MGI}/MGI_GenePheno.rpt", timeout=900).text.splitlines():
        p = line.split("\t")
        if len(p) < 7 or not p[4].startswith("MP:"):
            continue
        for acc in p[6].split(","):
            sym = mgi2sym.get(acc.strip())
            if not sym:
                continue
            term = mp2term.get(p[4], "")
            out.setdefault(sym, [])
            key = (p[4], p[5].strip())
            if not any(x["mp"] == key[0] and x["pmid"] == key[1] for x in out[sym]):
                out[sym].append({"mp": p[4], "term": term, "pmid": p[5].strip(),
                                 "allele": p[1].strip(), "genotype": p[0].strip(),
                                 "background": p[3].strip()})
    f.write_text(json.dumps(out))
    _pheno_cache.update(out)
    return _pheno_cache


def skeletal_phenotypes(symbol: str) -> list[dict]:
    return [x for x in mgi_phenotypes().get(symbol, [])
            if _any(SKELETAL_MP, x.get("term", ""))]


def length_direction(terms: list[str]) -> dict:
    joined = " | ".join(terms)
    return {k: sorted({_first([p], joined) for p in pats if _first([p], joined)})
            for k, pats in LENGTH_DIRECTION.items()}
