#!/usr/bin/env python3
"""Map Davis et al. 2020 cell types (PD-gene expression table) to FlyWire neurons.

Inputs
  --annotations  FlyWire neuron annotation table
                 (flyconnectome/flywire_annotations, supplemental_files/
                  Supplemental_file1_neuron_annotations (2).tsv)
  --expression   expression_table.csv (experiment, cell_type, expression, gene, relative)
  --clusters     clusterFBbt.csv (cluster, fbbt_id)
Outputs (in --outdir)
  celltype_flywire_mapping.csv   one row per expression cell type, with confidence
  flywire_roots_by_celltype.csv  one row per FlyWire neuron per matched cell type
  pd_gene_expression_head.csv    cell type x gene matrix (head experiment only)
  dopamine_seed_neurons.csv      FlyWire dopaminergic neurons (DAN class + PAM/PPL/PPM names)

Matching order: explicit rule (FlyWire label / family / class) -> exact label from
the cell type name -> no match. Every result is cross-checked against the fbbt_id
column of the FlyWire table.
"""
import argparse
import os
import re

import pandas as pd

# File paths. Relative paths are looked up in the folder this script lives in, so the
# script runs with PyCharm's Run button and no arguments. Edit the defaults to match
# where your files are, or override any of them on the command line.
ap = argparse.ArgumentParser()
ap.add_argument("--annotations", default="Supplemental_file1_neuron_annotations (2).tsv")
ap.add_argument("--expression", default="expression_table.csv")
ap.add_argument("--clusters", default="clusterFBbt.csv")
ap.add_argument("--outdir", default=".", help="where the four result files are written (relative to this script's folder)")
args = ap.parse_args()

BASE = os.path.dirname(os.path.abspath(__file__))


def resolve(p):
    """Use the path as given if it exists (or is absolute); otherwise look next to this script."""
    return p if os.path.isabs(p) or os.path.exists(p) else os.path.join(BASE, p)


for _name in ["annotations", "expression", "clusters"]:
    setattr(args, _name, resolve(getattr(args, _name)))
    if not os.path.exists(getattr(args, _name)):
        raise SystemExit(f"file not found: {getattr(args, _name)}\n"
                         f"Put it there, or change the default for --{_name} near the top of the script.")
args.outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(BASE, args.outdir)
os.makedirs(args.outdir, exist_ok=True)


def norm(s):
    s = s.replace("’", "'")
    s = re.sub(r"\s*/\s*", "/", s)
    s = re.sub(r"\s+", " ", s).strip().lower()
    return s.replace("centrigufal", "centrifugal")


# ---------------------------------------------------------------- load data
ann = pd.read_csv(
    args.annotations, sep="\t", low_memory=False,
    usecols=["root_id", "side", "super_class", "cell_class", "cell_type",
             "hemibrain_type", "top_nt", "fbbt_id"],
)
ann["label"] = ann["cell_type"].fillna(ann["hemibrain_type"])
label = ann["label"]

expr = pd.read_csv(args.expression)
expr["exp"] = expr["experiment"].str.extract(r"- (head|appendage)$")[0]
head = expr[expr["exp"] == "head"].copy()
cell_types = sorted(head["cell_type"].unique())

clusters = pd.read_csv(args.clusters)
clusters["key"] = clusters["cluster"].map(norm)
fbbt_by_key = dict(zip(clusters["key"], clusters["fbbt_id"]))

# -------------------------------------------------------------------- rules
rules = {}


def add(name, conf, method, note="", labels=None, regex=None, cell_class=None,
        super_class=None, top_nt=None, parent=False, kind="neuron"):
    rules[norm(name)] = dict(conf=conf, method=method, note=note, labels=labels,
                             regex=regex, cell_class=cell_class,
                             super_class=super_class, top_nt=top_nt,
                             parent=parent, kind=kind)


# Dm
add("distal medullary amacrine neuron Dm3", "medium", "label_family", regex=r"^Dm3[pqv]$", parent=True,
    note="FlyWire splits Dm3 into Dm3p/Dm3q/Dm3v; pooled here")
# Kenyon cells
add("alpha/beta Kenyon cell", "medium", "label_family", regex=r"^KCab", parent=True,
    note="pooled FlyWire alpha/beta types (KCab, KCab-p, ...)")
add("alpha/beta posterior Kenyon cell", "high", "label", labels=["KCab-p"])
add("alpha/beta surface Kenyon cell", "low", "label", labels=["KCab"],
    note="FlyWire does not separate surface/core/middle alpha/beta; KCab is the pooled alpha/beta type")
add("alpha/beta core Kenyon cell", "low", "label", labels=["KCab"],
    note="FlyWire does not separate surface/core/middle alpha/beta; KCab is the pooled alpha/beta type. Not in your cluster list")
add("gamma dorsal Kenyon cell", "high", "label", labels=["KCg-d"])
add("alpha'/beta' anterior-posterior Kenyon cell", "medium", "label_family", labels=["KCapbp-ap1", "KCapbp-ap2"],
    note="assumes FlyWire 'ap' = anterior-posterior; verify")
add("alpha'/beta' middle Kenyon cell", "medium", "label", labels=["KCapbp-m"],
    note="assumes FlyWire KCapbp-m = alpha'/beta' middle; verify")
# PAM (Davis/FBbt has 9 numbered types; FlyWire has 15 PAM types)
pam_note = "FlyWire has more PAM types than this table; numbering may not correspond one-to-one"
add("dopaminergic PAM neuron 1", "low", "label_numbering", labels=["PAM01"], note=pam_note + "; no FBbt support")
add("dopaminergic PAM neuron 2", "medium", "label_fbbt", labels=["PAM02"], note=pam_note + "; supported by FBbt id")
add("dopaminergic PAM neuron 3", "medium", "label_fbbt", labels=["PAM03"], note=pam_note + "; supported by FBbt id")
add("dopaminergic PAM neuron 4", "low", "label_numbering", labels=["PAM04"], note=pam_note + "; no FBbt support")
add("dopaminergic PAM neuron 5", "medium", "label_fbbt", labels=["PAM05"], note=pam_note + "; supported by FBbt id")
add("dopaminergic PAM neuron 6", "low", "label_numbering", labels=["PAM06"],
    note=pam_note + "; FlyWire's fbbt_id for this term points to PAM02, so the numbering is inconsistent")
add("dopaminergic PAM neuron 7", "medium", "label_fbbt", labels=["PAM07"], note=pam_note + "; supported by FBbt id")
add("dopaminergic PAM neuron 8", "medium", "label_fbbt", labels=["PAM08"], note=pam_note + "; supported by FBbt id")
add("dopaminergic PAM neuron 9", "medium", "label_fbbt", labels=["PAM09"],
    note=pam_note + "; FBbt id also covers a few PAM10 neurons")
add("mushroom body output neuron", "medium", "cell_class", cell_class="MBON", parent=True,
    note="pooled: one expression row for all ~96 FlyWire MBONs")
# optic lobe, naming differences
add("lobula columnar neuron LC10A", "low", "label", labels=["LC10a"],
    note="capital-A term; FlyWire has LC10a-f. Probably overlaps the LC10a row")
add("lobula plate intrinsic neuron LPi3-4", "low", "label_numbering", labels=["LPi03", "LPi04"],
    note="FlyWire zero-pads LPi numbers; assumed LPi3-4 = LPi03 + LPi04")
add("proximal medullary amacrine neuron Pm3", "low", "label_numbering", labels=["Pm03"],
    note="FlyWire zero-pads Pm numbers; assumed Pm3 = Pm03")
add("proximal medullary amacrine neuron Pm4", "low", "label_numbering", labels=["Pm04"],
    note="FlyWire zero-pads Pm numbers; assumed Pm4 = Pm04")
add("transmedullary neuron TM3a", "medium", "label", labels=["Tm3"], note="FlyWire does not subdivide Tm3")
add("transmedullary neuron Tm3Y", "low", "label", labels=["TmY3"], note="assumed older name for TmY3")
add("lamina intrinsic neuron", "high", "label", labels=["Lai"])
add("lamina wide-field 1 neuron", "high", "label", labels=["Lawf1"])
add("lamina wide-field 2 neuron", "high", "label", labels=["Lawf2"])
add("lamina monopolar neuron", "medium", "label_family", labels=["L1", "L2", "L3", "L4", "L5"], parent=True,
    note="generic parent of L1-L5")
add("centrifugal neuron", "low", "label_family", labels=["C2", "C3"], super_class="visual_centrifugal", parent=True,
    note="generic parent: C2, C3 plus FlyWire visual_centrifugal super class")
add("T neuron", "medium", "label_family", regex=r"^T[1-5][a-d]?$", parent=True,
    note="generic parent: T1-T5 families")
add("T neuron T4", "medium", "label_family", regex=r"^T4[a-d]$", parent=True, note="FlyWire splits T4 into T4a-d; pooled")
add("T neuron T5", "medium", "label_family", regex=r"^T5[a-d]$", parent=True, note="FlyWire splits T5 into T5a-d; pooled")
add("lobula plate tangential neuron", "low", "label_family", regex=r"^(HS[ENS]|VS\d|H1|H2|CH)$", parent=True,
    note="generic class; only HS/VS/H1/H2/CH labels that exist are included")
add("outer photoreceptor cell", "high", "label", labels=["R1-6"], note="FlyWire's R1-6 type is the outer photoreceptor class")
add("photoreceptor cell R8 of pale ommatidium", "medium", "label", labels=["R8"],
    note="FlyWire does not split R8 into pale/yellow; same neurons as the yellow row")
add("photoreceptor cell R8 of yellow ommatidium", "medium", "label", labels=["R8"],
    note="FlyWire does not split R8 into pale/yellow; same neurons as the pale row")
# clock / peptidergic
add("I-LNv neuron", "high", "label_fbbt", labels=["l-LNv"], note="FlyWire writes it with a lowercase l; supported by FBbt id")
add("Pdf neuron", "medium", "label_family", regex=r"^(l-LNv|s-LNv)", parent=True, note="PDF-expressing clock neurons")
add("corazonin neuron", "medium", "label_family", regex=r"^CRZ", note="FlyWire CRZ types")
add("Drosulfakinin neuron", "low", "label_family", regex=r"^DSK", note="name suggests Dsk neurons; verify")
add("adult NPF neuron", "low", "label_family", regex=r"^NPF",
    note="only a small lateral NPF subset is labelled in FlyWire. Not in your cluster list")
# neurotransmitter classes
for nm, nt in [("cholinergic neuron", "acetylcholine"), ("GABAergic neuron", "gaba"), ("glutamatergic neuron", "glutamate")]:
    add(nm, "class-level", "predicted_nt", top_nt=nt, parent=True, kind="nt_class",
        note="FlyWire top_nt is a computational prediction, not a cell type; kept out of the roots file (filter by top_nt)")
# no counterpart
add("adult extrinsic protocerebral bridge neuron", "none", "manual", kind="generic",
    note="generic class; would be several FlyWire central complex columnar types - choose manually")
add("adult CCAP neuron", "none", "no_match", kind="generic", note="no CCAP-labelled neurons in the FlyWire annotation table")
for nm in ["adult lamina epithelial glial cell", "adult lamina marginal glial cell", "pseudocartridge glial cell"]:
    add(nm, "none", "not_in_dataset", kind="non_neuron", note="glia are not segmented as neurons in this table")
add("muscle cell", "none", "not_in_dataset", kind="non_neuron", note="no muscle in the FlyWire brain dataset")


def default_rule(name):
    tok = name.split()[-1]
    if re.fullmatch(r"[A-Za-z]+\d+[a-z]?", tok) and (label == tok).any():
        return dict(conf="high", method="label", note="", labels=[tok], regex=None,
                    cell_class=None, super_class=None, top_nt=None, parent=False, kind="neuron")
    return dict(conf="none", method="no_match", note=f"no FlyWire label '{tok}'", labels=None,
                regex=None, cell_class=None, super_class=None, top_nt=None, parent=False, kind="neuron")


def select(r):
    m = pd.Series(False, index=ann.index)
    if r["labels"]:
        m |= label.isin(r["labels"])
    if r["regex"]:
        m |= label.fillna("").str.match(r["regex"])
    if r["cell_class"]:
        m |= ann["cell_class"].eq(r["cell_class"])
    if r["super_class"]:
        m |= ann["super_class"].eq(r["super_class"])
    if r["top_nt"]:
        m |= ann["top_nt"].eq(r["top_nt"])
    return m


rows, roots = [], []
for ct in cell_types:
    key = norm(ct)
    r = dict(rules.get(key) or default_rule(ct))
    sel = select(r)
    n = int(sel.sum())
    note = r["note"]
    conf = r["conf"]
    if conf == "high" and n < 10 and r["method"] == "label":
        conf = "low"
        note = (note + "; " if note else "") + f"only {n} FlyWire neurons carry this label, likely incomplete"
    if n == 0 and conf not in ("none",):
        conf = "none"
        note = (note + "; " if note else "") + "selection matched no FlyWire neurons"
    fb = fbbt_by_key.get(key, "")
    fb_n = int((ann["fbbt_id"] == fb).sum()) if fb else 0
    fb_cov = float((sel & (ann["fbbt_id"] == fb)).sum() / fb_n) if fb_n else None
    if fb_n and fb_cov is not None and fb_cov < 0.8 and r["kind"] == "neuron" and n > 0:
        note = (note + "; " if note else "") + (
            f"FBbt id {fb} is carried by {fb_n} FlyWire neurons, only {fb_cov:.0%} of them in this selection (shared id)")
    labs = sorted(label[sel].dropna().unique())
    rows.append(dict(
        cell_type=ct, fbbt_id=fb, kind=r["kind"], flywire_labels=";".join(labs) if len(labs) <= 25 else f"{len(labs)} labels",
        n_labels=len(labs), n_neurons=n, method=r["method"], confidence=conf, is_parent_class=r["parent"],
        fbbt_in_flywire_n=fb_n, fbbt_covered_frac=round(fb_cov, 3) if fb_cov is not None else None, notes=note))
    if n and r["kind"] != "nt_class":
        sub = ann.loc[sel, ["root_id", "label", "side", "super_class", "cell_class", "top_nt"]].copy()
        sub.insert(0, "confidence", conf)
        sub.insert(0, "cell_type", ct)
        roots.append(sub.rename(columns={"label": "flywire_label"}))

mapping = pd.DataFrame(rows)
mapping.to_csv(os.path.join(args.outdir, "celltype_flywire_mapping.csv"), index=False)
pd.concat(roots).to_csv(os.path.join(args.outdir, "flywire_roots_by_celltype.csv"), index=False)

# expression matrix (head only): raw expression and per-gene relative (min-max) values
wide = head.pivot_table(index="cell_type", columns="gene", values="expression", aggfunc="first")
rel = head.pivot_table(index="cell_type", columns="gene", values="relative", aggfunc="first").add_prefix("rel_")
wide.join(rel).to_csv(os.path.join(args.outdir, "pd_gene_expression_head.csv"))

# dopaminergic seed neurons
dan = ann[(ann["cell_class"] == "DAN") | label.fillna("").str.match(r"^(PAM|PPL|PAL|PPM)\d")]
dan[["root_id", "label", "side", "cell_class", "super_class", "top_nt", "fbbt_id"]].rename(
    columns={"label": "flywire_label"}).to_csv(os.path.join(args.outdir, "dopamine_seed_neurons.csv"), index=False)

# ------------------------------------------------------------------ summary
print("cell types:", len(mapping))
print(mapping.groupby(["confidence"]).size().to_string())
print()
print(mapping.groupby("kind").size().to_string())
print("\nneurons per confidence (leaf, non-parent rows):")
print(mapping[~mapping.is_parent_class].groupby("confidence").n_neurons.sum().to_string())
print("\nunique FlyWire neurons covered:", pd.concat(roots)["root_id"].nunique())
print("\nDAN seed counts by label prefix:")
print(dan["label"].fillna("").str.extract(r"^([A-Za-z]+\d?)")[0].value_counts().to_string())


#!/usr/bin/env python3
"""Downstream targets of dopaminergic neurons, scored by Parkinson's-gene expression.

Reads FlyWire's connections table (columns: pre_root_id, post_root_id, neuropil,
syn_count, nt_type), keeps connections that start at the dopamine seed neurons,
rolls the partners up to FlyWire cell types, and scores each partner by the PD-gene
expression of the Davis et al. cell type it maps to.

Outputs (in --outdir)
  downstream_by_flywire_type.csv  one row per partner FlyWire cell type
  seed_group_summary.csv          outgoing synapses per dopamine group (PAM, PPL1, ...)

Columns in downstream_by_flywire_type.csv
  syn_from_seeds       total synapses received from the dopamine seeds (after the --min-syn cutoff)
  syn_from_<group>     the same, split by seed group
  pd_score             synapse-weighted mean PD-gene score of the partner neurons that could be scored
  frac_syn_scored      share of this type's input that went to neurons with a Davis expression row
  impact               sum over partner neurons of (synapses from seeds x pd_score); unscored neurons add 0
  davis_rows           the Davis cell type row(s) the score came from
Types that have no Davis expression row are NOT scored; they appear with empty pd_score.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
# File paths. Relative paths are looked up in the folder this script lives in, so the
# script runs with PyCharm's Run button and no arguments. Edit the defaults to match
# where your files are, or override any of them on the command line.
ap.add_argument("--connections", default="classificationNeuronData/connections_princeton.csv")
ap.add_argument("--seeds", default="dopamine_seed_neurons.csv")
ap.add_argument("--annotations", default="Supplemental_file1_neuron_annotations (2).tsv",
                help="table with root_id and cell_type (and optionally hemibrain_type); .tsv or .csv")
ap.add_argument("--roots", default="flywire_roots_by_celltype.csv")
ap.add_argument("--mapping", default="celltype_flywire_mapping.csv")
ap.add_argument("--expression", default="pd_gene_expression_head.csv")
ap.add_argument("--outdir", default="results")
ap.add_argument("--min-syn", type=int, default=5, help="min synapses per neuron pair (summed over neuropils)")
ap.add_argument("--genes", default="pink1,park,Lrrk,DJ1alpha,DJ1beta",
                help="comma-separated genes that define the PD score")
ap.add_argument("--chunksize", type=int, default=2_000_000)
ap.add_argument("--pre-col", default="pre_root_id")
ap.add_argument("--post-col", default="post_root_id")
ap.add_argument("--syn-col", default="syn_count")
a = ap.parse_args()

BASE = os.path.dirname(os.path.abspath(__file__))


def resolve(p):
    """Use the path as given if it exists (or is absolute); otherwise look next to this script."""
    return p if os.path.isabs(p) or os.path.exists(p) else os.path.join(BASE, p)


for _name in ["connections", "seeds", "annotations", "roots", "mapping", "expression", "outdir"]:
    setattr(a, _name, resolve(getattr(a, _name)))
for _name in ["connections", "seeds", "annotations", "roots", "mapping", "expression"]:
    if not os.path.exists(getattr(a, _name)):
        sys.exit(f"file not found: {getattr(a, _name)}\n"
                 f"Put it there, or change the default for --{_name} near the top of the script.")
os.makedirs(a.outdir, exist_ok=True)

# ------------------------------------------------------------------- seeds
seeds = pd.read_csv(a.seeds)
seeds["group"] = seeds["flywire_label"].fillna("").str.extract(r"^(PAM|PPL1|PPL2|PPM|PAL)")[0].fillna("other")
seed_arr = np.array(sorted(seeds["root_id"].astype("int64").unique()), dtype=np.int64)
print(f"seed neurons: {len(seed_arr)}  groups: {seeds.group.value_counts().to_dict()}")

# ---------------------------------------------- stream the connections table
keep, seen = [], set()
n_rows = 0
for chunk in pd.read_csv(a.connections, usecols=[a.pre_col, a.post_col, a.syn_col],
                         dtype={a.pre_col: "int64", a.post_col: "int64"}, chunksize=a.chunksize):
    n_rows += len(chunk)
    sub = chunk[chunk[a.pre_col].isin(seed_arr)]
    if len(sub):
        keep.append(sub)
        seen.update(sub[a.pre_col].unique().tolist())
print(f"connection rows read: {n_rows:,}")
if not keep:
    sys.exit("None of the seed root IDs appear as presynaptic neurons in this connections file. "
             "The files are probably from different FlyWire releases, or the IDs are written differently "
             "(some releases drop the common 720575940 prefix). Compare a few IDs by eye.")
print(f"seed neurons with outgoing connections in this file: {len(seen)} of {len(seed_arr)}")
if len(seen) < 0.9 * len(seed_arr):
    print("WARNING: many seed IDs are missing from the connections file - possible release mismatch.")

edges = pd.concat(keep)
edges = edges.groupby([a.pre_col, a.post_col], as_index=False)[a.syn_col].sum()
edges.columns = ["pre", "post", "syn"]
edges = edges[edges["syn"] >= a.min_syn]
n_seed_targets = int(edges["post"].isin(seed_arr).sum())
edges = edges[~edges["post"].isin(seed_arr)]  # seeds are lost together, so they are not downstream targets
print(f"pairs with >= {a.min_syn} synapses: {len(edges):,} (dropped {n_seed_targets} seed-to-seed pairs)")

edges = edges.merge(seeds[["root_id", "group"]], left_on="pre", right_on="root_id").drop(columns="root_id")

# ------------------------------------------------------ partner cell types
sep = "\t" if a.annotations.lower().endswith(".tsv") else ","
cols = list(pd.read_csv(a.annotations, sep=sep, nrows=0).columns)
lab_cols = [c for c in ["cell_type", "hemibrain_type"] if c in cols]
if "root_id" not in cols or not lab_cols:
    sys.exit(f"annotation file needs a root_id column and a cell_type column; it has: {cols}")
ann = pd.read_csv(a.annotations, sep=sep, usecols=["root_id"] + lab_cols, low_memory=False)
ann["flywire_type"] = ann[lab_cols].bfill(axis=1).iloc[:, 0].fillna("unlabelled")
edges = edges.merge(ann[["root_id", "flywire_type"]], left_on="post", right_on="root_id", how="left").drop(columns="root_id")
edges["flywire_type"] = edges["flywire_type"].fillna("not_in_annotations")

# -------------------------- PD score per partner neuron (Davis expression rows)
roots = pd.read_csv(a.roots, usecols=["cell_type", "root_id"])
mp = pd.read_csv(a.mapping, usecols=["cell_type", "kind", "confidence", "is_parent_class", "n_neurons"])
mp = mp[(mp["kind"] == "neuron") & (mp["confidence"] != "none")]
r = roots.merge(mp, on="cell_type")
r["is_parent_class"] = r["is_parent_class"].astype(bool)
# prefer specific rows: leaf rows if the neuron has any, otherwise its most specific parent row
r["has_leaf"] = (~r["is_parent_class"]).groupby(r["root_id"]).transform("any")
r = r[(~r["is_parent_class"]) | (~r["has_leaf"])]
r["min_n"] = r.groupby("root_id")["n_neurons"].transform("min")
r = r[(~r["is_parent_class"]) | (r["n_neurons"] == r["min_n"])]

ex = pd.read_csv(a.expression, index_col=0)
genes = [g.strip() for g in a.genes.split(",") if g.strip()]
missing = [g for g in genes if f"rel_{g}" not in ex.columns]
if missing:
    sys.exit(f"genes not in the expression matrix: {missing}")
score = ex[[f"rel_{g}" for g in genes]].mean(axis=1).rename("pd_score")
r = r.merge(score, left_on="cell_type", right_index=True)
neuron_score = r.groupby("root_id")["pd_score"].mean()
neuron_rows = r.groupby("root_id")["cell_type"].agg(lambda s: ";".join(sorted(set(s))))

edges["pd_score"] = edges["post"].map(neuron_score)
edges["davis_rows"] = edges["post"].map(neuron_rows)
edges["impact"] = edges["syn"] * edges["pd_score"].fillna(0.0)
edges["syn_scored"] = np.where(edges["pd_score"].notna(), edges["syn"], 0)
edges["w_score"] = edges["syn"] * edges["pd_score"].fillna(0.0)

# ----------------------------------------------------------------- outputs
g = edges.groupby("flywire_type")
by_type = g.agg(syn_from_seeds=("syn", "sum"), n_partner_neurons=("post", "nunique"),
                n_seed_neurons=("pre", "nunique"), syn_scored=("syn_scored", "sum"),
                w_score=("w_score", "sum"), impact=("impact", "sum"))
by_type["frac_syn_scored"] = by_type["syn_scored"] / by_type["syn_from_seeds"]
by_type["pd_score"] = by_type["w_score"] / by_type["syn_scored"].replace(0, np.nan)
by_type["davis_rows"] = g["davis_rows"].agg(lambda s: s.dropna().mode().iat[0] if s.notna().any() else "")
per_group = edges.pivot_table(index="flywire_type", columns="group", values="syn", aggfunc="sum", fill_value=0)
by_type = by_type.join(per_group.add_prefix("syn_from_"))
by_type = by_type.drop(columns=["syn_scored", "w_score"]).sort_values("syn_from_seeds", ascending=False)
by_type.round(4).to_csv(os.path.join(a.outdir, "downstream_by_flywire_type.csv"))

grp = edges.groupby("group").agg(syn_out=("syn", "sum"), n_target_neurons=("post", "nunique"))
grp["n_seed_neurons"] = seeds[seeds["root_id"].isin(seen)].groupby("group")["root_id"].nunique()
grp.to_csv(os.path.join(a.outdir, "seed_group_summary.csv"))

tot = edges["syn"].sum()
print(f"\ntotal synapses from seeds to non-seed partners: {int(tot):,}")
print(f"share of those that land on neurons with a Davis expression row: {edges['syn_scored'].sum() / tot:.1%}")
print("\nseed groups:\n", grp.to_string())
print("\ntop 15 partner types by synapses:\n",
      by_type[["syn_from_seeds", "n_partner_neurons", "pd_score"]].head(15).to_string())
print("\ntop 15 partner types by impact (synapses x PD score):\n",
      by_type.sort_values("impact", ascending=False)[["syn_from_seeds", "pd_score", "impact", "davis_rows"]].head(15).to_string())