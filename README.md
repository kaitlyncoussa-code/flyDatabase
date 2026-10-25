# flyDatabase

# Parkinson's Gene Expression and FlyWire Neuron Database

A database linking Parkinson's disease genes to their fruit fly orthologs, the fly cell types where
those genes are expressed, and the FlyWire neuron types those cell types map to. It is a data
resource for the main project: [link to main repo].

**Final database: `combinedDatabase.csv`** (848 rows, 20 columns)

## What's in the final database

- **14 human Parkinson's genes**, mapped to **16 fly orthologs**.
- **53 FlyWire neuron types**. There is one row for every fly gene and FlyWire type pair
  (16 genes x 53 types = 848 rows).
- Expression for each gene in the matching cell type (Davis et al., 2020, via Expression Atlas).
- Synapse counts from dopaminergic seed neurons to each FlyWire type.

## Columns

| Column | Description |
|---|---|
| `human_gene` | Human Parkinson's gene |
| `fly_gene` | Fly ortholog |
| `flybase_id` | FlyBase ID of the fly gene |
| `ortholog_rating` | Ortholog confidence from the DRSC ortholog tool (high, moderate, or moderate/low) |
| `flywire_type` | FlyWire neuron type |
| `davis_cell_type` | Matching cell type in the Davis et al. 2020 expression data |
| `fbbt_id` | FlyBase anatomy (FBbt) ID of the cell type, blank where none was matched |
| `expression` | Expression of the fly gene in that cell type |
| `relative_expression` | Expression scaled 0 to 1 within each gene |
| `perturbation` | Perturbation value used for the type (equal to `relative_expression` in every row) |
| `syn_from_seeds` | Synapses from dopaminergic seed neurons to this type, counting only neuron pairs with at least 5 synapses (the sum of the five group columns below) |
| `n_partner_neurons` | Number of neurons of this FlyWire type that receive synapses from the seed neurons |
| `n_seed_neurons` | Number of dopaminergic seed neurons that connect to this type |
| `impact` | Sum over the type's partner neurons of synapses from the seeds x `pd_score`. Neurons without an expression match add 0 |
| `pd_score` | Score for the cell type: the mean `relative_expression` of pink1, park, Lrrk, DJ1alpha, and DJ1beta (the script's default gene set), weighted by synapses across the partner neurons |
| `syn_from_PAL`, `syn_from_PAM`, `syn_from_PPL1`, `syn_from_PPL2`, `syn_from_PPM` | Synapses from each dopaminergic neuron group |

## Genes

Orthologs were found with the DRSC ortholog tool (https://www.flyrnai.org/cgi-bin/DRSC_orthologs.pl,
paper: https://pubmed.ncbi.nlm.nih.gov/21880147/). Only the highest-rated ortholog(s) per human gene
were kept.

| Human gene | UniProt | Fly gene | FlyBase ID | Rating |
|---|---|---|---|---|
| LRRK2 | Q5S007 | Lrrk | FBgn0038816 | moderate |
| VPS35 | Q96QK1 | Vps35 | FBgn0034708 | high |
| PRKN | O60260 | park | FBgn0041100 | high |
| PINK1 | Q9BXM7 | pink1 | FBgn0029891 | high |
| PARK7 | Q99497 | DJ1beta | FBgn0039802 | high |
| PARK7 | Q99497 | DJ1alpha | FBgn0033885 | high |
| ATP13A2 | Q9NQ11 | anne | FBgn0052000 | moderate/low |
| PLA2G6 | O60733 | iPLAVIA | FBgn0036053 | high |
| FBXO7 | Q9Y3I1 | ntc | FBgn0035461 | high |
| DNAJC6 | O75061 | aux | FBgn0037218 | moderate/low |
| SYNJ1 | O43426 | Synj | FBgn0034691 | high |
| VPS13C | Q709C8 | Vps13 | FBgn0033194 | high |
| TAF1 | P21675 | taf1 | FBgn0010355 | high |
| RAB39B | Q96DA2 | Rab39 | FBgn0029959 | moderate |
| GBA1 | P04062 | Gba1a | FBgn0051148 | high |
| GBA1 | P04062 | Gba1b | FBgn0051414 | high |

## How it was built

1. **Expression table** (`expressionDatabase.py`): reads one Expression Atlas file per fly gene from
   `expressionRawData/`, combines them into one table (experiment, cell type, expression, gene,
   relative), and saves it as `expression_table.csv` and the SQLite database `pd_expression.db`.
2. **Cell types to neurons** (`searchData.py`, first part): matches each expression cell type to
   FlyWire neurons using the FlyWire annotation table, checking each match against its FBbt ID.
   Outputs `celltype_flywire_mapping.csv`, `flywire_roots_by_celltype.csv`,
   `pd_gene_expression_head.csv`, and `dopamine_seed_neurons.csv`.
3. **Downstream of dopaminergic neurons** (`searchData.py`, second part): reads the FlyWire
   connections table, keeps connections starting at the dopaminergic seed neurons (PAL, PAM, PPL1,
   PPL2, PPM), rolls the partner neurons up to FlyWire cell types, and scores each type by the
   expression of the cell type it maps to. Outputs `results/downstream_by_flywire_type.csv`
   (1,560 types, 53 of them matched to an expression cell type) and `results/seed_group_summary.csv`.
4. The expression, mapping, and synapse tables were combined into `combinedDatabase.csv`, which keeps
   the 53 matched types.

## Input files

- `expressionRawData/`: Expression Atlas downloads, one TSV per fly gene
- `Supplemental_file1_neuron_annotations (2).tsv`: FlyWire neuron annotations
- `clusterFBbt.csv`: cell type clusters with their FBbt IDs, looked up by hand on Virtual Fly Brain
- `classificationNeuronData/connections_princeton.csv`: FlyWire connections table (about 270 MB,
  over GitHub's 100 MB limit, so it is not stored in this repo). Download it from the FlyWire Codex link
  below and place it in `classificationNeuronData/`, or add it with Git LFS.

## Other data files

- `downstream_by_flywire_type.csv`: synapses from the dopaminergic seed neurons to each FlyWire type,
  split by seed group
- `seed_group_summary.csv`: for each seed group, total outgoing synapses, number of target neurons,
  and number of seed neurons

## Setup

```
pip install pandas
```

The scripts read and write relative to the folder they are run from, so run them from the
project root.

## Data sources

- Davis et al., 2020 (basis of the expression data): https://europepmc.org/article/MED/31939737
- EMBL-EBI Expression Atlas (readable version of that data): https://www.ebi.ac.uk/gxa/home
- Virtual Fly Brain: https://v2.virtualflybrain.org
- FlyWire Codex downloads: https://codex.flywire.ai/api/download?dataset=fafb
- FlyWire neuron annotations: https://github.com/flyconnectome/flywire_annotations/blob/main/supplemental_files/Supplemental_file1_neuron_annotations.tsv
- Ortholog lookup (DRSC): https://www.flyrnai.org/cgi-bin/DRSC_orthologs.pl
- Human gene records: UniProt (IDs above) and NCBI Gene

## Notes

The database-search scripts were written with AI assistance.
