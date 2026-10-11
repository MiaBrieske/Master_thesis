# TF annotation and TF–TG correlation pipeline

Scripts from my master thesis for identifying bHLH and MYB transcription factors (TFs) across species of different organism groups (plants, fungi, bacteria, animals) and for correlating the expression of TF orthologues with their target genes (TGs).

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Input data](#input-data)
- [1. Bait files](#1-bait-files)
- [2. Header correction](#2-header-correction)
- [3. TF annotation](#3-tf-annotation)
- [4. Orthologue matrix](#4-orthologue-matrix)
- [5. Filtering TF orthologues](#5-filtering-tf-orthologues)
- [6. TF–TG pair list](#6-tftg-pair-list)
- [7. Spearman correlation](#7-spearman-correlation)
- [8. Plots](#8-plots)

## Overview

<img width="827" height="900" alt="Flow_chart_2_github" src="https://github.com/user-attachments/assets/8ba23729-c6cc-4ffa-a547-1aabb41f3322" />

The bHLH-Annotator is used for all organism groups, the MYB-Annotator for plants only. Some steps are done in Excel; these are marked below.

## Requirements

**Python 3** with:

| Package | Used for |
|---|---|
| pandas, NumPy | Tables and arrays |
| SciPy | Spearman correlation, hierarchical clustering |
| statsmodels | p-value correction |
| Matplotlib, seaborn | Plots |
| Biopython | Reading the OrthoFinder species tree |

Install them, for example with conda:

```bash
conda create -n tf_pipeline -c conda-forge python=3.12 pandas numpy scipy statsmodels matplotlib seaborn biopython
conda activate tf_pipeline
```

**External tools:**

- [bHLH-Annotator](https://github.com/bpucker/bHLH_annotator)
- [MYB-Annotator](https://github.com/bpucker/MYB_annotator)
- [OrthoFinder](https://github.com/davidemms/OrthoFinder)
- MAFFT (alignment mode used by the bHLH-Annotator)

## Input data

One file of each type per species:

- **Pep files:** protein sequences (FASTA). Input for the bHLH-Annotator, the MYB-Annotator and OrthoFinder.
- **CDS files:** coding sequences (FASTA). Used to map transcripts to genes.
- **Count tables:** TPM expression values. Used for the Spearman correlation.

## 1. Bait files

The bHLH-Annotator needs bait sequences, a BaitsInfo file and a reference TF list from a model organism.

**Bacteria:**

```bash
python3 bait_files_bacteria.py \
    --pep_file <PEP FILE FOR BAITS> \
    --cds_file <CDS FILE FOR BAITS> \
    --translation_file <TRANSLATION FILE> \
    --out_path <PATH FOR RESULT FILES> \
    --out_prefix <PREFIX FOR BAIT FILES>
```

**Animals:**

```bash
python3 bait_files_animals.py \
    --TF_dir <DIRECTORY CONTAINING <MODEL ORGANISM>_TF_protein.fasta AND TFs_<MODEL ORGANISM>>
```

**Fungi:** bait files are created in Excel.

The BaitsInfo file is created in Excel for all organism groups.

## 2. Header correction

The bHLH-Annotator does not accept `|` and `_` in sequence IDs. This script splits each ID at `|` and replaces `_` with `.`.

```bash
python3 header_correction_final.py \
    --pep_dir <DIRECTORY WITH PEP FILES> \
    --out_dir <OUTPUT DIRECTORY>
```

The correction is reversed later for the downstream analysis (see [step 6](#6-tftg-pair-list)).

## 3. TF annotation

The TFs are annotated with external tools that are not part of this repository: the bHLH-Annotator for all organism groups and the MYB-Annotator for plants only (see [Requirements](#requirements)). How to run them is described in their own documentation.

Inputs from the previous steps:

- **bHLH-Annotator:** the header-corrected pep files from [step 2](#2-header-correction) and the bait files from [step 1](#1-bait-files).
- **MYB-Annotator:** the pep files.

The bHLH-Annotator output is expected with one result folder per species, named after the species. The scripts in [step 4](#4-orthologue-matrix) read these results.

### Optional: split pep files into batches

For large numbers of species (used here for the plants), the pep files can be split into batches that are annotated one after another:

```bash
python3 extract_pep_to_batches_final.py \
    --folder_path <PATH TO PEP FOLDER> \
    [--batch_out <BATCH OUTPUT FOLDER; default: pep_batch>] \
    [--batch_size <BATCH SIZE; default: 50>]
```

## 4. Orthologue matrix

There are separate scripts for the bHLH- and MYB-Annotator output because their folder structures differ.

### bHLH

Creates a TF orthologue matrix with the header

```
RefMember	Count/Total	Percentage [%]	species1	species2	…	speciesN
```

and optionally a presence/absence heatmap.

```bash
python3 ./scripts/bHLH_ortholog_matrix_2.py \
    --in <DIRECTORY WITH bHLH RESULTS> \
    --ref <BAITS REFERENCE FILE (RefTFs)> \
    --out <OUTPUT DIRECTORY> \
    [--heatmap]
```

### MYB

Creates a TF orthologue matrix with the header

```
RefMYB	Count/Total	Percentage [%]	species1	species2	…	speciesN
```

and optionally a presence/absence heatmap.

```bash
python3 myb_ortholog_matrix.py \
    --in <DIRECTORY WITH MYB RESULTS> \
    --ref <BAITS REFERENCE FILE (RefTFs)> \
    --out <OUTPUT DIRECTORY> \
    [--heatmap]
```

### Optional: filter species by number of TF orthologues

Works for both bHLH and MYB results. At least one of the two limits must be set.

```bash
python3 filter_num_myb.py \
    --in <DIRECTORY WITH bHLH OR MYB RESULTS> \
    --out <OUTPUT DIRECTORY> \
    [--lower_limit <MINIMUM NUMBER OF TFs; default: 1>] \
    [--upper_limit <MAXIMUM NUMBER OF TFs; default: no limit>]
```

## 5. Filtering TF orthologues

Keeps TFs that occur in at least the given percentage of species:

```bash
python3 find_mybs_final.py \
    --TF_orthologue_list <LIST WITH SPECIES AND NUMBER OF TF ORTHOLOGUES FOUND> \
    --TF_orthologue_matrix_file <TF ORTHOLOGUE MATRIX FROM bHLH OR MYB> \
    --out_dir <OUTPUT DIRECTORY> \
    --percentage_cutoff <TF OCCURRENCE ACROSS SPECIES>
```

Creates a presence/absence heatmap sorted by the OrthoFinder species tree and coloured by family or order. The classification list has the header `species	family` or `species	order`.

```bash
python3 sorte_heatmap_phylogeny_tree_cropped_final.py \
    --species_tree <SPECIES TREE FILE FROM ORTHOFINDER> \
    --csv_heatmap <CSV FROM INITIAL HEATMAP> \
    --species_classification_list <SPECIES ORDER OR FAMILY LIST> \
    --out_file <OUTPUT FILE> \
    --classification_type <family OR order>
```

## 6. TF–TG pair list

**Plants:** the TF–TG list is created in Excel.

### Bacteria, fungi, animals: reverse the header correction

The downstream analysis needs the original IDs:

```bash
python3 header_correction_back_final.py \
    --bHLH_dir <DIRECTORY OF bHLH-ANNOTATOR RESULTS>
```

The corrected file is saved as `04b_candidates_group_around_reference.txt`, the original as `04b_candidates_group_around_reference_original.txt`, both in the same folder.

### Animals only

The animal IDs do not have a uniform prefix. To match them with the bHLH-Annotator IDs, the OrthoFinder IDs are brought into a uniform structure:

```bash
python3 header_correction_OrthoFinder_final.py \
    --OrthoFinder_file_model <ORTHOFINDER OUTPUT FILE OF MODEL ORGANISM> \
    --out <CORRECTED OUTPUT FILE>
```

The TFLink file is too large for Excel, so filtering is done with scripts. First, keep only TF–TG pairs with evidence:

```bash
python3 filter_evidence_pairs_final.py \
    --TFLink_file <FULL TFLink FILE FOR MODEL ORGANISM> \
    --out <FILTERED OUTPUT FILE>
```

Then find the TGs of the filtered TFs:

```bash
python3 find_TGs_matching_TFs_final.py \
    --TFLink_file <TFLink FILE AFTER filter_evidence_pairs_final.py> \
    --TF_list <TF LIST FROM find_mybs_final.py> \
    --out_file <OUTPUT FILE>
```

Example TFLink file: `TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0_altered.tsv`.

### Bacteria, fungi, animals: map transcripts to genes

`TG_list_unique` is created in Excel.

```bash
python3 gene_from_header_cds_final.py \
    --TG_list_unique <LIST WITH UNIQUE TGs> \
    --cds_model <CDS FILE OF MODEL ORGANISM, BEFORE CLEANING> \
    --out_file <OUTPUT FILE>
```

The TF–TG list is created in Excel for bacteria and fungi, and comes from `find_TGs_matching_TFs_final.py` for animals. Its header is `Transcription factor	Transcript ID	Target gene`, with the `Transcript ID` column left empty; this script fills it:

```bash
python3 transcript_to_gene_final.py \
    --TF_TG_list <TF–TG LIST> \
    --gene_from_header <OUTPUT FROM gene_from_header_cds_final.py> \
    --out <OUTPUT FILE>
```

## 7. Spearman correlation

Calculates the Spearman correlation coefficient for each TF–TG pair in each species.

To check which correlation method fits the data, `correlation_scatter_plot.py` creates scatter plots of the expression data.

The species file contains one species per line, spelled as in the result file of `myb_ortholog_matrix.py`.

```bash
python3 spearman_correlation.py \
    --orthologe_TF <DIRECTORY WITH TF ORTHOLOGUE FILES (e.g. MYB)> \
    --species <FILE WITH SPECIES FOR ANALYSIS> \
    --TF_TG_list <TF–TG LIST> \
    --orthoFinder_file <SINGLE ORTHOFINDER RESULT FILE OF MODEL ORGANISM> \
    --counts_dir <DIRECTORY WITH COUNT TABLES> \
    --ref <REFERENCE SPECIES, e.g. Arabidopsis_thaliana> \
    --out_dir <OUTPUT DIRECTORY> \
    --organism <plants | fungi | bacteria | animals> \
    [--orthologe_TF_2 <SECOND TF ORTHOLOGUE DIRECTORY (e.g. bHLH)>] \
    [--plot <NUMBER OF SCATTER PLOTS FOR THE REFERENCE SPECIES>] \
    [--p_correction <bonferroni (default) | fdr_bh>]
```

**Output:**

- `{species}_compressed.tpms.txt`: count table with only the TF and TG transcripts
- `{species}_spaerman_corr.txt`: Spearman correlation for each TF–TG pair, with p-value correction. Header:
  ```
  TF_ref	TF_id	TG_ref	TG_id	spearman_r	p_value	n_samples	p_value_adj	significant
  ```
- With `--plot`: a folder with scatter plots

## 8. Plots

### Bar plots

One bar plot per species, showing the Spearman correlation coefficient for each TF orthologue – TG orthogroup member pair.

```bash
python3 bar_plots.py \
    --files <COMMA-SEPARATED LIST OF TSV/CSV FILES> \
    [--input-dir <FOLDER TO READ ALL MATCHING FILES FROM>] \
    [--pattern <FILE PATTERN FOR --input-dir; default: *.tsv>] \
    [--outdir <OUTPUT FOLDER; default: ./barplots>] \
    [--strong-corr-cutoff <CUTOFF>] \
    [--variable-corr-cutoff <CUTOFF>] \
    [--corr-filter-combine <or (default) | and>] \
    [--pair-variable-cutoff <CUTOFF>] \
    [--top-n-tfs <N>] \
    [--tg-batch-size <N>] \
    [--tf-batch-size <N>]
```

| Option | Effect |
|---|---|
| `--strong-corr-cutoff` | Keep a `TF_ref` or `TG_ref` (each a gene) only if the mean of \|spearman_r\| over all its transcripts and correlation partners is ≥ cutoff. |
| `--variable-corr-cutoff` | Keep a `TF_ref` or `TG_ref` only if the range (max − min) of its `spearman_r` values over all transcripts is ≥ cutoff. |
| `--corr-filter-combine` | If both cutoffs above are set: a gene must meet one (`or`, default) or both (`and`). |
| `--pair-variable-cutoff` | Keep a single `TF_ref`–`TG_ref` pair only if the range of its own `spearman_r` values (over all `TF_id` × `TG_id` transcript combinations of that pair) is ≥ cutoff. Applied after the two filters above, which can serve as a cheap pre-filter. Off by default. |
| `--top-n-tfs` | Keep only the N `TF_ref` that occur in the most species (files) and create one extra bar plot per TF and species in `<outdir>/top_tfs`. Off by default. |
| `--tg-batch-size` | Maximum number of `TG_id` (bars on the x-axis) per plot; larger plots are split into several files (`..._tgbatch1-of-N.png`). Also applies to the `--top-n-tfs` plots. No limit by default. |
| `--tf-batch-size` | Maximum number of `TF_id` (bars per `TG_id` group) per plot; larger plots are split (`..._tfbatch1-of-N.png`). No limit by default. |

### Violin plots per TF–TG pair

Violin plots of the Spearman correlation coefficients for each TF–TG pair, with the coefficients of all species in one domain.

```bash
python3 violin_plot_optimized_2.py \
    --TF_TG_list <TF–TG LIST> \
    --spearman_corr_dir <DIRECTORY WITH SPEARMAN CORRELATION FILES> \
    --out_dir <OUTPUT DIRECTORY> \
    --organism <plants | fungi | bacteria | animals> \
    [--n <NUMBER OF RANDOM SPECIES TO ANALYSE>] \
    [--sig-only] \
    [--p-col <p_value_adj (default) | p_value>] \
    [--alpha <SIGNIFICANCE THRESHOLD; default: 0.05>]
```

`--sig-only` keeps only significant TF–TG pairs; `--p-col` and `--alpha` set which p-value and threshold are used for this.

**Output:**

- Violin plots with 10 violins per plot (several files)
- `{domain}_means.txt`: mean and median Spearman correlation coefficient for each TF–TG pair, also as absolute values. `m` is the number of species containing the pair, `n` the number of TF–TG orthologue/orthogroup member pairs. Header:
  ```
  TF-TG pair	mean	median	n	m	mean_abs	median_abs
  ```

### Violin plots per domain

One violin per domain. `--mean_dir` contains one `{domain}_means.txt` file from `violin_plot_optimized_2.py` per domain.

```bash
python3 violin_plot_domains.py \
    --mean_dir <DIRECTORY WITH MEAN FILES> \
    --out_dir <OUTPUT DIRECTORY>
```

**Output:**

- One violin plot for each of the columns `mean`, `median`, `mean_abs` and `median_abs`
- One violin plot of `mean` without self-regulating TF–TG pairs
- One violin plot combining `mean` with and without self-regulating TF–TG pairs
