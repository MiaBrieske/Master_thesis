# TF annotation and TF–TG correlation pipeline

Scripts from my master thesis for identifying bHLH and MYB transcription factors (TFs) across species of four organism groups (plants, fungi, bacteria, animals) and for correlating the expression of TF orthologues with the expression of their target genes (TGs).

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Input data](#input-data)
- [Data sources](#data-sources)
- [Running the scripts](#running-the-scripts)
- [1. Bait files](#1-bait-files)
- [2. Header correction](#2-header-correction)
- [3. TF annotation](#3-tf-annotation)
- [4. Orthologue matrix](#4-orthologue-matrix)
- [5. Filtering TF orthologues](#5-filtering-tf-orthologues)
- [6. TF–TG pair list](#6-tftg-pair-list)
- [7. Choosing the correlation method](#7-choosing-the-correlation-method)
- [8. Spearman correlation](#8-spearman-correlation)
- [9. Plots](#9-plots)

## Overview

<img width="827" height="900" alt="Flow_chart_2_github" src="https://github.com/user-attachments/assets/fcc94b62-6869-4fee-93b5-544cf9a694cc" />

The bHLH-Annotator is used for all organism groups, the MYB-Annotator for plants only. bHLH-Annotator, MYB-Annotator and OrthoFinder are external tools and not part of this repository. Some steps are done in Excel; these are marked below.

## Requirements

**Python 3** with:

| Package | Used for |
|---|---|
| pandas, NumPy | Tables and arrays |
| SciPy | Spearman correlation, hierarchical clustering |
| statsmodels | Multiple-testing correction of p-values |
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

## Input data

One file of each type per species:

- **Pep files:** protein sequences (FASTA). Input for the bHLH-Annotator, the MYB-Annotator and OrthoFinder.
- **CDS files:** coding sequences (FASTA, NCBI headers with `[gene=...]` / `[locus_tag=...]`). Used to map gene names to transcript IDs.
- **Count tables:** TPM expression values. Used for the correlation.

## Data sources

| Organism group | TF list | TF–TG pairs |
|---|---|---|
| Plants | Reference sets of the bHLH- and MYB-Annotator | Compiled from the literature (references in the thesis) |
| Fungi | YEASTRACT+ | YEASTRACT+ |
| Bacteria | RegulonDB | RegulonDB |
| Animals | AnimalTFDB | TFLink |

The data are not included in this repository. Please cite the original resources when using them:

- **AnimalTFDB 4.0:** Shen, W.-K., Chen, S.-Y., Gan, Z.-Q., Zhang, Y.-Z., Yue, T., Chen, M.-M., Xue, Y., Hu, H., & Guo, A.-Y. (2023). AnimalTFDB 4.0: A comprehensive animal transcription factor database updated with variation and expression annotations. Nucleic Acids Research, 51(D1), D39–D45. [doi:10.1093/nar/gkac907](https://doi.org/10.1093/nar/gkac907)
- **TFLink:** Liska, O., Bohár, B., Hidas, A., Korcsmáros, T., Papp, B., Fazekas, D., & Ari, E. (2022). TFLink: An integrated gateway to access transcription factor–target gene interactions for multiple species. Database, 2022, baac083. [doi:10.1093/database/baac083](https://doi.org/10.1093/database/baac083)
- **YEASTRACT+:** Teixeira, M. C., Viana, R., Palma, M., Oliveira, J., Galocha, M., Mota, M. N., Couceiro, D., Pereira, M. G., Antunes, M., Costa, I. V., Pais, P., Parada, C., Chaouiya, C., Sá-Correia, I., & Monteiro, P. T. (2023). YEASTRACT+: A portal for the exploitation of global transcription regulation and metabolic model data in yeast biotechnology and pathogenesis. Nucleic Acids Research, 51(D1), D785–D791. [doi:10.1093/nar/gkac1041](https://doi.org/10.1093/nar/gkac1041)
- **RegulonDB v12.0:** Salgado, H., Gama-Castro, S., Lara, P., Mejia-Almonte, C., Alarcón-Carranza, G., López-Almazo, A. G., Betancourt-Figueroa, F., Peña-Loredo, P., Alquicira-Hernández, S., Ledezma-Tejeida, D., Arizmendi-Zagal, L., Mendez-Hernandez, F., Diaz-Gomez, A. K., Ochoa-Praxedis, E., Muñiz-Rascado, L. J., García-Sotelo, J. S., Flores-Gallegos, F. A., Gómez, L., Bonavides-Martínez, C., … Collado-Vides, J. (2024). RegulonDB v12.0: A comprehensive resource of transcriptional regulation in E. coli K-12. Nucleic Acids Research, 52(D1), D255–D264. [doi:10.1093/nar/gkad1072](https://doi.org/10.1093/nar/gkad1072)

**External tools:**

- **bHLH-Annotator:** Thoben, C., & Pucker, B. (2023). Automatic annotation of the bHLH gene family in plants. BMC Genomics, 24(1), 780. [doi:10.1186/s12864-023-09877-](https://doi.org/10.1186/s12864-023-09877-2)
- **MYB-Annotator:** Pucker, B. (2022). Automatic identification and annotation of MYB gene family members in plants. BMC Genomics, 23(1), 220. [doi:10.1186/s12864-022-08452-5](https://doi.org/10.1186/s12864-022-08452-5)
- **OrthoFinder:** Emms, D. M., Liu, Y., Belcher, L., Holmes, J., & Kelly, S. (2025). OrthoFinder: Scalable phylogenetic orthology inference for comparative genomics. Bioinformatics. [doi:10.1101/2025.07.15.664860](https://doi.org/10.1101/2025.07.15.664860)

## Running the scripts

- All scripts are called as `python3 <script>.py --option value ...`.
- Running a script without its required options prints its usage.
- Output folders are created automatically.
- All tables are tab-separated.

## 1. Bait files

The bHLH-Annotator needs a set of known TFs of a model organism: bait sequences, a BaitsInfo file and a reference TF list.

### Bacteria

```bash
python3 bait_files_bacteria.py \
    --pep_file <PEP FILE OF THE MODEL ORGANISM> \
    --cds_file <CDS FILE OF THE MODEL ORGANISM> \
    --translation_file <TRANSLATION FILE> \
    --out_path <OUTPUT DIRECTORY> \
    --out_prefix <PREFIX, e.g. Ecoli>
```

The translation file has a header line and one TF per line: TF name in column 1, gene name(s) in column 3 and locus tag(s) in column 4 (several separated by `|`). The script finds the matching CDS entries via `[gene=...]` or `[locus_tag=...]` and writes the corresponding protein sequences with the TF name as header.

**Output** (same content in all three files): `<prefix>TFs_baits.fasta`, `TFs_<prefix>.fasta`, `TFs_<prefix>.tsv`.

The BaitsInfo file and the reference TF list for bacteria are created in Excel.

### Animals

```bash
python3 bait_files_animals.py \
    --TF_dir <DIRECTORY WITH TF FILES> \
    [--out_dir <OUTPUT DIRECTORY; default: --TF_dir>]
```

`--TF_dir` contains two files per model organism:

- `TFs_<organism>.txt`: TF list with a header line; gene symbol in column 2 (column 3 if empty), TF family in column 4, protein IDs in column 5 (separated by `;`).
- `<organism>_TF_protein.fasta`: protein sequences with tab-separated headers that contain the gene symbol in field 2 or 3.

**Output** in `<out_dir>/bHLH_MYB/`:

| File | Content |
|---|---|
| `<organism>RefTFs.tsv` | Reference TF list |
| `<organism>BaitsInfo.txt` | BaitsInfo file |
| `<organism>TFs_baits.fasta` | Bait sequences |
| `TFs_<organism>.fasta` | Bait IDs only (headers) |

The characters `(`, `)` and `_` in gene symbols are replaced by `.`, because the bHLH-Annotator does not accept them.

### Fungi

Bait files are created in Excel.

## 2. Header correction

The bHLH-Annotator does not accept `|` and `_` in sequence IDs. For every pep file in `--pep_dir`, the script shortens the IDs and replaces `_` with `.`:

| Original header | Corrected header |
|---|---|
| `>lcl\|NZ_CP009072.1_prot_WP_000005.1_1` | `>NZ.CP009072.1.prot.WP.000005.1.1` |
| `>lcl_...` | everything after `lcl_`, with `_` → `.` |
| any other | `_` → `.` |

```bash
python3 header_correction_final.py \
    --pep_dir <DIRECTORY WITH PEP FILES> \
    --out_dir <OUTPUT DIRECTORY>
```

The corrected files keep their names. Duplicate IDs after the correction are printed. The correction is reversed later for the downstream analysis ([step 6](#6-tftg-pair-list)).

## 3. TF annotation

TFs are annotated with the bHLH-Annotator (all organism groups) and the MYB-Annotator (plants only). How to run them is described in their own documentation.

Inputs from the previous steps:

- **bHLH-Annotator:** header-corrected pep files ([step 2](#2-header-correction)) and bait files ([step 1](#1-bait-files)).
- **MYB-Annotator:** pep files.

The following steps expect one result folder per species:

| Tool | Folder name | Files read |
|---|---|---|
| bHLH-Annotator | `<species>` or `<species>.<anything>` | `RESULTS/04b_candidate_2_reference_mapping_file.txt`, `RESULTS/04b_candidates_group_around_reference.txt` |
| MYB-Annotator | `<5-digit number>_<species>` | `RESULTS/03b_*` |

### Optional: split pep files into batches

For large numbers of species (used here for the plants). The script reads all `*.tar.gz` archives in `--folder_path` (one per species), extracts the pep file from each archive and distributes them over the folders `<batch_out>_1`, `<batch_out>_2`, ... in the current directory. The batches can then be annotated one after another.

```bash
python3 extract_pep_to_batches_final.py \
    --folder_path <DIRECTORY WITH .tar.gz ARCHIVES> \
    [--batch_out <BATCH FOLDER PREFIX; default: pep_batch>] \
    [--batch_size <SPECIES PER BATCH; default: 50>]
```

## 4. Orthologue matrix

Collects the annotation results of all species in one matrix: one row per reference TF, one column per species, with the IDs of the TF orthologues found in that species (`-` if none).

There are separate scripts for bHLH and MYB because the annotators' folder structures differ.

### bHLH

```bash
python3 bHLH_ortholog_matrix_2.py \
    --in <DIRECTORY WITH bHLH-ANNOTATOR RESULTS> \
    --ref <REFERENCE TF LIST> \
    --out <OUTPUT PREFIX> \
    [--heatmap]
```

The reference TF names are taken from column 3 of `--ref` (column 1 if column 3 is empty).

### MYB

```bash
python3 myb_ortholog_matrix_2.py \
    --in <DIRECTORY WITH MYB-ANNOTATOR RESULTS> \
    --ref <REFERENCE MYB LIST> \
    --out <OUTPUT PREFIX> \
    [--heatmap]
```

The reference MYB names are taken from column 2 of `--ref`.

### Output (both scripts)

| File | Content |
|---|---|
| `<prefix>_RESULTS.txt` | Matrix sorted by occurrence. Header: `RefMember` (bHLH) or `RefMYB` (MYB), `Count/Total`, `Percentage [%]`, then one column per species |
| `<prefix>_heatmap.csv` | With `--heatmap`: presence/absence matrix (1/0), input for the sorted heatmap in [step 5](#5-filtering-tf-orthologues) |
| `<prefix>_heatmap.jpg` | With `--heatmap`: presence/absence heatmap |

### Number of TF orthologues per species

Counts for each species how many different reference TFs have an orthologue, and keeps only species within the given limits. The output list is the input for [step 5](#5-filtering-tf-orthologues).

```bash
python3 filter_num_bhlh.py \
    --in <DIRECTORY WITH bHLH-ANNOTATOR RESULTS> \
    --out <OUTPUT FILE> \
    [--lower_limit <MINIMUM; default: 0>] \
    [--upper_limit <MAXIMUM; default: no limit>]
```

```bash
python3 filter_num_myb.py \
    --in <DIRECTORY WITH MYB-ANNOTATOR RESULTS> \
    --out <OUTPUT FILE> \
    [--lower_limit <MINIMUM; default: 0>] \
    [--upper_limit <MAXIMUM; default: no limit>]
```

At least one of the two limits must be given. Output header: `Species	Number bHLHs` or `Species	Number MYBs`.

## 5. Filtering TF orthologues

Keeps only the reference TFs that have an orthologue in at least the given percentage of species:

```bash
python3 find_mybs_final.py \
    --TF_orthologue_list <SPECIES LIST FROM filter_num_bhlh.py / filter_num_myb.py> \
    --TF_orthologue_matrix_file <<prefix>_RESULTS.txt FROM STEP 4> \
    --out_dir <OUTPUT DIRECTORY> \
    --percentage_cutoff <PERCENTAGE, e.g. 80>
```

Only the species in `--TF_orthologue_list` are considered. Works for bHLH and MYB.

**Output:**

- `<percentage>%_found_TF_orthologues.txt`: filtered matrix (reference TF and the orthologue IDs per species)
- `<percentage>%_found_TF.txt`: list of the filtered TFs, input for [step 6](#6-tftg-pair-list)

Both files start with a comment line.

### Phylogenetically sorted heatmap

Draws the presence/absence heatmap with the species sorted by the OrthoFinder species tree, the tree above it, and the cells coloured by family or order.

```bash
python3 sorte_heatmap_phylogeny_tree_cropped_final.py \
    --species_tree <ORTHOFINDER SPECIES TREE (NEWICK)> \
    --csv_heatmap <<prefix>_heatmap.csv FROM STEP 4> \
    --species_classification_list <SPECIES CLASSIFICATION LIST> \
    --classification_type <family | order> \
    --out_file <OUTPUT FILE (.svg)> \
    [--species_filter_list <SPECIES TO KEEP, ONE PER LINE>] \
    [--TF_filter <TFs TO KEEP, e.g. <percentage>%_found_TF.txt>]
```

The classification list has the header `species	family` or `species	order`; `--classification_type` must match the second column. With `--TF_filter`, the first line of the file is skipped and the TFs are shown in the order of the file.

## 6. TF–TG pair list

The TF–TG pairs of the model organism, with the transcript IDs of the target genes. The final list has the header `Transcription factor	Transcript ID	Target gene` and is the input for [step 8](#8-spearman-correlation).

**Plants:** the TF–TG list is created in Excel.

### Bacteria, fungi, animals: reverse the header correction

The downstream analysis needs the original IDs. For every species folder, the script restores the original IDs in `RESULTS/04b_candidates_group_around_reference.txt`:

```bash
python3 header_correction_back_final.py \
    --bHLH_dir <DIRECTORY WITH bHLH-ANNOTATOR RESULTS>
```

The corrected file replaces `04b_candidates_group_around_reference.txt`; the original is kept as `04b_candidates_group_around_reference_original.txt`. Species that already have an `_original` file are skipped, so the script can safely be run again.

Only IDs in the NCBI format (`lcl|NZ_CP009072.1_prot_WP_000005.1_1`, at least 8 parts after the correction) can be restored; other IDs are reported as skipped and left out.

### Animals only

**OrthoFinder IDs.** The animal IDs do not have a uniform prefix. To match them with the bHLH-Annotator IDs, the first `_` in each ID of the OrthoFinder file of the model organism is replaced with `|`:

```bash
python3 header_correction_OrthoFinder_final.py \
    --OrthoFinder_file_model <ORTHOFINDER FILE OF MODEL ORGANISM> \
    --out <OUTPUT FILE>
```

**TF–TG pairs from TFLink.** The TFLink file is too large for Excel. First, keep only interactions with small-scale evidence (column 11 = `Yes`):

```bash
python3 filter_evidence_pairs_final.py \
    --TFLink_file <TFLink FILE OF MODEL ORGANISM> \
    --out <OUTPUT FILE>
```

Then keep the interactions of the filtered TFs from [step 5](#5-filtering-tf-orthologues):

```bash
python3 find_TGs_matching_TFs_final.py \
    --TFLink_file <OUTPUT OF filter_evidence_pairs_final.py> \
    --TF_list <<percentage>%_found_TF.txt FROM STEP 5> \
    --out_file <OUTPUT FILE>
```

Example TFLink file: `TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0_altered.tsv`. The output already has the header of the TF–TG list.

### Bacteria, fungi, animals: add transcript IDs

First, find the CDS IDs of the target genes:

```bash
python3 gene_from_header_cds_final.py \
    --TG_list_unique <LIST OF UNIQUE TGs> \
    --cds_model <CDS FILE OF MODEL ORGANISM, BEFORE HEADER CORRECTION> \
    --out_file <OUTPUT FILE>
```

`TG_list_unique` is created in Excel and has no header: either the gene name only, or locus tag(s) (separated by `;`) in column 1 and the gene name in column 2. Output: one line per match with `CDS ID	gene name / locus tag`.

Then fill the `Transcript ID` column of the TF–TG list:

```bash
python3 transcript_to_gene_final.py \
    --TF_TG_list <TF–TG LIST> \
    --gene_from_header <OUTPUT OF gene_from_header_cds_final.py> \
    --out <OUTPUT FILE>
```

The TF–TG list comes from Excel (bacteria, fungi) or from `find_TGs_matching_TFs_final.py` (animals). The target genes in column 3 are matched with the gene names from `gene_from_header_cds_final.py`.

## 7. Choosing the correlation method

Before the correlation is calculated, scatter plots of TF against TG expression in the reference species are used to decide between Pearson and Spearman correlation. Pearson measures linear relationships and assumes roughly normally distributed values without strong outliers. Spearman uses ranks, so it also captures monotonic but non-linear relationships and is robust against outliers and skewed values, as is typical for TPM data.

```bash
python3 correlation_scatter_plot.py \
    --orthologe_TF <DIRECTORY TO THE TF ORTHOLOGE FILES (z.B. MYB)> \
    [--orthologe_TF_2 <OPTIONAL: ZWEITES TF ORTHOLOGE DIRECTORY, z.B. bHLH>] \
    --species <FILE WITH SPECIES FOR ANALYSIS> \
    --TF_TG_list <LIST WITH TF AND RELATED TG> \
    --orthoFinder_file <SINGLE ORTHOFINDER RESULTS FILE> \
    --counts_dir <DIRECTORY TO COUNTTABLES> \
    --ref <REFERENCE SPECIES NAME e.g. Arabidopsis_thaliana> \
    --out_dir <OUTPUT DIRECTORY> \
    [--plot] <OPTIONAL: erstellt zusätzlich Scatter-Plots für die Referenz-Spezies>
```

## 8. Spearman correlation

For every species, the script collects the TPM values of the TF orthologues and TG orthologues, and calculates the Spearman correlation coefficient for each combination of TF transcript and TG transcript.

```bash
python3 spearman_correlation.py \
    --orthologe_TF <DIRECTORY WITH TF-ANNOTATOR RESULTS (e.g. MYB)> \
    --species <SPECIES FILE> \
    --TF_TG_list <TF–TG LIST FROM STEP 6> \
    --orthoFinder_file <ORTHOFINDER FILE OF MODEL ORGANISM> \
    --counts_dir <DIRECTORY WITH COUNT TABLES> \
    --ref <REFERENCE SPECIES, e.g. Arabidopsis_thaliana> \
    --out_dir <OUTPUT DIRECTORY> \
    --organism <plants | fungi | bacteria | animals> \
    [--orthologe_TF_2 <SECOND TF-ANNOTATOR RESULT DIRECTORY (e.g. bHLH)>] \
    [--p_correction <fdr_bh (default) | bonferroni>] \
    [--plot <NUMBER OF TF–TG PAIRS>]
```

**Inputs:**

- `--species`: one species per line, spelled like the species columns of the orthologue matrix.
- `--TF_TG_list`: header line; TF in column 1, TG transcript ID in column 2, TG name in column 3 (if empty, the ID is used). For plants, TG IDs are converted to upper case.
- `--orthologe_TF`: annotator result folders as described in [step 3](#3-tf-annotation); the TF orthologues are read from `RESULTS/*group_around*`.
- `--orthoFinder_file`: species name in column 2, model organism IDs in column 3, orthologue IDs in column 4 (several separated by `, `). Used to find the TG orthologues in all species except the reference species.
- `--counts_dir`: one folder per species with `<species>.tpms.txt` or `<species>.tpms.tsv`; first column transcript ID, then one column per sample.

Pairs with fewer than 5 samples with values are skipped. p-values are corrected for multiple testing per species (Benjamini–Hochberg FDR by default).

**Output** in `--out_dir`:

- `<species>_compressed.tpms.txt`: count table with only the TF and TG transcripts (columns `RefMember`, `Orthologue`, samples)
- `<species>_spearman_corr.txt`: one line per combination of TF transcript and TG transcript. Header:
  ```
  TF_ref	TF_id	TG_ref	TG_id	spearman_r	p_value	n_samples	p_value_adj	significant
  ```
- With `--plot <n>`: scatter plots of the first n TF–TG pairs of the reference species in `scatter_plots/`

## 9. Plots

### Bar plots

One bar plot per species: one group of bars per TG transcript, one bar per TF transcript, bar height = Spearman correlation coefficient. TG transcripts are grouped by gene and ordered by hierarchical clustering.

```bash
python3 bar_plots.py \
    --input-dir <DIRECTORY WITH _spearman_corr.txt FILES> \
    --pattern "*_spearman_corr.txt" \
    --strong-corr-cutoff <CUTOFF> \
    [--outdir <OUTPUT DIRECTORY; default: ./barplots>]
```

Instead of `--input-dir`, single files can be given with `--files <FILE1,FILE2,...>`. The default `--pattern` is `*.tsv`, so set it when reading the output of [step 8](#8-spearman-correlation).

At least one of the two correlation filters must be set. The statistics for these filters are calculated from significant values only (`p_value_adj` ≤ 0.05).

| Option | Effect |
|---|---|
| `--strong-corr-cutoff` | Keep a TF or TG gene only if the mean of \|spearman_r\| over all its transcripts and correlation partners is ≥ cutoff. |
| `--variable-corr-cutoff` | Keep a TF or TG gene only if the range (max − min) of its `spearman_r` values is ≥ cutoff. |
| `--corr-filter-combine` | If both cutoffs are set: a gene must meet one (`or`, default) or both (`and`). |
| `--pair-variable-cutoff` | Keep a single TF–TG pair only if the range of its own `spearman_r` values (over all its transcript combinations) is ≥ cutoff. Applied after the filters above. Off by default. |
| `--top-n-tfs` | Keep only the N TFs that occur in the most species, and create one extra bar plot per TF and species in `<outdir>/top_tfs`. Off by default. |
| `--tg-batch-size` | Maximum number of TG transcripts per plot; larger plots are split into several files (`..._tgbatch1-of-N.png`). No limit by default. |
| `--tf-batch-size` | Maximum number of TF transcripts per plot; larger plots are split (`..._tfbatch1-of-N.png`). No limit by default. |

**Output:** `<outdir>/<file name>_barplot.png` per input file.

### Violin plots per TF–TG pair

One violin per TF–TG pair, containing the correlation coefficients of all its transcript combinations in all species of one organism group. Combinations of a transcript with itself are excluded.

```bash
python3 violin_plot_optimized_2.py \
    --TF_TG_list <TF–TG LIST FROM STEP 6> \
    --spearman_corr_dir <OUTPUT DIRECTORY OF STEP 8> \
    --out_dir <OUTPUT DIRECTORY> \
    --organism <plants | fungi | bacteria | animals> \
    [--n <NUMBER OF RANDOMLY CHOSEN SPECIES>] \
    [--sig-only] \
    [--p-col <p_value_adj (default) | p_value>] \
    [--alpha <SIGNIFICANCE THRESHOLD; default: 0.05>] \
    [--cropped <LOWER Y-AXIS LIMIT; default: -1>]
```

`--sig-only` keeps only significant values; `--p-col` and `--alpha` set which p-value and threshold are used.

**Output:**

- `violinplot_001.png`, `violinplot_002.png`, ...: 10 TF–TG pairs per plot; above each violin, `n` (number of values) and the number of species
- `<organism>_means.tsv`: statistics per TF–TG pair. Header:
  ```
  TF-TG pair	mean	median	n	m	mean_abs	median_abs
  ```
  `n` = number of correlation values, `m` = number of species containing the pair, `_abs` = calculated from absolute values.

### Violin plots per domain

Compares the organism groups: one violin per group, built from the per-pair statistics of the previous step. Only TF–TG pairs found in more than one species are used.

```bash
python3 violin_plot_domains.py \
    --mean_dir <DIRECTORY WITH <organism>_means.tsv FILES> \
    --out_dir <OUTPUT DIRECTORY> \
    [--cropped <LOWER Y-AXIS LIMIT; default: -1>]
```

The organism group is taken from the file name (part before the first `_`).

**Output:**

- `violinplot_comparison_<column>.png` for `Mean`, `Median`, `Mean_abs` and `Median_abs`
- `violinplot_comparison_Mean_without_self_regulation.png`: mean without pairs in which TF and TG are the same gene
- `violinplot_comparison_Mean_combined.png`: both mean plots side by side
