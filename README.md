# TF-TG correlation across biological domains (master thesis)

<img width="827" height="900" alt="Flow_chart_2" src="https://github.com/user-attachments/assets/c363c726-57be-426a-951c-1bf595fcf7df" />


## Baits_bHLH:


Creating baits files for the bHLH-Annotator for the domains bacteria, fungi and animals. 
```
Usage:

python3 bait_files_bacteria.py 
--pep_file <PEP FILE FOR BAITS> 
--cds_file <CDS FILE FOR BAITS> 
--translation_file <TRANSLATION FILE> 
--out_path <PATH FOR RESULT FILES> 
--out_prefix <PREFIX FOR BAIT FILES>
```

```
Usage:

python3 bait_files_animals.py 
--TF_dir <DIRECTORY CONTAINING MODEL ORGANISM_TF_protein.fasta AND TFs_MODEL ORGANISM>
```

BaitsInfo in Excel

Baits files for fungi in Excel



## Header_correction:



For the bHLH-Annotator no '|' and '\_' -> splits the ID at '|' and replaces '\_' with '.'

```
Usage:
 
python3 header_correction_final.py 
--pep_dir <DIRECTORY WITH PEP FILES> 
--out_dir <PATH TO DIRECTORY FOR RESULTS>
```

## TF_annotation



#### Bacteria, fungi, animals:

#### bHLH-Annotator



example:

```
nohup bash -c '
# ---- Configurations ----
PEP\_DIR="/vol/data/results/bacteria/bHLH_pep"
OUT\_DIR="/vol/data/results/bacteria/bHLH_out_3"
LOG\_DIR="/vol/data"
BAITS="./data_bacteria/EcoliTFs_baits.fasta"
BAITSINFO="./data_bacteria/EcoliBaitsInfo.txt"
OPTIMISEDBAITS="./data_bacteria/EcoliTFs_baits.fasta"
REFERENCE="./data_bacteria/EcoliRefTFs_2.txt"
ATH="./data_plants/bHLH_at.fasta"
MODE_ALN="mafft"
CPU=16

for f in "$PEP\_DIR"/\*.fasta; do
  filename=$(basename "$f")
  species="${filename%.pep.fasta}"

  echo "Starte Annotation für $species ..."

  nohup python3 bHLH_annotator.py \
    --subject "$f" \
    --out "$OUT_DIR/$species" \
    --baits "$BAITS" \
    --baitsinfo "$BAITSINFO" \
    --optimisedbaits "$OPTIMISEDBAITS" \
    --reference "$REFERENCE" \
    --ath "$ATH" \
    --mode_aln "$MODE_ALN" \
    --cpu "$CPU" \
     > "$LOG_DIR/${species}_bHLH.out" 2>&1
done

echo "Alle Jobs abgeschlossen."
' > /vol/data/bHLH_annotator_bacteria.out 2>&1 &
```

#### plants:

#### bHLH-Annotator



```
Usage: 

python3 extract_pep_to_batches_final.py 
--folder_path <PATH TO PLANT FOLDER>
OPTIONAL:
--batch_out <BATCH OUTPUT FOLDER; Default: pep_batch>
--batch_size <BATCH SIZE; Default: 50>
```

example:

```
./run_bHLH_loop.sh 2
```
--> the number after ./run_bHLH_loop.sh correlates with the batch number


#### MYB-Annotator



example:

```
python3 ./MYB_annotator/MYB_annotator.py 
--baits ./plants/MYB_files/MYB_baits.fasta 
--info ./plants/MYB_files/MYB_baits.txt 
--out ./plants/myb_data_2/ 
--subjectdir ./plants/myb_pep_batch_2
```

## Orthologue_matrix



different scripts for bHLH- and MYB-Anntoator Output, because of different folder structures



#### bHLH-files:



Creates a transcription factor orthologue matrix with the header:

RefMember	Count/Total	Percentage [%]	species1	species2	…	speciesN

and optinal a presence/absence plot (heatmap)



```
Usage:

python3 ./scripts/bHLH_ortholog_matrix_2.py 
--in <DIRECTORY WITH bHLH RESULTS> 
--ref <BAITS REFERENCE FILE; RefTFs> 
--out <OUTPUT_DIRECTORY>
OPTIONAL:
--heatmap <NO_INPUT_NEEDED>
```

optional filtering species for number of TF orthologues found:

```
Usage:

python3 filter_num_myb.py 
--in <DIRECTORY_WITH_bHLH_RESULTS> 
--out <OUTPUT_DIRECTORY> 
OPTIONAL (one of the following needed)
--lower_limit <LOWER_LIMIT_NUMBER_bHLHs> if not set -> 1
--upper_limit <UPPER_LIMIT_NUMBER_bHLHs> if not set -> infinity
```

#### MYB-files:



Creates a transcription factor orthologue matrix with the header:

RefMYB	Count/Total	Percentage [%]	species1	species2	…	speciesN

and optinal a presence/absence plot (heatmap)

```
Usage:

python3 myb_ortholog_matrix.py 
--in <DIRECTORY_WITH_MYB_RESULTS> 
--ref <BAITS REFERENCE FILE; RefTFs> 
--out <OUTPUT_DIRECTORY>
OPTIONAL
--heatmap <NO_INPUT_NEEDED>
```

optional filtering species for number of TF orthologues found:

```
Usage:

python3 filter_num_myb.py 
--in <DIRECTORY_WITH_MYB_RESULTS> 
--out <OUTPUT_DIRECTORY>
OPTIONAL (one of the following needed)
--lower_limit <LOWER_LIMIT_NUMBER_MYBS> if not set -> 1
--upper_limit <UPPER_LIMIT_NUMBER_MYBS> if not set -> infinity
```

## Filtering_TF_orthologues



Filtering step 1 -> TF occurrence across the species

```
Usage:

python3 find_mybs_final.py 
--TF_orthologue_list <LIST WITH SPECIES AND NUMBER OF TF ORTHOLOGUES FOUND> 
--TF_orthologue_matrix_file <TF ORTHOLOGUE MATRIX FILE FROM bHLH OR MYB> 
--out_dir <OUTPUT DIRECTORY> 
--percentage_cutoff <TF OCCURENCE ACROSS SPECIES>
```

Creates a presence/absence plot phylogenetically sorted and coloured regarding family or order, depending on spieces_classification_list with the header:

species	family

OR

species order

```
Usage:

python3 sorte_heatmap_phylogeny_tree_cropped_final.py 
--species_tree <SPECIES TREE FILE FROM ORTHOFINDER> 
--csv_heatmap <CSV FROM INITIAL HEATMAP> 
--species_classification_list <SPECIES ORDER OR FAMILY LIST> 
--out_file <OUTPUT FILE> 
--classification_type <family OR order; HEADER FROM species_classification_list APART FROM species>
```

## List_TF_TG_pairs

TF-TG list plants created in Excel!



Bacteria, fungi, animals:

For the Downstream analysis, the original IDs are necessary, therefore the header correction for the bHLH-Annotator must be reversed

```
Usage:

python3 header_correction_back_final.py 
--bHLH_dir <DIRECTORY OF BHLH ANNOTATOR>
```

saved in the same folder as '04b_candidates_group_around_reference.txt'

original file saved in the same folder as '04b_candidates_group_around_reference_original.txt'





**Animals:**

The IDs from the animal data does not have a uniform prefix. To find the IDs from the bHLH-Annotator, the IDs in OrthoFinder must be corrected to uniform structure

```
Usage:

python3 header_correction_OrthoFinder_final.py 
--OrthoFinder_file_model <OrthoFinder OUT FILE OF MODEL ORGANISM> 
--out <OUT FILE FOR CORRECTED OrthoFinder OUT FILE OF MODEL ORGANISM>
```	

Filtering TF-TG pairs for evidence (only in animals necessary, because the file is too large for Excel)

```
Usage:	

python3 filter_evidence_pairs_final.py
--TFLink_file <FULL TFLink FILE FOR MODEL ORGANISM> 
--out <OUT FILE FILTERED>
```

Finding the TGs for the filtered TFs (only in animals necessary, because the file is too large for Excel)

```
Usage:

python3 find_TGs_matching_TFs_final.py 
--TFLink_file <TFLink file -> TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0_altered.tsv -> AFTER filter_evidence_pairs_final.py> 
--TF_list <TF LIST FROM find_mybs_final.py> 
--out_file <OUT FILE>
```



**Bacteria, fungi, animals:**

To create a list with TFs and their corresponding TG_IDs and TGs

TG_list_unique -> created in Excel

```
Usage:

python3 gene_from_header_cds_final.py 
--TG_list_unique <LIST WITH UNIQUE TG> 
--cds_model <CDS FILE FROM MODEL ORGANSIM BEFORE CLEAN> 
--out_file <PATH TO OUT FILE>
```

TF_TG_list created in Excel for Bacteria and fungi or for animals from find_TGs_matching_TFs_final.py

```
Usage:

python3 transcript_to_gene_final.py 
--TF_TG_list <LIST WITH HEADER Transcription factor	Transcript ID	Target gene -> Transcript ID COLUMN EMPTY> 
--gene_from_header <OUT FILE FROM gene_from_header_cds_final.py> 
--out <OUTPUT FILE>
```

## Spearman_correlation



Calculation of Spearman correlation coefficients.

species file: One species per row, same spelling as in RESULT file from myb_ortholog_matrix.py

Output:

{species}_compressed.tpms.txt: tpm count tables only containing TF and TG transcripts

{species}_spaerman_corr.txt: Spearman correlation for each TF-TG pair, with p-value correction

&#x09;-> header: TF_ref	TF_id	TG_ref	TG_id	spearman_r	p_value	n_samples	p_value_adj	significant

OPTIONAL, if --plot used:

Folder: with scatter_plots

```
Usage:

python3 spearman_correlation.py 
--orthologe_TF <DIRECTORY TO THE TF ORTHOLOGE FILES (e.g. MYB)> 
--species <FILE WITH SPECIES FOR ANALYSIS> 
--TF_TG_list <LIST WITH TF AND RELATED TG> 
--orthoFinder_file <SINGLE ORTHOFINDER RESULTS FILE FROM MODEL ORGANISM> 
--counts_dir <DIRECTORY TO COUNTTABLES> 
--ref <REFERENCE SPECIES NAME e.g. Arabidopsis_thaliana> 
--out_dir <OUTPUT DIRECTORY> 
--organism <plants OR fungi OR bacteria OR animals> 
OPTIONAL:
--orthologe_TF_2 <SECOND TF ORTHOLOGE DIRECTORY (e.g. bHLH)> 
--plot <OPTIONAL: creates Scatter-Plots for reference species; NUMBER OF SCATTER PLOTS> 
--p_correction <CORRECTION p-value DEFAULT: bonferroni; ALTERNATIVE: fdr_bh>
```

Creates bar plots for each species -> Spearman correlation coefficient per TF orthologue - TG orthogroup member pair

```
Usage:

python3 bar_plots.py 
--files <COMMA-SEPARATED LIST OF TSV/CSV FILES> 
--input-dir <OPTIONAL: FOLDER TO READ ALL MATCHING FILES FROM> 
--pattern <OPTIONAL: FILE PATTERN FOR 
--input-dir, DEFAULT: *.tsv> 
--outdir <OPTIONAL: OUTPUT FOLDER, DEFAULT: ./barplots>
--strong-corr-cutoff <ONLY KEEP TF_ref AND TG_ref (EACH A GENE) IF mean(|spearman_r|) OVER ALL THEIR TRANSCRIPTS AND CORRELATION PARTNERS IS >= CUTOFF> 
--variable-corr-cutoff <ONLY KEEP TF_ref AND TG_ref (EACH A GENE) IF THE RANGE (max - min) OF THEIR spearman_r VALUES (OVER ALL TRANSCRIPTS) IS >= CUTOFF> 
--corr-filter-combine <OPTIONAL: 'or' (DEFAULT) OR 'and' -- IF BOTH CUTOFFS ABOVE ARE SET, A TF_ref/TG_ref MUST FULFILL ONE (or) OR BOTH (and) CRITERIA> 
--pair-variable-cutoff <OPTIONAL: ONLY KEEP A SINGLE TF_ref-TG_ref PAIR IF THE RANGE (max - min) OF ITS OWN spearman_r VALUES (OVER ALL TF_id x TG_id TRANSCRIPT COMBINATIONS OF \*THAT PAIR ONLY\*) IS >= CUTOFF. UNLIKE --variable-corr-cutoff, THIS DOES NOT LOOK AT THE WHOLE TF_ref/TG_ref (ACROSS ALL ITS PARTNERS), ONLY AT THIS ONE PAIR. APPLIED \*AFTER\* --strong-corr-cutoff/--variable-corr-cutoff, WHICH CAN SERVE AS A CHEAP PRE-FILTER TO SHRINK THE DATA FIRST. DEFAULT: FEATURE OFF> 
--top-n-tfs <OPTIONAL: FROM THE FILTERED DATA OF ALL FILES, ONLY KEEP THE top-n-tfs TF_ref THAT OCCUR IN THE MOST SPECIES (= FILES); FOR THESE, ONE EXTRA BAR PLOT PER TF AND SPECIES IS CREATED IN <outdir>/top_tfs. DEFAULT: FEATURE OFF> 
--tg-batch-size <OPTIONAL: MAX AMOUNT OF TG_id (BARS ON THE X-AXIS) PER PLOT -- IF EXCEEDED, THE PLOT IS SPLIT INTO SEVERAL FILES (..._tgbatch1-of-N.png etc). ALSO APPLIED TO THE --top-n-tfs PLOTS. DEFAULT: NO LIMIT> 
--tf-batch-size <OPTIONAL: MAX AMOUNT OF TF_id (BARS PER TG_id GROUP) PER PLOT -- IF EXCEEDED, THE PLOT IS SPLIT INTO SEVERAL FILES (..._tfbatch1-of-N.png etc). DEFAULT: NO LIMIT>
```

Creates violin plots with spearman correlation coefficients for each TF-TG pair with coefficients from each species in one domain

Output:

10 violins per plot -> multiple violin plots

{domain}\_means.txt -> mean, median spearman correlation coefficient for each TF-TG pair, including absolute mean and absolute median; m: number of species containing pair; n: number TF-TG orthologue/orthogroup member pairs

&#x09;-> header: TF-TG pair	mean	median	n	m	mean\_abs	median\_abs

```
Usage:

python3 violin_plot_optimized_2.py 
--TF_TG_list <LIST WITH TF AND RELATED TG> 
--spearman_corr_dir <PATH TO SPEARMAN CORRELATION DIR> 
--out_dir <OUTPUT DIRECTORY> 
--organism <plants OR fungi OR bacteria OR animals>
OPTIONAL:
--n <NUMBER OF RANDOM SPECIES TO ANALYZE>
--sig-only <ONLY SIGNIFICANT TF-TG PAIRS, NO VALUE NEEDED>
--p-col <p_value_adj (DEFAULT) OR p_value, FOR --sig-only>
--alpha <SIGNIFICANCE THRESHOLD; DEFAULT: 0.05>
```

Creates violin plots with one violin per domain

-> mean_dir: one {domain}_mean.txt files from violin_plot_optimized_2.py per domain

Output:

-> one violin plot for each column: mean, median, mean\_abs, median\_abs

-> one violin plot: mean without self regulation TF-TG pairs

-> one violin plot: mean and mean without self regulation TF-TG pairs combined

```
Usage:

python3 violin_plot_domains.py 
--mean_dir <PATH TO DIR WITH MEAN FILES> 
--out_dir <OUTPUT DIRECTORY>
```


