# TF-TG correlation across biological domains (master thesis)
### Baits\_bHLH:



python3 bait\_files\_bacteria.py \\
\--pep\_file <PEP FILE FOR BAITS> \\
\--cds\_file <CDS FILE FOR BAITS> \\
\--translation\_file <TRANSLATION FILE> \\
\--out\_path <PATH FOR RESULT FILES> \\
\--out\_prefix <PREFIX FOR BAIT FILES> 



python3 bait\_files\_animals.py \\
\--TF\_dir <DIRECTORY CONTAINING MODEL ORGANISM\_TF\_protein.fasta AND TFs\_MODEL ORGANISM>



BaitsInfo in Excel
Baits files for fungi in Excel 



### Header\_correction: 



For the bHLH-Annotator no '|' and '\_' -> splits the ID at '|' and replaces '\_' with '.'



python3 header\_correction\_final.py \\
\--pep\_dir <DIRECTORY WITH PEP FILES> \\
\--out\_dir <PATH TO DIRECTORY FOR RESULTS>



### TF\_annotation



#### Bacteria, fungi, animals:

#### bHLH-Annotator



nohup bash -c '

\# ---- Configurations ----
PEP\_DIR="/vol/data/results/bacteria/bHLH\_pep"
OUT\_DIR="/vol/data/results/bacteria/bHLH\_out\_3"
LOG\_DIR="/vol/data"
BAITS="./data\_bacteria/EcoliTFs\_baits.fasta"
BAITSINFO="./data\_bacteria/EcoliBaitsInfo.txt"
OPTIMISEDBAITS="./data\_bacteria/EcoliTFs\_baits.fasta"
REFERENCE="./data\_bacteria/EcoliRefTFs\_2.txt"
ATH="./data\_plants/bHLH\_at.fasta"
MODE\_ALN="mafft"
CPU=16
for f in "$PEP\_DIR"/\*.fasta; do
&#x20;   filename=$(basename "$f")
&#x20;   species="${filename%.pep.fasta}"

&#x20;   echo "Starte Annotation für $species ..."

&#x20;   nohup python3 bHLH\_annotator.py \\
&#x20;       --subject "$f" \\
&#x20;       --out "$OUT\_DIR/$species" \\
&#x20;       --baits "$BAITS" \\
&#x20;       --baitsinfo "$BAITSINFO" \\
&#x20;       --optimisedbaits "$OPTIMISEDBAITS" \\
&#x20;       --reference "$REFERENCE" \\
&#x20;       --ath "$ATH" \\
&#x20;       --mode\_aln "$MODE\_ALN" \\
&#x20;       --cpu "$CPU" \\
&#x20;       > "$LOG\_DIR/${species}\_bHLH.out" 2>\&1
done

echo "Alle Jobs abgeschlossen."
' > /vol/data/bHLH\_annotator\_bacteria.out 2>\&1 \&



#### plants:

#### bHLH-Annotator



python3 extract\_pep\_to\_batches\_final.py \\
\--folder\_path <PATH TO PLANT FOLDER> 
OPTIONAL:
\--batch\_out <BATCH OUTPUT FOLDER; Default: pep\_batch>
\--batch\_size <BATCH SIZE; Default: 50>

./run\_bHLH\_loop.sh 2 
\--> the number after ./run\_bHLH\_loop.sh correlates with the batch number



#### MYB-Annotator



python3 ./MYB\_annotator/MYB\_annotator.py \\
\--baits ./plants/MYB\_files/MYB\_baits.fasta \\
\--info ./plants/MYB\_files/MYB\_baits.txt \\
\--out ./plants/myb\_data\_2/ \\
\--subjectdir ./plants/myb\_pep\_batch\_2


### Orthologue\_matrix


different scripts for bHLH- and MYB-Anntoator Output, because of different folder structures


#### bHLH-files:


Creates a transcription factor orthologue matrix with the header:
RefMember	Count/Total	Percentage \[%]	species1	species2	…	speciesN
and optinal a presence/absence plot (heatmap)


python3 ./scripts/bHLH\_ortholog\_matrix\_2.py \\
\--in <DIRECTORY WITH bHLH RESULTS> \\
\--ref <BAITS REFERENCE FILE; RefTFs> \\
\--out <OUTPUT\_DIRECTORY> 
OPTIONAL:
\--heatmap <NO\_INPUT\_NEEDED>


optional filtering species for number of TF orthologues found:


python3 filter\_num\_myb.py \\
\--in <DIRECTORY\_WITH\_bHLH\_RESULTS> \\
\--out <OUTPUT\_DIRECTORY> \\
OPTIONAL (one of the following needed)
\--lower\_limit <LOWER\_LIMIT\_NUMBER\_bHLHs> if not set -> 1
\--upper\_limit <UPPER\_LIMIT\_NUMBER\_bHLHs> if not set -> infinity


#### MYB-files:


Creates a transcription factor orthologue matrix with the header:
RefMYB	Count/Total	Percentage \[%]	species1	species2	…	speciesN
and optinal a presence/absence plot (heatmap)


python3 myb\_ortholog\_matrix.py \\
\--in <DIRECTORY\_WITH\_MYB\_RESULTS> \\
\--ref <BAITS REFERENCE FILE; RefTFs> \\
\--out <OUTPUT\_DIRECTORY>
OPTIONAL
\--heatmap <NO\_INPUT\_NEEDED>


optional filtering species for number of TF orthologues found:


python3 filter\_num\_myb.py \\
\--in <DIRECTORY\_WITH\_MYB\_RESULTS> \\
\--out <OUTPUT\_DIRECTORY>
OPTIONAL (one of the following needed)
\--lower\_limit <LOWER\_LIMIT\_NUMBER\_MYBS> if not set -> 1
\--upper\_limit <UPPER\_LIMIT\_NUMBER\_MYBS> if not set -> infinity


### Filtering\_TF\_orthologues


Filtering step 1 -> TF occurrence across the species


python3 find\_mybs\_final.py \\
\--TF\_orthologue\_list <LIST WITH SPECIES AND NUMBER OF TF ORTHOLOGUES FOUND> \\
\--TF\_orthologue\_matrix\_file <TF ORTHOLOGUE MATRIX FILE FROM bHLH OR MYB> \\
\--out\_dir <OUTPUT DIRECTORY> \\
\--percentage\_cutoff <TF OCCURENCE ACROSS SPECIES>


Creates a presence/absence plot phylogenetically sorted and coloured regarding family or order, depending on spieces\_classification\_list with the header:
species	family
OR
species order


python3 sorte\_heatmap\_phylogeny\_tree\_cropped\_final.py \\
\--species\_tree <SPECIES TREE FILE FROM ORTHOFINDER> \\
\--csv\_heatmap <CSV FROM INITIAL HEATMAP> \\
\--species\_classification\_list <SPECIES ORDER OR FAMILY LIST> \\
\--out\_file <OUTPUT FILE> \\
\--classification\_type <family OR order; HEADER FROM species\_classification\_list APART FROM species>


### List\_TF\_TG\_pairs 
TF-TG list plants created in Excel!


Bacteria, fungi, animals:
For the Downstream analysis, the original IDs are necessary, therefore the header correction for the bHLH-Annotator must be reversed


python3 header\_correction\_back\_final.py \\
\--bHLH\_dir <DIRECTORY OF BHLH ANNOTATOR>
&#x09;				
saved in the same folder as '04b\_candidates\_group\_around\_reference.txt'
original file saved in the same folder as '04b\_candidates\_group\_around\_reference\_original.txt'

\------------------------------------------------------------------------------------------------------------------------------------------------------------------------
animals:
The IDs from the animal data does not have a uniform prefix. To find the IDs from the bHLH-Annotator, the IDs in OrthoFinder must be corrected to uniform structure


python3 header\_correction\_OrthoFinder\_final.py \\
\--OrthoFinder\_file\_model <OrthoFinder OUT FILE OF MODEL ORGANISM> \\
\--out <OUT FILE FOR CORRECTED OrthoFinder OUT FILE OF MODEL ORGANISM>
&#x09;																			
Filtering TF-TG pairs for evidence (only in animals necessary, because the file is too large for Excel)									
&#x09;
python3 filter\_evidence\_pairs\_final.py \\
\--TFLink\_file <FULL TFLink FILE FOR MODEL ORGANISM> \\
\--out <OUT FILE FILTERED>


Finding the TGs for the filtered TFs (only in animals necessary, because the file is too large for Excel)


python3 find\_TGs\_matching\_TFs\_final.py \\
\--TFLink\_file <TFLink file -> TFLink\_Mus\_musculus\_interactions\_All\_simpleFormat\_v1.0\_altered.tsv -> AFTER filter\_evidence\_pairs\_final.py> \\
\--TF\_list <TF LIST FROM find\_mybs\_final.py> \\
\--out\_file <OUT FILE>


\------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Bacteria, fungi, animals:
To create a list with TFs and their corresponding TG\_IDs and TGs 
TG\_list\_unique -> created in Excel 


python3 gene\_from\_header\_cds\_final.py \\
\--TG\_list\_unique <LIST WITH UNIQUE TG> \\
\--cds\_model <CDS FILE FROM MODEL ORGANSIM BEFORE CLEAN> \\
\--out\_file <PATH TO OUT FILE>


TF\_TG\_list created in Excel for Bacteria and fungi or for animals from find\_TGs\_matching\_TFs\_final.py


python3 transcript\_to\_gene\_final.py \\
\--TF\_TG\_list <LIST WITH HEADER Transcription factor	Transcript ID	Target gene -> Transcript ID COLUMN EMPTY> \\
\--gene\_from\_header <OUT FILE FROM gene\_from\_header\_cds\_final.py> \\
\--out <OUTPUT FILE>


### Spearman\_correlation


Calculation of Spearman correlation coefficients.
species file: One species per row, same spelling as in RESULT file from myb\_ortholog\_matrix.py
Output:
{species}\_compressed.tpms.txt: tpm count tables only containing TF and TG transcripts
{species}\_spaerman\_corr.txt: Spearman correlation for each TF-TG pair, with p-value correction 
&#x09;-> header: TF\_ref	TF\_id	TG\_ref	TG\_id	spearman\_r	p\_value	n\_samples	p\_value\_adj	significant
OPTIONAL, if --plot used:
Folder: with scatter\_plots


python3 spearman\_correlation.py \\
\--orthologe\_TF <DIRECTORY TO THE TF ORTHOLOGE FILES (e.g. MYB)> \\
\--species <FILE WITH SPECIES FOR ANALYSIS> \\
\--TF\_TG\_list <LIST WITH TF AND RELATED TG> \\
\--orthoFinder\_file <SINGLE ORTHOFINDER RESULTS FILE FROM MODEL ORGANISM> \\
\--counts\_dir <DIRECTORY TO COUNTTABLES> \\
\--ref <REFERENCE SPECIES NAME e.g. Arabidopsis\_thaliana> \\
\--out\_dir <OUTPUT DIRECTORY> \\
\--organism <plants OR fungi OR bacteria OR animals> \\
OPTIONAL:
\--orthologe\_TF\_2 <SECOND TF ORTHOLOGE DIRECTORY (e.g. bHLH)> \\
\--plot <OPTIONAL: creates Scatter-Plots for reference species; NUMBER OF SCATTER PLOTS> \\
\--p\_correction <CORRECTION p-value DEFAULT: bonferroni; ALTERNATIVE: fdr\_bh>




Creates bar plots for each species -> Spearman correlation coefficient per TF orthologue - TG orthogroup member pair 


&#x20;python3 bar\_plots.py \\
\--files <COMMA-SEPARATED LIST OF TSV/CSV FILES> \\
\--input-dir <OPTIONAL: FOLDER TO READ ALL MATCHING FILES FROM> \\
\--pattern <OPTIONAL: FILE PATTERN FOR \\
\--input-dir, DEFAULT: \*.tsv> \\
\--outdir <OPTIONAL: OUTPUT FOLDER, DEFAULT: ./barplots> \\
\--strong-corr-cutoff <ONLY KEEP TF\_ref AND TG\_ref (EACH A GENE) IF mean(|spearman\_r|) OVER ALL THEIR TRANSCRIPTS AND CORRELATION PARTNERS IS >= CUTOFF> \\
\--variable-corr-cutoff <ONLY KEEP TF\_ref AND TG\_ref (EACH A GENE) IF THE RANGE (max - min) OF THEIR spearman\_r VALUES (OVER ALL TRANSCRIPTS) IS >= CUTOFF> \\
\--corr-filter-combine <OPTIONAL: 'or' (DEFAULT) OR 'and' -- IF BOTH CUTOFFS ABOVE ARE SET, A TF\_ref/TG\_ref MUST FULFILL ONE (or) OR BOTH (and) CRITERIA> \\
\--pair-variable-cutoff <OPTIONAL: ONLY KEEP A SINGLE TF\_ref-TG\_ref PAIR IF THE RANGE (max - min) OF ITS OWN spearman\_r VALUES (OVER ALL TF\_id x TG\_id TRANSCRIPT COMBINATIONS OF \*THAT PAIR ONLY\*) IS >= CUTOFF. UNLIKE --variable-corr-cutoff, THIS DOES NOT LOOK AT THE WHOLE TF\_ref/TG\_ref (ACROSS ALL ITS PARTNERS), ONLY AT THIS ONE PAIR. APPLIED \*AFTER\* --strong-corr-cutoff/--variable-corr-cutoff, WHICH CAN SERVE AS A CHEAP PRE-FILTER TO SHRINK THE DATA FIRST. DEFAULT: FEATURE OFF> \\
\--top-n-tfs <OPTIONAL: FROM THE FILTERED DATA OF ALL FILES, ONLY KEEP THE top-n-tfs TF\_ref THAT OCCUR IN THE MOST SPECIES (= FILES); FOR THESE, ONE EXTRA BAR PLOT PER TF AND SPECIES IS CREATED IN <outdir>/top\_tfs. DEFAULT: FEATURE OFF> \\
\--tg-batch-size <OPTIONAL: MAX AMOUNT OF TG\_id (BARS ON THE X-AXIS) PER PLOT -- IF EXCEEDED, THE PLOT IS SPLIT INTO SEVERAL FILES (...\_tgbatch1-of-N.png etc). ALSO APPLIED TO THE --top-n-tfs PLOTS. DEFAULT: NO LIMIT> \\
\--tf-batch-size <OPTIONAL: MAX AMOUNT OF TF\_id (BARS PER TG\_id GROUP) PER PLOT -- IF EXCEEDED, THE PLOT IS SPLIT INTO SEVERAL FILES (...\_tfbatch1-of-N.png etc). DEFAULT: NO LIMIT>


Creates violin plots with spearman correlation coefficients for each TF-TG pair with coefficients from each species in one domain
Output:
10 violins per plot -> multiple violin plots
{domain}\_means.txt -> mean, median spearman correlation coefficient for each TF-TG pair, including absolute mean and absolute median; m: number of species containing pair; n: number TF-TG orthologue/orthogroup member pairs 
&#x09;-> header: TF-TG pair	mean	median	n	m	mean\_abs	median\_abs


python3 violin\_plot\_optimized\_2.py \\
\--TF\_TG\_list <LIST WITH TF AND RELATED TG> \\
\--spearman\_corr\_dir <PATH TO SPEARMAN CORRELATION DIR> \\
\--out\_dir <OUTPUT DIRECTORY> \\
\--organism <plants OR fungi OR bacteria OR animals> 
OPTIONAL: 
\--n <NUMBER OF RANDOM SPECIES TO ANALYZE> 
\--sig-only <ONLY SIGNIFICANT TF-TG PAIRS, NO VALUE NEEDED> 
\--p-col <p\_value\_adj (DEFAULT) OR p\_value, FOR --sig-only> 
\--alpha <SIGNIFICANCE THRESHOLD; DEFAULT: 0.05>


Creates violin plots with one violin per domain
\-> mean\_dir: one {domain}\_mean.txt files from violin\_plot\_optimized\_2.py per domain
Output:
\-> one violin plot for each column: mean, median, mean\_abs, median\_abs
\-> one violin plot: mean without self regulation TF-TG pairs
\-> one violin plot: mean and mean without self regulation TF-TG pairs combined


python3 violin\_plot\_domains.py \\
\--mean\_dir <PATH TO DIR WITH MEAN FILES> \\
\--out\_dir <OUTPUT DIRECTORY>
