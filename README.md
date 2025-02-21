# Pipeline for the Creation of an Orthologous Gene Database

## 1) Genome Annotation

Three different programs were used to annotate the genomes:

### BUSCO (Benchmarking Universal Single-Copy Orthologs)

Used to assess the integrity of genome assemblies and for gene prediction using **AUGUSTUS**. Detailed operation:

-   **BUSCO database**: Uses a library of single-copy orthologs, conserved and present in a unique form in most of the genomes of a taxonomic group. The library adapted to your group must be defined (`GROUPE_odb10`).
-   **Alignment and ortholog search**: BUSCO compares the sequences provided with the ortholog database using tools such as **BLAST+** or **HMMER**.
-   **Classification of orthologs** :
    -   *Complete*: orthologs found in their entirety.
    -   *Complete and duplicated*: orthologs present in several copies.
    -   *Fragmented*: orthologs partially found.
    -   *Missing*: orthologs not found.
-   **Integrity assessment**: Generation of reports showing the percentage of complete, duplicated, fragmented and missing genes.

BUSCO generates one FASTA file per gene. All the FASTA files for the complete classified genes (`single_copy_busco_sequences`) are compiled into a single FASTA file for each species.

### Miniprot

Used to annotate genomes based on protein sequences. Detailed operation:

-   **Protein sequence pre-processing**: Miniprot uses well-annotated protein sequences for mapping.
-   **Protein mapping**: The sequences are aligned with the target genome to identify the corresponding regions.
-   **Exon and intron identification**: Miniprot identifies exons (coding) and introns (non-coding).
-   **Annotation transfer**: Miniprot transfers protein annotations to the target genome.
-   **Report generation** : Results are produced in GFF format.

### Scipio

Scipio is used to identify coding genes and transfer functional annotations. Detailed operation:

-   **Protein alignment**: Scipio aligns protein sequences to the target genome.
-   **Exon and intron identification**: Identification of exons and introns in aligned regions.
-   **Gene annotation**: Scipio transfers functional annotations of proteins to the target genome.
-   **Results**: Scipio produces files in YAML format with gene annotations.

The genomes are annotated in three phases, with each phase using protein sequences from a different species. Each genome is therefore annotated three times.

The YAML files from Scipio are converted into GFF, and then several steps are carried out to process and clean up the results of the Miniprot and Scipio annotations:

### Processing GFF files

1.  **Copy** of the three GFF files.
2.  **Extraction of CDS** from GFF files (`listAlternativeTranscript2.py`).
3.  **Cleaning of GFF files** to eliminate alternative transcripts (`ExcludeCDSFromGFF.py`).
4.  **Reconciliation of cleaned CDS GFF** files.
5.  **Conversion of GFF files to FASTA** to obtain exonic sequences (`putIDinNameGFF.py`).
6.  **Concatenation** of exons into full transcripts (`concateExonFromBedtoolsGetFasta2_forCat.py`).

### Gene Selection and Translation

1.  **Gene sorting** to eliminate overlapping exons (`SelectOverlappingExonGFFsort_forCat.py`).
2.  **Recovery** of single transcripts (`copySeqList.py`).
3.  **Translation of transcripts** into protein sequences.
4.  **List of redundant** transcripts (`RenameTranscriptRedundants_forList.py`).
5.  **Renaming of redundant** transcripts (`RenameTranscriptRedundants_forFasta.py`).
6.  **Selection of final** transcripts (`SelectSeq.py`).
7.  **Transcript renaming** with species-specific information.

Bash scripts have been created to automate these steps for the outputs of **Miniprot** and **Scipio** : - `Script_concatenateGFF_Miniprot.bash` - `Script_concatenateGFF_Scipio.bash`

After annotation, we processed the results independently for each set of genes produced by **BUSCO**, **Miniprot** and **Scipio**.

------------------------------------------------------------------------

## 2)  Search for Orthologous Genes
