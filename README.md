# Pipeline pour la Création d'une Base de Données de Gènes Orthologues

## 1) Annotation des Génomes

Trois programmes différents ont été utilisés pour annoter les génomes :

### BUSCO (Benchmarking Universal Single-Copy Orthologs)

Utilisé pour évaluer l'intégrité des assemblages de génomes et pour la prédiction de gènes en utilisant **AUGUSTUS**. Fonctionnement détaillé :

-   **Base de données BUSCO** : Utilise une bibliothèque d'orthologues à copie unique, conservés et présents sous une forme unique dans la plupart des génomes d'un groupe taxonomique. La bibliothèque adaptée à votre groupe doit être définit (`GROUPE_odb10`).
-   **Alignement et recherche des orthologues** : BUSCO compare les séquences fournies à la base de données des orthologues avec des outils comme **BLAST+** ou **HMMER**.
-   **Classification des orthologues** :
    -   *Complets* : orthologues trouvés dans leur intégralité.
    -   *Complets et dupliqués* : orthologues présents en plusieurs copies.
    -   *Fragmentés* : orthologues trouvés partiellement.
    -   *Manquants* : orthologues non trouvés.
-   **Évaluation de l'intégrité** : Génération de rapports indiquant le pourcentage de gènes complets, dupliqués, fragmentés et manquants.

BUSCO génère un fichier FASTA par gène. L'ensemble des fichiers FASTA des gènes classés complets (`single_copy_busco_sequences`) est compilé sous un seul fichier FASTA par espèce.

### Miniprot

Utilisé pour l'annotation des génomes en se basant sur des séquences de protéines. Fonctionnement détaillé :

-   **Préprocessing des séquences de protéines** : Miniprot utilise des séquences de protéines bien annotées pour le mappage.
-   **Mappage des protéines** : Les séquences sont alignées sur le génome cible pour identifier les régions correspondantes.
-   **Identification des exons et des introns** : Miniprot identifie les exons (codants) et les introns (non codants).
-   **Transfert des annotations** : Miniprot transfère les annotations de protéines au génome cible.
-   **Génération des rapports** : Les résultats sont produits au format GFF.

### Scipio

Scipio est utilisé pour identifier les gènes codants et transférer des annotations fonctionnelles. Fonctionnement détaillé :

-   **Alignement de protéines** : Scipio aligne les séquences de protéines sur le génome cible.
-   **Identification des exons et introns** : Identification des exons et introns dans les régions alignées.
-   **Annotation des gènes** : Scipio transfère les annotations fonctionnelles des protéines vers le génome cible.
-   **Résultats** : Scipio produit des fichiers au format YAML avec les annotations de gènes.

L'annotation des génomes est réalisée en trois phases, chaque phase utilisant les séquences protéiques d'une espèce différente. Chaque génome est donc annoté trois fois.

Les fichiers YAML issus de Scipio sont convertis en GFF, puis plusieurs étapes sont réalisées pour traiter et nettoyer les résultats des annotations de Miniprot et Scipio :

### Traitement des Fichiers GFF

1.  **Copie** des trois fichiers GFF.
2.  **Extraction des CDS** des fichiers GFF (`listAlternativeTranscript2.py`).
3.  **Nettoyage des fichiers GFF** pour éliminer les transcrits alternatifs (`ExcludeCDSFromGFF.py`).
4.  **Concaténation des fichiers GFF** des CDS nettoyés.
5.  **Conversion des fichiers GFF en FASTA** pour obtenir les séquences exoniques (`putIDinNameGFF.py`).
6.  **Concatenation** des exons en transcrits complets (`concateExonFromBedtoolsGetFasta2_forCat.py`).

### Sélection et Traduction des Gènes

1.  **Tri des gènes** pour éliminer les exons chevauchants (`SelectOverlappingExonGFFsort_forCat.py`).
2.  **Récupération** des transcrits uniques (`copySeqList.py`).
3.  **Traduction des transcrits** en séquences protéiques.
4.  **Liste des transcrits** redondants (`RenameTranscriptRedundants_forList.py`).
5.  **Renommage des transcrits** redondants (`RenameTranscriptRedundants_forFasta.py`).
6.  **Sélection des transcrits** finaux (`SelectSeq.py`).
7.  **Renommage des transcrits** avec des informations spécifiques à l'espèce.

Des scripts bash ont été créés pour automatiser ces étapes pour les sorties de **Miniprot** et **Scipio** : - `Script_concatenateGFF_Miniprot.bash` - `Script_concatenateGFF_Scipio.bash`

Après l'annotation, nous avons traité les résultats de manière indépendante pour chaque jeu de gènes produits par **BUSCO**, **Miniprot** et **Scipio**.

------------------------------------------------------------------------

## 2) Recherche de Gènes Orthologues
