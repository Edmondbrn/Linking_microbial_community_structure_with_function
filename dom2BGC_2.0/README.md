
# Dom2BGC update

This part is dedicated to the update of the dom2BGC pipeline.



- First you have to run the pre_process.sh script to prepare your data from DADA2 to the pipeline. The script takes two arguments the name of the file containing the forward ASV fasta sequences (from DADA2) and the frame to translate the ASV into amino acids. (-frame 6 will translate into the six frames).
It will translate the ASVs, aligned them to the AMP-binding HMM profile and will compoute the phylogeny tree using FastTreeMP

```bash
bash pre_process.sh -sequence dna_fasta -f 6
```

- Then you can launch the boot_pipeline.sh. It will ask the path to the folder withh all the data such as the aligned amino acid sequences, the phylogeny tree and the feature table from DADA2. Then the pipeline will be launched.


```bash
bash boot_pipeline.sh
```

- The modifications of the pipeline involve the addition of 8 new functions:

- `divide_amplicon_list`
- `sub_cluster` (perform the clusterisation in sub groups of sequences)
- `regroup_cluster` (select all the unique clusters)
- `good_data` (format some pandas table)
- `merge_cluster` (create a structure to gather all the information)
- `merge_sequence` (create the output table with the clusters and the corresponding sequences)
- `merge_amplicon` (create the file with all the clusters and their sequences)
- `merge_feature` (create the count table for all the sub groups)

- The update of the database consists to the parsing of the genbank files from Mibig and Antismash databases. You can find the Mibig database v3.1 on their website and for ANtismash you need to perform a query with AMP-binding protein profil.
You can launch the java script by commenting the two last lines. You will get the organism list file. You have to put this file in NCBI website to get the ouput with taxaID. You can place this file in the ressources/data folders. Then you can comment the three first lines and uncommenting the two last ones.

```java
javac Parsing.java && java Parsing
```