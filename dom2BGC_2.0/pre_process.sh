#!/usr/bin/bash -

help(){
    echo
    echo "This program aims to process data from DADA2 in order to perform dom2BCG analysis"
    echo "It will translate the ASV sequences into protein in one specific frame"
    echo "If the dna sequence file is not in the same directory as this bash script, you have to specify the path to the dna sequence file in a second time"
    echo "Then it will execute a hmmsearch with an AMP binding protein domain as reference"
    echo "Finally it will create a phylogeny tree thanks to FastTree (mulithreading)"
    echo
    echo "Parameters: -s DNA sequences"
    echo "            -f frame for the translation (1,2,3,-1,-2,-3,F,R,6) see transeq -h section -frame for the details"
    echo 
    echo "Example: bash pre_process.sh -s dna_sequence.fasta -f 1"
    exit
}

if [ "$1" == "-h" ] || [ "$#" -ne 4 ]; then
    help

elif [ "$1" != "-s" ] || [ "$3" != "-f" ]; then
    echo
    echo "Please enter the parameters as shown in the example"
    help
else 
    read -p "Enter the absolute path to your dna sequence (without the file at the end). Exampe: /home/edmond/data-dna : " path
    read -p "Enter the name of the output file (example: protein.faa)" out
    concatenated_dna="${path}/${2}"
    concatenated_out_prot="${path}/${out}"
    transeq -sequence "$concatenated_dna" -frame "$4" -outseq "$concatenated_out_prot"
    echo "Translation over"
    hmm="hmm_file"
    concatenated_hmm="${path}/${hmm}"
    hmmsearch -o "$concatenated_hmm" ./hmm_profiles/PF00501.27.hmm "$concatenated_out_prot"
    echo "hmmsearch over"
    out="processed_sequence.faa"
    concatenated_out_processed="${path}/${out}"
    python3.10 ./hmm_profiles/parse_hmm.py -i "$concatenated_hmm" -o "$concatenated_out_processed"
    echo "Parsing over. Beginning of the phylogeny tree creation this step can take a moment"
    export OMP_NUM_THREADS=16 # fix the number of threads used by FastTree
    out_tree="sequence.tree"
    concatenated_out_tree="${path}/${out_tree}"
    ./FastTreeMP "$concatenated_out_processed" > "$concatenated_out_tree"
    echo "Operation is over. You can now launch amplicon_pipeline.py with the processed_protein sequence and the pylogeny tree."
fi