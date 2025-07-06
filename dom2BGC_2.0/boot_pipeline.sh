#!/usr/bin/bash -

# python3.10 amplicon_pipeline.py --featureTable "data/table.from_biom.tsv" --phylogenyTree "data/protein.faa.tree" --amplicons "data/protein.faa" --MiBIGAmplicons databases/parsed_AMP_MIBIG.faa --antismashDBAmplicons databases/parsed_AMP_antismashdb2.faa --metagenomeAmplicons databases/parsed_AMP_metagenome.faa --parsedAsdbGenbank databases/asdb_parsed_product_list.txt --MiBIGOutputDir OUT_testt/mibig/ --antismashOutputDir OUT_testt/asdb/ --metagenomeOutputDir OUT_testt/metagenome/ --outputDir OUT_testt/test/ --parsedTaxonomy databases/antismashdb2_taxa_dict.txt --sampleNames I23-1089,I23-1093,I23-1095 --replicateList I23-1089-A01,I23-1089-A02,I23-1089-A03,I23-1089-A04,I23-1089-A05,I23-1089-B01,I23-1089-B02,I23-1089-B03,I23-1089-B04,I23-1089-B05,I23-1089-C01,I23-1089-C02,I23-1089-C03,I23-1089-C04,I23-1089-D01,I23-1089-D02,I23-1089-D03,I23-1089-D04,I23-1089-E01,I23-1089-E02,I23-1089-E03,I23-1089-E04,I23-1089-F01,I23-1089-F02,I23-1089-F03,I23-1089-F04,I23-1089-G01,I23-1089-G02,I23-1089-G03,I23-1089-G04,I23-1089-H01,I23-1089-H02,I23-1089-H03,I23-1089-H04,I23-1093-A01,I23-1093-A02,I23-1093-A03,I23-1093-A04,I23-1093-B01,I23-1093-B02,I23-1093-B03,I23-1093-B04,I23-1093-C01,I23-1093-C02,I23-1093-C03,I23-1093-C04,I23-1093-D01,I23-1093-D02,I23-1093-D03,I23-1093-D04,I23-1093-E01,I23-1093-E02,I23-1093-E03,I23-1093-E04,I23-1093-F01,I23-1093-F02,I23-1093-F03,I23-1093-F04,I23-1093-G01,I23-1093-G02,I23-1093-G03,I23-1093-H01,I23-1093-H02,I23-1093-H03,I23-1095-A01,I23-1095-A03,I23-1095-A04,I23-1095-A05,I23-1095-A06,I23-1095-A07,I23-1095-A08,I23-1095-B01,I23-1095-B03,I23-1095-B04,I23-1095-B05,I23-1095-B06,I23-1095-B07,I23-1095-B08,I23-1095-C01,I23-1095-C03,I23-1095-C04,I23-1095-C05,I23-1095-C06,I23-1095-C07,I23-1095-C08,I23-1095-D01,I23-1095-D03,I23-1095-D04,I23-1095-D05,I23-1095-D06,I23-1095-D07,I23-1095-D08,I23-1095-E01,I23-1095-E03,I23-1095-E04,I23-1095-E05,I23-1095-E06,I23-1095-E07,I23-1095-E08,I23-1095-F03,I23-1095-F04,I23-1095-F05,I23-1095-F06,I23-1095-F07,I23-1095-F08,I23-1095-G03,I23-1095-G05,I23-1095-G06,I23-1095-G07,I23-1095-H02,I23-1095-H03,I23-1095-H05,I23-1095-H06,I23-1095-H07 --separator '-' --rarefactionThreshold 500 --minReplicates 3 --verbose


echo "Please run this bash file in the dom2BGC-master folder."
echo "What do you want to do ?"
echo "Run the program ? (type run)"
echo "Install dependencies ? (type install)"
read -p "Your choice (run/install)" choice

run_process(){
    echo "$2"
    echo "Current directory : $(data)"
    read -p ".tsv name file (with extension):" tsv
    read -p "Tree file (with extension):" Tree
    read -p "Sequence (protein file) (with extension):" protein
    concatenated_tsv="${2}/${tsv}"
    concatenated_seq="${2}/${protein}"
    concatenated_tree="${2}/${Tree}"
    python3.10 amplicon_pipeline.py \
    --featureTable "$concatenated_tsv" --phylogenyTree "$concatenated_tree" --amplicons "$concatenated_seq" \
    --MiBIGAmplicons database2/parsed_AMP_MIBIG.faa --antismashDBAmplicons database2/parsed_AMP_antismashdb2.faa \
    --metagenomeAmplicons database2/parsed_AMP_metagenome.faa --parsedAsdbGenbank database2/asdb_parsed_product_list.txt \
    --MiBIGOutputDir "$1"/mibig/ --antismashOutputDir "$1"/asdb/ --metagenomeOutputDir "$1"/metagenome/ --outputDir "$1"/test/ \
    --parsedTaxonomy database2/antismashdb2_taxa_dict2.txt --sampleNames I23-1089,I23-1093,I23-1095 \
    --replicateList I23-1089-A01,I23-1089-A02,I23-1089-A03,I23-1089-A04,I23-1089-A05,I23-1089-B01,I23-1089-B02,I23-1089-B03,I23-1089-B04,I23-1089-B05,I23-1089-C01,I23-1089-C02,I23-1089-C03,I23-1089-C04,I23-1089-D01,I23-1089-D02,I23-1089-D03,I23-1089-D04,I23-1089-E01,I23-1089-E02,I23-1089-E03,I23-1089-E04,I23-1089-F01,I23-1089-F02,I23-1089-F03,I23-1089-F04,I23-1089-G01,I23-1089-G02,I23-1089-G03,I23-1089-G04,I23-1089-H01,I23-1089-H02,I23-1089-H03,I23-1089-H04,I23-1093-A01,I23-1093-A02,I23-1093-A03,I23-1093-A04,I23-1093-A05,I23-1093-A06,I23-1093-A07,I23-1093-A08,I23-1093-B01,I23-1093-B02,I23-1093-B03,I23-1093-B04,I23-1093-B05,I23-1093-B06,I23-1093-B07,I23-1093-B08,I23-1093-C01,I23-1093-C02,I23-1093-C03,I23-1093-C04,I23-1093-C05,I23-1093-C06,I23-1093-C07,I23-1093-C08,I23-1093-D01,I23-1093-D02,I23-1093-D03,I23-1093-D04,I23-1093-D05,I23-1093-D06,I23-1093-D07,I23-1093-D08,I23-1093-E01,I23-1093-E02,I23-1093-E03,I23-1093-E04,I23-1093-E05,I23-1093-E06,I23-1093-E07,I23-1093-E08,I23-1093-F01,I23-1093-F02,I23-1093-F03,I23-1093-F04,I23-1093-F05,I23-1093-F06,I23-1093-F07,I23-1093-F08,I23-1093-G01,I23-1093-G02,I23-1093-G03,I23-1093-G05,I23-1093-G06,I23-1093-G07,I23-1093-H01,I23-1093-H02,I23-1093-H03,I23-1093-H05,I23-1093-H06,I23-1093-H07,I23-1095-A01,I23-1095-A02,I23-1095-A03,I23-1095-A04,I23-1095-A05,I23-1095-A06,I23-1095-A07,I23-1095-A08,I23-1095-B01,I23-1095-B02,I23-1095-B03,I23-1095-B04,I23-1095-B05,I23-1095-B06,I23-1095-B07,I23-1095-B08,I23-1095-C01,I23-1095-C02,I23-1095-C03,I23-1095-C04,I23-1095-C05,I23-1095-C06,I23-1095-C07,I23-1095-C08,I23-1095-D01,I23-1095-D02,I23-1095-D03,I23-1095-D04,I23-1095-D05,I23-1095-D06,I23-1095-D07,I23-1095-D08,I23-1095-E01,I23-1095-E02,I23-1095-E03,I23-1095-E04,I23-1095-E05,I23-1095-E06,I23-1095-E07,I23-1095-E08,I23-1095-F01,I23-1095-F02,I23-1095-F03,I23-1095-F04,I23-1095-F05,I23-1095-F06,I23-1095-F07,I23-1095-F08,I23-1095-G01,I23-1095-G02,I23-1095-G03,I23-1095-G05,I23-1095-G06,I23-1095-G07,I23-1095-H01,I23-1095-H02,I23-1095-H03,I23-1095-H05,I23-1095-H06,I23-1095-H07 \
    --separator '-' --rarefactionThreshold 500 --minReplicates 3 --verbose
    read -p "Operation over. Please press enter to close the program" end
    exit
}


if [ "$choice" == "run" ]; then
    read -p "Where are your data folder (.tsv file, protein sequence etc...) ?" data
    while true; do
        read -p "Output folder name : " name
        if [ -d "$name" ]; then
            read -p "The directory already exists, do you want to overwrite it or start analysis from it? (overwrite/start)" answer
            if [ "$answer" == "overwrite" ]; then
                rm -rf "$name"
                run_process "$name" "$data"
            else
                run_process "$name" "$data"
            fi
        else
            run_process "$name" "$data"
        fi
    done


elif [ "$choice" == "install" ]; then
    echo "Module installation"
    python3.10 -m pip install scikit-learn
    python3.10 -m pip install scikit-bio
    python3.10 -m pip install numpy
    python3.10 -m pip install pandas
    python3.10 -m pip install matplotlib
    python3.10 -m pip install numba
    python3.10 -m pip install seaborn
    python3.10 -m pip install ete3
    python3.10 -m pip install networkx
    echo
    echo "Installation successfull"
    echo "If any packages are missing please run the command :"
    echo "python3.10 -m pip install <package_name> (without the <>)"

else
    echo "Please enter a correct answer (run/install/tree)"
fi