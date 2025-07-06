#/usr/bin/python

'''
Author: Vittorio Tracanna
Use:    Script to calculate distance between in-silico and regular amplicons. 
        group in-silico amplicons based on BGC. 
        calculate average distance for all BGC.
        evaluate if BGC is in data or not.
'''
#imports
import pickle
from sys import getsizeof
import argparse
import os
import subprocess
import time
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
from ete3 import NCBITaxa
from io import StringIO
from skbio import TreeNode
from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio.stats import subsample_counts
from skbio.stats.ordination import pcoa
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pairwise_distance_calculation import *
import csv

def commandLineParser():
    '''
    Parses input and holds default values.
    Returns a dictionary with all the parameters for the wrapper
    '''
    parser = argparse.ArgumentParser(description='Annotates amplicons from a feature table with available databases like antiSMASHdb, MiBIG and paired metagenomes')
    local = os.getcwd()
    #General argnuments to specify where to put/find files
    parser.add_argument('-f', '--featureTable', type=str, required=True, help='path to the tab separated feature table with the amino acid amplicon sequences and counts')
    parser.add_argument('-amp', '--amplicons', type=str, required=True, help='path to the amplicons amino acid sequences for annotation')
    parser.add_argument('-asdb', '--antismashDBAmplicons', type=str, required=True, help='path to the antismashdb in silico amplicons for annotation')
    parser.add_argument('-mibig', '--MiBIGAmplicons', type=str, required=True, help='path to the MiBIG in silico amplicons for annotation')
    parser.add_argument('-meta', '--metagenomeAmplicons', type=str, required=False, help='path to the metagenome in silico amplicons')
    parser.add_argument('-l', '--lipopeptideList', type=str, required=False, help='txt file with the list of BGCs to consider lipopeptides [if available]')
    parser.add_argument('-tree', '--phylogenyTree', type=str, required=False, default=None, help='newick file to calculate phylogeny based metrics [UniFrac]')
    parser.add_argument('-mibigout', '--MiBIGOutputDir', type=str, required=False, help='path to the directory with the MiBIG output files')
    parser.add_argument('-asdbout', '--antismashOutputDir', type=str, required=False, help='path to the directory for the antismashdb output files')
    parser.add_argument('-asdbtaxa', '--parsedTaxonomy', type=str, required=False, help='path to the parsed antismashdb taxonomy, if not provided it will create the needed file, long runtime')
    parser.add_argument('-metaout', '--metagenomeOutputDir', type=str, required=False, help='path to the directory for the shotgun metagenome output files')
    parser.add_argument('-asdbgbk', '--asdbGenbank', type=str, required=False, help='path to the directory with the antismashdb genbank files')
    parser.add_argument('-parsedProd', '--parsedAsdbGenbank', type=str, required=False, help='path to the preparsed list of antismashdb products')
    parser.add_argument('-mibiggbk', '--MiBIGGenbank', type=str, required=False, help='path to the directory with the MiBIG genbank files')
    parser.add_argument('-o', '--outputDir', type=str, required=False, default=local, help='Main output directory, generates different subfolders to organize output')
    parser.add_argument('-s', '--sampleNames', type=str, default='', required=False, help='Comma separated list of sample names, for sample replicates do specify the separator. ex: --sampleNames S01,S02,S03')
    parser.add_argument('-r', '--replicateList', type=str, required=False,help='Comma separated list of replicates, for sample replicates do specify the separator. ex: --replicateList S01_1,S01_2,S02_1,S02_2,S03_1,S03_2')
    parser.add_argument('-sep', '--separator', type=str, required=False, default='_', help='Separator between sample name and replicate number. Default=_. ex: --separator _')
    parser.add_argument('-t', '--networkThreshold', type=float, required=False, default=0.75, help='Correlation threshold to apply to the cooccurence matrix to generate the network. Default=0.75')
    parser.add_argument('-rt', '--rarefactionThreshold', type=int, required=False, default=30000, help='Rarefaction level for all the replicates, Default=30000')
    parser.add_argument('-min', '--minReplicates', type=int, required=False, default=3, help='Minimal number of replicates an amplicon should appear in the feature table for its inclusion in the cooccurrence network. Default=3.')
    parser.add_argument('-v', '--verbose', required=False, default=False, action='store_true', help='Set verbosity while running')
    return vars(parser.parse_args())

def load_amplicons(feature_table_index, argOptions):
    '''
    Loads amplicons for both normal amplicon data and the in silico amplicons from reference dbs
    Sequences should be translated to aa sequence and prealigned using the same hmm profile as described in [add reference]
    '''
    amplicons_header_len = len(feature_table_index[0])
    amplicons = open(argOptions['amplicons']).read()
    amplicons_list = [[x.split('\n')[0][:-4], x.split('\n')[1]] for x in amplicons.split('>')[1:] if x.split('\n')[0][:-4] in feature_table_index]
    mibig_amplicons = open(argOptions['MiBIGAmplicons']).read()
    mibig_amplicons_list = [[x.split('\n')[0], x.split('\n')[1]] for x in mibig_amplicons.split('>')[1:]]
    antismashdb_amplicons = open(argOptions['antismashDBAmplicons']).read()
    antismashdb_amplicons_list = [[x.split('\n')[0], x.split('\n')[1]] for x in antismashdb_amplicons.split('>')[1:]]
    if argOptions['metagenomeAmplicons']:
        metagenome_amplicons = open(argOptions['metagenomeAmplicons']).read()
        metagenome_amplicons_list = [[x.split('\n')[0], x.split('\n')[1]] for x in metagenome_amplicons.split('>')[1:]]
    else:
        metagenome_amplicons_list = []
    return amplicons_list, mibig_amplicons_list, antismashdb_amplicons_list, metagenome_amplicons_list

def group_similar_amplicons(feature_table, amplicon_list, argOptions, version :str = "", is_sub : bool = False):
    '''
    Groups identical amplicons with agglomerative clustering. Reads are summed and assigned to the new amplicon.
    Renames amplicons in the format amplicon_cluster_{} where {} is the amplicon number. Helps avoiding confusion in the next steps
    '''

    if not (os.path.isfile('{}grouped_amps.tsv'.format(argOptions['outputDir'])) and os.path.isfile('{}amplicons_list.p'.format(argOptions['outputDir'])) and os.path.isfile('{}grouped-feature-table.tsv'.format(argOptions['outputDir']))):
        from sklearn.cluster import AgglomerativeClustering

        new_soil_amplicon_list = []
        soil_amplicons_seq = [x[1] for x in amplicon_list]
        soil_amplicons_len = [len([x for x in y if x!='-']) for y in soil_amplicons_seq]
        soil_amplicons_header = [x[0] for x in amplicon_list]
        # calculates pairwise sequence identity
        pw_identity = pd.DataFrame(pairwise_distances(soil_amplicons_seq, soil_amplicons_seq, soil_amplicons_len, soil_amplicons_len), index=soil_amplicons_header, columns=soil_amplicons_header)

        clt = AgglomerativeClustering(n_clusters=None, distance_threshold=0.01, linkage='single')
        cal = np.subtract(pw_identity, np.ones((len(pw_identity), len(pw_identity)))) * -1
        model = clt.fit(cal)
        new_feature_table = pd.DataFrame(columns=feature_table.columns)
        out_report = ''
        for similar_domain_group in set(model.labels_):
            positions = [x for x in range(len(model.labels_)) if model.labels_[x] == similar_domain_group]
            if len(positions):
                amplicons_to_group = pw_identity.index[positions]
                new_feature_table.loc['amplicon_cluster_{}'.format(similar_domain_group)] = feature_table.loc[amplicons_to_group].sum()
                out_report += 'amplicon_cluster_{}\t{}\n'.format(similar_domain_group, '\t'.join(amplicons_to_group))
                new_soil_amplicon_list.extend([['amplicon_cluster_{}'.format(similar_domain_group), soil_amplicons_seq[positions[0]]]])
            else:
                amplicon = pw_identity.index[positions]
                out_report += 'amplicon_cluster_{}\t{}\n'.format(similar_domain_group, amplicon)
                new_feature_table.loc['amplicon_cluster_{}'.format(similar_domain_group)] = feature_table.loc[amplicon]
                new_soil_amplicon_list.extend(['amplicon_cluster_{}'.format(similar_domain_group), soil_amplicons_seq[positions]])
        
        if is_sub: # change the ouput dir if we have to process sub amplicons
            outputdir = f"{argOptions['outputDir']}amplicon_processing/output/"
        else:
            outputdir = f"{argOptions['outputDir']}"
        new_feature_table.to_csv(f'{outputdir}grouped-feature-table{version}.tsv', sep='\t')
        header = '# Constructed from biom file\n'
        outcsv = open(f'{outputdir}grouped-feature-table{version}.tsv').read()
        out_csv_txt = header+outcsv
        outCsvFile = open(f'{outputdir}grouped-feature-table{version}.tsv', 'w')
        outCsvFile.write(out_csv_txt)
        outCsvFile.close()
        outputFile = open(f'{outputdir}grouped_amps{version}.tsv','w')
        outputFile.write(out_report)
        outputFile.close()
        pickle.dump(new_soil_amplicon_list, open(f'{outputdir}amplicons_list{version}.p', 'wb'))

    else:
        new_soil_amplicon_list = pickle.load(open(f'{outputdir}amplicons_list{version}.p', 'rb'))
        new_feature_table = pd.read_csv(f'{outputdir}grouped-feature-table{version}.tsv', sep='\t', index_col=0, header=1)

    return new_feature_table, new_soil_amplicon_list

def extract_amplicon_name(feature_table_index : pd.DataFrame, amplicon_liste_path : str):
    amplicons = open(amplicon_liste_path).read()
    amplicons_list = [[x.split('\n')[0][:-4], x.split('\n')[1]] for x in amplicons.split('>')[1:] if x.split('\n')[0][:-4] in feature_table_index]
    return amplicons_list


def divide_amplicon_list(amplicon_path : str, n : int) -> None:
    """
    Function to divide the amplicon list into n subfiles
    """
    amp = open(amplicon_path, "r")
    lines = amp.readlines()
    length = len(lines) // n # get the number of lines per sub file
    cpt = 1
    index = 0
    if length % 2 != 0: # we have to have an even number of lines (header + sequence)
        length += 1
    for line in lines:
        if index % length == 0: # if we reach the number of lines per sub file
            if os.path.exists(f"{argOptions['outputDir']}amplicon_processing/sub_amplicon/amplicon_{cpt}.faa"):
                fh.close()
                cpt += 1
            fh = open(f"{argOptions['outputDir']}amplicon_processing/sub_amplicon/amplicon_{cpt}.faa", "w")
        fh.write(line)
        index += 1
    fh.close()
    amp.close()
    print("Sub files created")

def sub_cluster(list_ampl : list[list[str]], df : pd.DataFrame) -> None:
    """
    Function to cluster all the sub amplicons
    """
    tot = len(list_ampl)
    for i, element in enumerate(list_ampl):
        if not os.path.isfile(f"{argOptions['outputDir']}amplicon_processing/output/amplicons_list_{i+1}.p"):
            _, _ = group_similar_amplicons(df, element, argOptions ,f"_{i+1}", True)
        print(f"Sub amplicon \r{i+1} / {tot} done", sep = "", end = "")
    print(f"Sub amplicon {tot} / {tot} done", sep = "", end = "")

def regroup_cluster(path_sub :str) -> list[str]:
    """
    Function to regroup all the clusters in one structure
    """
    list_cluster = [] # get the list of all clusters
    for element in os.listdir(path_sub):
        if element.endswith(".p"):
            clus = pd.read_pickle(f"{path_sub}/{element}")
            for seq in clus:
                list_cluster.append(seq[1])
    # remove duplicates with a set
    list_cluster = list(set(list_cluster))
    return list_cluster


def good_data(path: str) -> list[pd.DataFrame]:
    """
    Function that add NaN in blank columns if the dataframe doesn't have the same number of columns
    """
    df_list = []
    list_file = os.listdir(f"{path}")
    nb = len([x for x in list_file if x.startswith("grouped_amps_")])
    for i in range(1, nb +1):  # Browse all the files
        file_path = f"{path}grouped_amps_{i}.tsv"
        max_cols = 0
        with open(file_path, newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                max_cols = max(max_cols, len(row))  # find the maximum number of columns
        col_names = [f"col_{j}" for j in range(max_cols)]  # Create the column names (add NaN in blank columns)
        df_list.append(pd.read_csv(file_path, sep="\t", header=None, names=col_names))
    return df_list  # Liste de DataFrame corrects

def merge_cluster(path_sub:str, n : int) -> dict[str : dict[str : list]]:
    """
    Function that merge all the sub clusters in one dictionnary
    Output : dict with the cluster as key and a dictionnary as value. 
            This dictionnary contains as keys the ID of the amplicon in its file and as value the list of associated sequences (ASVs / OTUs)
    """
    list_cluster = regroup_cluster(path_sub)
    dico_cluster = {}
    cpt = 1
    data_list = [pd.read_pickle(f"{path_sub}amplicons_list_{i}.p") for i in range(1, n+1)] # we have a list for each files, and the list is divided in sub list of 2 elements [ID_filename, sequences]
    seq_ass_list = good_data(path_sub)
    for cluster in list_cluster: # browse the unique clusters
        dico_sub = {}
        print(f"Cluster \r{cpt} / {len(list_cluster)} done", sep ="\n", end = "", flush = True)
        for i in range(1, n+1): # browse the files
            data = data_list[i-1]
            seq_ass = seq_ass_list[i-1]
            for sub_clustr in data: # browse the cluster of the file
                if sub_clustr[1] == cluster: # if the studied sequence correponds to the current unique cluster
                    number = int(sub_clustr[0].split("_")[-1])
                    try:
                        filtered_result = seq_ass.loc[seq_ass["col_0"] == sub_clustr[0]].values.tolist() # take the line in the table
                        filtered_result = [value for value in filtered_result[0] if not pd.isna(value)] # remove NaN if introduced by good_data()
                    except:
                        print("Here the number" ,number)
                        print(i)
                        print(seq_ass.shape)
                        quit()
                    else:
                        dico_sub[f"{i}_{sub_clustr[0]}"] = filtered_result[1:] # remove the fisrt one which is the cluster name
                        break # break because the clusters are unique in each files
              
        dico_cluster[cluster] = dico_sub
        cpt += 1
    print("Saving the dictionnary...")
    with open(f"{argOptions['outputDir']}amplicon_processing/output/dico_cluster.p", 'wb') as file:
        pickle.dump(dico_cluster, file)
    return dico_cluster

def merge_sequence(dico_cluster : dict[str : dict[str : list]]) -> None:
    """
    Function that extract the name of unique clusters and that create a dictionnary with the cluster and the associated sequences
    """
    dico_final = {}
    tot = len(dico_cluster.keys())
    for index, cluster in enumerate(dico_cluster.keys()): # browse all the clusters
        list_seq = []
        for seq in dico_cluster[cluster].values(): # get the list of sequences that matchs with the cluster
            list_seq.extend(seq) # regroup theses sequences
        dico_final[f"amplicon_cluster_{index}"] = list_seq # creation of a new cluster with the sequences
        print(f"Cluster \r{index} /{tot} done", sep = "", end = "")
    print(f"Cluster {tot} /{tot} done", sep = "", end = "")
    print("Creating the final sequence dataframe ...")
    df = pd.DataFrame({key : pd.Series(value) for key, value in dico_final.items()})
    df = df.transpose()
    df.to_csv(f"{argOptions['outputDir']}grouped_amps.tsv", sep = "\t")

def merge_amplicon(dico_cluster : dict[str : dict[str : list]]) -> None:
    """
    Function that merge all the amplicons in one list to serialise it (to get it more easily in the next steps)
    """
    list_amplicon = []
    for index, cluster in enumerate(dico_cluster.keys()):
        liste = [f"amplicon_cluster_{index}", cluster]
        list_amplicon.append(liste)
    pickle.dump(list_amplicon, open(f"{argOptions['outputDir']}amplicon_list.p", "wb"))


def merge_feature(dico_cluster : dict[str : dict[str : list]], n : int) -> None:
    """
    Function that merges all the feature table in one dataframe
    """
    liste_df_feature = [pd.read_csv(f"{argOptions['outputDir']}amplicon_processing/output/grouped-feature-table_{i}.tsv", sep = "\t", low_memory=False, index_col=0, header = 1) for i in range(1, n+1)]
    df = pd.DataFrame(0, index = [f"amplicon_cluster_{i}" for i in range(len(dico_cluster.keys()))], 
                      columns = liste_df_feature[0].columns)
    max_ = len(dico_cluster.keys())
    for index, cluster in enumerate(dico_cluster.keys()): # browse all the clusters
        for id in dico_cluster[cluster].keys(): # browse all the ID of the cluster in sub file
            id_file = id.split("_") # get the number of the file
            df_sub = liste_df_feature[int(id_file[0])-1] # get the dataframe of the file
            line_df = "_".join(id_file[1:]) # get the name of the cluster in the sub files
            line = np.array(df_sub.loc[line_df]) # get the line of the cluster
            df.loc[f"amplicon_cluster_{index}"] += line # add the line to the final dataframe
        print(f"Cluster \r{index} /{max_} done", sep = "", end = "")
    print(f"Cluster {max} /{max} done", sep = "", end = "")
    print("Sums verification")
    somme_sous_dfs = sum([df_sub.sum().sum() for df_sub in liste_df_feature])
    somme_df_final = df.sum().sum()
    print(f"Sum of the final dataframe : {somme_df_final} and sum of alll the sub datframe : {somme_sous_dfs}")
    df.to_csv(f"{argOptions['outputDir']}grouped-feature-table.tsv", sep = "\t")


def group_metagenome_amplicons(in_silico_amplicons):
    '''groups BGCs that start with the same code'''
    amplicon_headers = [x[0] for x in in_silico_amplicons]
    bgc_set = list(set([x.split('#')[0] for x in amplicon_headers]))
    bgc_dict = {}
    for bgc in bgc_set:
        bgc_dict[bgc] = [x for x in amplicon_headers if x.startswith(bgc+"_")]
    return bgc_dict

def make_rare_table(filtered_table, treshold):
    '''
    Creates a rare table to normalize all samples, set the threshold throught the --rarefactionThreshold flag
    :param filtered_table: 
    :param treshold: 
    :return: 
    '''
    rare_table = pd.DataFrame(0, index=filtered_table.index, columns=filtered_table.columns)
    
    for sample in filtered_table.columns:
        try:
            rare_table.loc[:, sample] = rarefy(filtered_table.loc[:, sample], treshold)
        except ValueError:
            continue
    rare_table = rare_table.fillna(0)
    return rare_table

def rarefy(df, size):
    df = df.astype(int)
    names = df[df!=0].index
    prob = df[df!=0]
    rarefied = pd.Series(subsample_counts(prob, size, replace=False), index=names)
    return rarefied

def group_insilico_amplicons(in_silico_amplicon_list):
    import Grouping_BCG
    result = Grouping_BCG.group_insilico_amplicons(in_silico_amplicon_list)
    return result
    # amplicon_headers = [x[0] for x in in_silico_amplicon_list]
    # bgc_set = list(set(['_'.join(x.split('_')[:-1]) for x in amplicon_headers]))
    # bgc_dict = {}
    # i = 1
    # for bgc in bgc_set:
    #     bgc_dict[bgc] = [x for x in amplicon_headers if x.startswith(bgc)]
    #     print("Grouping BGCs {}/{}".format(i, len(bgc_set)))
    #     i += 1
    # return bgc_dict

def calculate_similarity(bgc_dict, soil_amplicons, db_amplicons, outputFolder, dbtype):
    '''
    Calculates similarity [pairwise identity] between db in silico amplicons and soil amplicons

    :param bgc_dict: 
    :param soil_amplicons: 
    :param mibig_amplicons: 
    :param filtered_table: 
    :return: 
    '''
    db_amplicons_seq = [x[1] for x in db_amplicons];
    soil_amplicons_seq = [x[1] for x in soil_amplicons];
    soil_amplicons_len = [len([x for x in y if x != '-']) for y in soil_amplicons_seq]
    db_amplicons_header = [x[0] for x in db_amplicons];
    soil_amplicons_header = [x[0] for x in soil_amplicons];
    db_amplicons_len = [len([x for x in y if (x != '-' and x != 'X')]) for y in db_amplicons_seq]

    # calculates pairwise sequence identity
    pw = pd.DataFrame(pairwise_distances(db_amplicons_seq, soil_amplicons_seq, db_amplicons_len, soil_amplicons_len),
                      index=db_amplicons_header, columns=soil_amplicons_header)
    grouped_bgc = {}
    tot = len(set(['_'.join(x.split('_')[:2]) for x in bgc_dict.keys()]))
    cpt = 0
    # creates a dict with all the in silico amplicons that belong to a BGC
    print("Creation of the dictionnary for amplicons that matchs with a BGC")
    for bgc_group in list(set(['_'.join(x.split('_')[:2]) for x in bgc_dict.keys()])):
        grouped_bgc[bgc_group] = []
        cpt+=1
        for bgc in pw.index:
            if bgc.startswith(bgc_group + '_'):
                grouped_bgc[bgc_group].append(bgc)
        print(f"BGC number \r{cpt} out of {tot} done", sep = "", end = "", flush = True)
    # large function that produces plots
    bgc_extended_dict = {}

    # creates dict to store all the matches for a db domain
    count = 0
    print("Creation of the dictionnary to store all the matches for a db domain")
    tot = len(grouped_bgc.keys())
    for bgc in grouped_bgc.keys():
        count += 1
        if len(grouped_bgc[bgc]) > 2:
            for bgc_dom in grouped_bgc[bgc]:
                matches = pw.loc[bgc_dom]
                if max(pw.loc[bgc_dom]) > 0.9:
                    matches = matches.loc[matches > 0.9].index
                    if bgc in bgc_extended_dict:
                        bgc_extended_dict[bgc][bgc_dom] = list(matches)
                    else:
                        bgc_extended_dict[bgc] = {bgc_dom: list(matches)}
        print(f"BGC number \r{count} out of {tot} done", sep = "", end = "", flush = True)
    # groups all the matches for all domains of a BGC
    BGC_to_soil = {}
    print("Grouping all the matches together")
    tot = len(bgc_extended_dict.keys())
    cpt = 0
    for k1, v in bgc_extended_dict.items():
        cpt += 1
        amps = []
        for k2 in v.keys():
            amps.extend(v[k2])
        print(f"BGC number \r{cpt} out of {tot} done", sep = "", end = "", flush = True)
        BGC_to_soil[k1] = amps

    # creates plots
    soil_annotation = {}

    # useful, restore later
    print("Restoration creation")
    tot = len(BGC_to_soil.keys())
    cpt = 0
    for k in BGC_to_soil.keys():
        cpt += 1
        for dom in BGC_to_soil[k]:
            if dom in soil_annotation.keys():
                soil_annotation[dom].extend([k])
            else:
                soil_annotation[dom] = [k]
        print(f"BGC number \r{cpt} out of {tot} done", sep = "", end = "", flush = True)
    # clusters = create_clusters(log_table, BGC_to_soil, bgc_extended_dict, soil_annotation)

    pickle.dump(soil_annotation, open('{}{}_soil_annotation.p'.format(outputFolder, dbtype), 'wb'))
    pickle.dump(BGC_to_soil, open('{}{}_BGC_to_soil.p'.format(outputFolder, dbtype), 'wb'))
    pickle.dump(bgc_extended_dict, open('{}{}_bgc_extended_dict.p'.format(outputFolder, dbtype), 'wb'))

    return soil_annotation, BGC_to_soil, bgc_extended_dict

def parse_antismashdb_bgc_family(argOptions):
    '''
    Parses the preparsed antismashdb compound file. 
    Currently disregards hrslactone.
    '''
    annotation = open(argOptions['parsedAsdbGenbank']).read().split('\n')[:-1]
    product_dict = {x.split('\t')[0] : x.split('\t')[1] for x in annotation}
    pickle.dump(product_dict, open('{}antismashdb_product_dict.p'.format(argOptions['outputDir']), 'wb'))
    return product_dict

def assign_taxonomy(BGC_antismashdb_dict, argOptions):
    '''
    Downloads taxonomy for each entry from ncbi and gets the taxid
    :param BGC_antismashdb_dict: 
    :param argOptions: 
    :return: 
    '''
    taxa = {}
    if not os.path.isfile(argOptions['parsedTaxonomy']):
        for BGC in BGC_antismashdb_dict.keys():
            acc = BGC.split('.')[0]
            cmd = "curl -s 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=xml' |   grep TSeq_taxid |   cut -d '>' -f 2 |   cut -d '<' -f 1 |   tr -d '\n';    echo;".format(acc)
            time.sleep(0.05)
            taxid = subprocess.check_output(cmd, shell=True)
            taxid = str(taxid)
            taxa[BGC] = taxid
        pickle.dump(taxa, open('{}antismashdb2_taxa_dict.p'.format(argOptions['antismashOutputDir']), 'wb'))
    else:
        taxa_file = open(argOptions['parsedTaxonomy']).read().split('\n')[:-1]
        for entry in taxa_file:
            BGC = entry.split('\t')[0]
            taxid = entry.split('\t')[1]
            taxa[BGC] = taxid
        pickle.dump(taxa, open('{}antismashdb2_taxa_dict.p'.format(argOptions['antismashOutputDir']), 'wb'))
    return taxa

def write_annotation_to_csv(feature_table, product_amplicons, amp_lca):
    '''
    Writes amplicon annotation to file
    :param feature_table: 
    :param product_amplicons: 
    :return: 
    '''

    NCBI = NCBITaxa()
    node_annot = pd.DataFrame(index=feature_table.index, columns=['product', 'product_class', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus','species'])

    if not os.path.isfile('{}amplicon_annotation.csv'.format(outputFolder)):
        for id1 in feature_table.index:
            BGC_class_id1 = [x for x in product_dict.keys() if id1 in product_dict[x]]
            BGC_id = [x for x in product_amplicons.keys() if id1 in product_amplicons[x]]
            if BGC_class_id1:
                id1_prod = '-'.join(BGC_class_id1)
                id1_prod = '-'.join(list(set(id1_prod.split('-'))))
                node_annot.loc[id1, 'product_class'] = id1_prod
            else:
                node_annot.loc[id1, 'product_class'] = 'unassigned'

            if BGC_id:
                id1_prod = '-'.join(BGC_id)
                id1_prod = '-'.join(list(set(id1_prod.split('-'))))
                node_annot.loc[id1, 'product'] = id1_prod
            else:
                node_annot.loc[id1, 'product'] = 'unassigned'

        canonical_groups = ['superkingdom','phylum','class','order','family','genus','species']
        for id in feature_table.index:
            for group in canonical_groups:
                try:
                    node_annot.loc[id, group] = NCBI.get_taxid_translator([amp_lca.loc[id, group]])[amp_lca.loc[id, group]]
                except (IndexError, KeyError) as e:
                    node_annot.loc[id, group] = 'unassigned'
                    continue

    node_annot.to_csv('{}amplicon_annotation.csv'.format(outputFolder), sep='\t')

def filter_correlation_table_and_annotate(corr_df, product_dict, outputFolder, amp_lca, BGC_dict, gbk_amps, argOptions, T=''):
    '''
    takes the correlation table and applies a set threshold to make a network with the strongest interactions
    :param T: Threshold [user set or 0.75 default]
    :return: 
    '''
    NCBI = NCBITaxa()

    if not T:
        arr = corr_df.values.flatten()
        arr = np.asarray([x for x in arr if x])
        arr = np.asarray([x for x in arr if not np.isnan(x)])
        print ('95th percentile {}. 99th percentile {}.'.format(np.percentile(arr, 95), np.percentile(arr, 99)))
        fig = sns.distplot(arr, bins=200)
        plt.savefig('{}distplot_corr.pdf'.format(outputFolder))

    else:
        if not os.path.isfile('{}{}_corr_network_annot_table.csv'.format(outputFolder, T)):
            filt_corr_df = corr_df>T
            out_nw_txt = 'amplicon_id1\tamplicon_id2\tcorrelation\tproduct_type_id1\tproduct_type_id2\ttaxonomy_id1\ttaxonomy_id2\tasdb_matches_id1\tasdb_matches_id2\n'
            for id1 in range(len(corr_df)):
                BGC_class_id1 = [x for x in product_dict.keys() if corr_df.index[id1] in product_dict[x]]
                for id2 in range(id1+1, len(corr_df)):
                    if filt_corr_df.iloc[id1, id2]:
                        id1_prod = 'unassigned'; id2_prod = 'unassigned'

                        if BGC_class_id1:
                            id1_prod = '-'.join(BGC_class_id1)
                            id1_prod = '-'.join(list(set(id1_prod.split('-'))))

                        BGC_class_id2 = [x for x in product_dict.keys() if corr_df.index[id2] in product_dict[x]]
                        if BGC_class_id2:
                            id2_prod = '-'.join(BGC_class_id2)
                            id2_prod = '-'.join(list(set(id2_prod.split('-'))))

                        try:
                            taxa1 = amp_lca.loc[corr_df.index[id1], 'genus']
                        except KeyError:
                            taxa1 = '0'
                        try:
                            taxa2 = amp_lca.loc[corr_df.index[id2], 'genus']
                        except KeyError:
                            taxa2 = '0'

                        try:
                            gbk1 = '-'.join(gbk_amps[corr_df.index[id1]])
                        except KeyError:
                            gbk1 = ''

                        try:
                            gbk2 = '-'.join(gbk_amps[corr_df.index[id2]])
                        except KeyError:
                            gbk2 = ''

                        out_nw_txt += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(corr_df.index[id1], corr_df.index[id2], corr_df.iloc[id1, id2], id1_prod, id2_prod, taxa1, taxa2, gbk1, gbk2)

            out_nw = open('{}{}_corr_network_annot_table.csv'.format(outputFolder, T), 'w')
            out_nw.write(out_nw_txt)
            out_nw.close()

    node_annot = pd.DataFrame(index=corr_df.index, columns=['product', 'product_class', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
    sampleNames = argOptions['sampleNames']
    if sampleNames:
        samples = argOptions['sampleNames'].split(',')
        replicates = argOptions['replicateList'].split(',')
        sample_dict = {}
        for sample in samples:
            sample_dict[sample] = []
            for replicate in replicates:
                if replicate.startswith(sample):
                    sample_dict[sample].append(replicate)

        node_annot = pd.DataFrame(index=corr_df.index, columns=['product', 'product_class', 'superkingdom','phylum','class','order','family','genus','species']+['{}_associated'.format(x) for x in samples])

    if not os.path.isfile('{}{}_node_annotation.csv'.format(outputFolder, T)):
        for id1 in corr_df.index:
            BGC_class_id1 = [x for x in product_dict.keys() if id1 in product_dict[x]]
            BGC_id = [x for x in BGC_dict.keys() if id1 in BGC_dict[x]]
            if BGC_class_id1:
                id1_prod = '-'.join(BGC_class_id1)
                id1_prod = '-'.join(list(set(id1_prod.split('-'))))
                node_annot.loc[id1, 'product_class'] = id1_prod
            else:
                node_annot.loc[id1, 'product_class'] = 'unassigned'

            if BGC_id:
                id1_prod = '-'.join(BGC_id)
                id1_prod = '-'.join(list(set(id1_prod.split('-'))))
                node_annot.loc[id1, 'product'] = id1_prod
            else:
                node_annot.loc[id1, 'product'] = 'unassigned'

        canonical_groups = ['superkingdom','phylum','class','order','family','genus','species']

        for id in corr_df.index:
            samples_associated = []
            if argOptions['sampleNames']:
                for sample in samples:
                    if feature_table.loc[id, sample_dict[sample]].sum() > len(sample_dict[sample]) - 1:
                        samples_associated.extend([sample])
            if samples_associated:
                for s in samples_associated:
                    node_annot.loc[id, '{}_associated'.format(s)] = 'Yes'
            for group in canonical_groups:
                try:
                    node_annot.loc[id, group] = NCBI.get_taxid_translator([amp_lca.loc[id, group]])[amp_lca.loc[id, group]]
                except (IndexError, KeyError) as e:
                    node_annot.loc[id, group] = 'unassigned'
                    continue
        node_annot.to_csv('{}{}_node_annotation.csv'.format(outputFolder, T), sep='\t')
    return node_annot

def parse_taxonomy(amp_taxonomy, outputFolder):
    '''
    Make taxonomy table for assigned amplicons
    :param amp_taxonomy: 
    :return: 
    '''
    if not os.path.isfile('{}antismashdb2_taxa_df.p'.format(outputFolder)):
        NCBI = NCBITaxa()
        canonical_groups = ['superkingdom','phylum','class','order','family','genus','species']
        taxa_df = pd.DataFrame(0, index=amp_taxonomy.keys(), columns=canonical_groups)
        for amp in amplicon_taxonomy:
            try:
                ptaxa = amplicon_taxonomy[amp][2:-3]
                ancestors = NCBI.get_lineage(ptaxa)
                #names = NCBI.get_taxid_translator(ancestors)
                ranks = NCBI.get_rank(ancestors)
                inv_ranks = {v: k for k, v in ranks.items()}
                _canonical = [x for x in canonical_groups if x in inv_ranks.keys()]
                taxa_df.loc[amp, _canonical] = [inv_ranks[x] for x in _canonical]
            except (TypeError, ValueError) as error:
                continue
        pickle.dump(taxa_df, open('{}antismashdb2_taxa_df.p'.format(outputFolder),'wb'))
    else:
        taxa_df = pickle.load(open('{}antismashdb2_taxa_df.p'.format(outputFolder),'rb'))
    return taxa_df

def assign_taxonomy_to_amplicons(matches, taxa_df, outputFolder):
    '''
    Function to get all the taxonomy for a given amplicon, returns a dict with all the information
    :param matches: dictionary with the format {BGC1:{dom1: [amp1, amp2], dom2: [amp3, amp4]}, BGC2:{etc..}}
    :param taxa_df: dataframe with taxonomy of the BGCs
    :param outputFolder: 
    :return: 
    '''
    print(taxa_df.head())
    taxa_df.to_csv('taxa_df.csv', sep='\t')
    amplicon_taxonomy_dict = {}
    for BGC in matches.keys():
        amplicon_matches = []
        for domain in matches[BGC].keys():
            amplicon_matches.extend(matches[BGC][domain])

        for amplicon in amplicon_matches:
            if amplicon in amplicon_taxonomy_dict.keys():
                amplicon_taxonomy_dict[amplicon].append(taxa_df.loc[BGC].values)
            else:
                print("Here the amplicon : {} ".format(amplicon))
                print("Here the BGC : {} ".format(BGC))
                amplicon_taxonomy_dict[amplicon] = [taxa_df.loc[BGC].values]
    pickle.dump(amplicon_taxonomy_dict, open('{}antismashdb_amplicon_taxonomy_dict.p'.format(outputFolder), 'wb'))
    return amplicon_taxonomy_dict

def assign_lca_taxonomy(amplicon_taxonomy_dict, outputFolder):
    canonical_groups = ['superkingdom','phylum','class','order','family','genus','species']
    amp_lca = pd.DataFrame(0, index=amplicon_taxonomy_dict.keys(), columns=canonical_groups)
    for amp in amplicon_taxonomy_dict.keys():
        if len(set([str(x) for x in amplicon_taxonomy_dict[amp]])) == 1:
            amp_lca.loc[amp, canonical_groups] = amplicon_taxonomy_dict[amp][0]
        else:
            for level in range(7):
                LCA = [x[level] for x in amplicon_taxonomy_dict[amp]]
                if LCA.count(most_frequent(LCA)) > 0.5*len(LCA):
                    amp_lca.loc[amp, canonical_groups[level]] = most_frequent(LCA)
    pickle.dump(amp_lca, open('{}antismashdb_amp_lca.p'.format(outputFolder), 'wb'))
    return amp_lca

def make_upset_product_subgroups(product_dict, isamplicon_matches, argOptions):
    '''
    Makes the groups of compounds, here it's set up to use AMP domain so it looks only for NRPS hybrids  
    :param product_dict: 
    :param isamplicon_matches: 
    :param argOptions: 
    :return: 
    '''

    valid_groups = ['nrps', 't1pks', 'terpene', 't2pks', 't3pks', 'hybrids']
    product_amplicons = {}

    for product_group in set(product_dict.values()):
        subcomponents = set(product_group.split('-'))
        #for each existing product combination
        #if it it's one of the most relevant combinations
        product_clusters = []
        if len(set(valid_groups).intersection(subcomponents)) == len(subcomponents):
            for bgc in product_dict.keys():
                #for each bgc
                if product_dict[bgc] == product_group:
                    #if it's product == to the combination
                    #get all domains associated with the bgc
                    subclusters = [x for x in isamplicon_matches if x.startswith(bgc)]
                    for cluster in subclusters:
                        for domain in isamplicon_matches[cluster].keys():
                            for match in isamplicon_matches[cluster][domain]:
                                product_clusters.append(match)
        if (product_clusters):
            product_amplicons[product_group] = set(product_clusters)

    amp_gbk_dict = {}
    for cluster in isamplicon_matches.keys():
        for domain in isamplicon_matches[cluster]:
            for match in isamplicon_matches[cluster][domain]:
                if match in amp_gbk_dict.keys():
                    amp_gbk_dict[match].append(domain)
                else:
                    amp_gbk_dict[match] = [domain]
    pickle.dump(amp_gbk_dict, open('{}amp_gbk_dict.p'.format(argOptions['outputDir']) , 'wb'))
    return product_amplicons, amp_gbk_dict

def create_cooccurrence_network(rare_table, outputFolder, argOptions):
    '''calculates pairwise occurrence between all features that appear in at least 3 replicates. 
    Uses pearsonr of presence absence patterns across samples'''
    bool_df = rare_table.astype(bool)*1
    abundant_features = bool_df.sum(axis=1)
    abundant_features_index = abundant_features[abundant_features>=argOptions['minReplicates']].index
    bool_df = bool_df.loc[abundant_features_index]
    correlation_df = pd.DataFrame(0, index=abundant_features_index, columns=abundant_features_index)
    pval_df = pd.DataFrame(0, index=abundant_features_index, columns=abundant_features_index)
    length_index = len(abundant_features_index)
    for id1 in range(length_index):
        for id2 in range(id1+1, len(abundant_features_index)):
            correlation_df.iloc[id1, id2], pval_df.iloc[id1,id2] = scipy.stats.pearsonr(bool_df.iloc[id1].values, bool_df.iloc[id2].values)
            correlation_df.iloc[id1, id2] = round(correlation_df.iloc[id1, id2], 4)
        print(f"Feature \r{id1} out of {length_index} done", sep = "", end = "", flush = True)
    pval_df.to_csv('{}{}_pval_abuntant_amps.csv'.format(outputFolder, argOptions['minReplicates']), sep='\t')
    correlation_df.to_csv('{}{}_correlation_abuntant_amps.csv'.format(outputFolder, argOptions['minReplicates']), sep='\t')
    return correlation_df

def alpha_diversity_plot(df, argOptions):
    '''makes alpha diversity plot. in the future can be extended to include metadata'''
    alpha_df = pd.DataFrame(0, index=df.columns, columns=['sample','shannon', 'amplicon_counts', 'simpson_e'])
    for s in df.columns:
        alpha_df.loc[s, 'amplicon_counts'] = alpha_diversity(metric='observed_otus', counts=df.loc[:, s])[0]
        alpha_df.loc[s, 'shannon'] = alpha_diversity(metric='shannon', counts=df.loc[:, s])[0]
        alpha_df.loc[s, 'simpson_e'] = alpha_diversity(metric='simpson_e', counts=df.loc[:, s])[0]
        alpha_df.loc[s, 'sample'] = s.split(argOptions['separator'])[0]

    fig = sns.barplot(alpha_df.loc[:, 'shannon'], orient='h')
    fig.axes.tick_params(axis='y', labelsize=6)
    plt.savefig('{}shannon_diversity.pdf'.format(argOptions['outputDir']))
    plt.close()

    fig = sns.barplot(alpha_df.loc[:, 'amplicon_counts'], orient='h')
    fig.axes.tick_params(axis='y', labelsize=6)
    plt.savefig('{}amplicon_counts.pdf'.format(argOptions['outputDir']))
    plt.close()

    fig = sns.barplot(alpha_df.loc[:, 'simpson_e'], orient='h')
    fig.axes.tick_params(axis='y', labelsize=6)
    plt.savefig('{}simpson_e.pdf'.format(argOptions['outputDir']))
    plt.close()

    fig = sns.boxplot(data=alpha_df, x='sample', y='amplicon_counts')
    fig.axes.tick_params(axis='y', labelsize=6)
    plt.savefig('{}amplicon_counts_swarmplot.pdf'.format(argOptions['outputDir']))
    plt.close()

def beta_diversity_plot(normalized_df, argOptions):
    '''makes beta diversity plots, in the future it can be extended to include metadata'''
    if argOptions['phylogenyTree']:
        tree = TreeNode.read(StringIO(open(argOptions['phylogenyTree']).read()))
        _ids = []
        for node in tree.tips():
            _ids.append(node.name)
        beta_div = beta_diversity(metric='unweighted_unifrac', counts=feature_table.transpose(), ids=feature_table.columns, otu_ids=_ids, tree=tree, validate=False)
    else:
        beta_div = beta_diversity(metric='braycurtis', counts=feature_table.transpose(), ids=feature_table.columns, validate=False)

    pcoa_results = pcoa(beta_div)
    if argOptions['sampleNames']:
        samples = argOptions['sampleNames'].split(',')
        meta_df = pd.DataFrame(0, index=normalized_df.columns, columns=['Sample'])
        for c in list(normalized_df.columns):
            for s in samples:
                if c.startswith(s):
                    meta_df.loc[c, ['Sample']] = s
        fig = pcoa_results.plot(df=meta_df, column='Sample', cmap='Set1', s=50)
        plt.savefig('{}beta_diversity_pcoa.pdf'.format(argOptions['outputDir']))
        plt.close()

        pca = PCA(n_components=2)
        X = feature_table.transpose().values
        X = StandardScaler().fit_transform(X)
        principalComponents = pca.fit_transform(X)
        principalDf = pd.DataFrame(data=principalComponents, columns=['principal component 1', 'principal component 2'])
        finalDf = pd.concat([principalDf, meta_df.reset_index().loc[:, 'Sample']], axis=1)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel('Principal Component 1', fontsize=15)
        ax.set_ylabel('Principal Component 2', fontsize=15)
        ax.set_title('2 component PCA', fontsize=20)
        targets = argOptions['sampleNames'].split(',')
        for target, color in zip(targets, sns.color_palette("husl", len(targets))):
            indicesToKeep = finalDf['Sample'] == target
            ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1'],
                       finalDf.loc[indicesToKeep, 'principal component 2'], c=np.array([color]), s=50)
            ax.legend(targets)
            ax.grid()
        plt.savefig('{}beta_div_pca.pdf'.format(argOptions['outputDir']))
        plt.close()

        mds = MDS(n_components=2, dissimilarity='precomputed')
        X_r = mds.fit_transform(beta_div.data)
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlabel('MDS 1', fontsize=15)
        ax.set_ylabel('MDS 2', fontsize=15)
        ax.set_title('MDS plot', fontsize=20)
        targets = argOptions['sampleNames'].split(',')
        for target, color in zip(targets, sns.color_palette("husl", len(targets))):
            indicesToKeep = finalDf[finalDf['Sample'] == target]
            for i in indicesToKeep.index:
                ax.scatter(X_r[i, 0], X_r[i, 1], c=np.array([color]), s=50)
            ax.legend(targets)
            ax.grid()
        plt.savefig('{}beta_diversity_mds.pdf'.format(argOptions['outputDir']))

    beta_div.plot()
    plt.savefig('{}beta_div.pdf'.format(argOptions['outputDir']))


def correlation_table_to_network(corr_network_annot, outputFolder, argOptions):
    import networkx
    from sklearn.cluster import KMeans, DBSCAN, AffinityPropagation

    #create network
    graph = networkx.convert_matrix.from_pandas_edgelist(corr_network_annot.reset_index(), source='amplicon_id1', target ='amplicon_id2', edge_attr='correlation')
    nw_matrix = networkx.convert_matrix.to_numpy_array(graph, weight='correlation')
    dist_matrix = np.subtract(nw_matrix, np.ones((len(nw_matrix), len(nw_matrix)))) * -1
    kmeans = KMeans(n_clusters=len(argOptions['sampleNames'].split(',')), init='k-means++', max_iter=300, n_init=10,
                    random_state=0)
    aff = AffinityPropagation()

    pred_aff = aff.fit_predict(dist_matrix)
    pred_y = kmeans.fit_predict(dist_matrix)

    cluster_annot_txt = 'node_id\tkmeans\tcluster_affinity_propagation\n'
    cluster_annot_dict = {}

    for x in range(len(pred_y)):
        cluster_annot_txt += '{}\t{}\t{}\n'.format(list(graph.nodes)[x], pred_y[x], pred_aff[x])
        if pred_y[x] in cluster_annot_dict:
            cluster_annot_dict[pred_y[x]].append(list(graph.nodes)[x])
        else:
            cluster_annot_dict[pred_y[x]] = [list(graph.nodes)[x]]

    outfile = open('{}network_clustering.csv'.format(outputFolder), 'w')
    outfile.write(cluster_annot_txt)
    outfile.close()
    clustering_results = pd.read_csv('{}network_clustering.csv'.format(outputFolder), sep='\t', index_col=0)
    return clustering_results

def network_output_tables(corr_network_annot, node_annot, kmeans_clusters, amp_gbk_dict, outputFolder):
    import networkx

    graph = networkx.convert_matrix.from_pandas_edgelist(corr_network_annot.reset_index(), source='amplicon_id1', target='amplicon_id2', edge_attr='correlation')
    for cluster in set(kmeans_clusters.loc[:, 'kmeans']):
        if cluster >= 0:

            header = 'antismashdb_cluster\tcooccurring_domains\tcluster_amplicons\tputative_bgc_taxonomy\n';
            separator = '\n-------------\nbgc\tamplicons\n'
            out_txt = ''; multi_out_txt = ''; appendix_dom = ''; multi_appendix_dom = ''
            nodes = []; nodes_plus_neigh = set([])
            cluster_tab = kmeans_clusters.loc[kmeans_clusters.loc[:, 'kmeans'] == cluster]
            cluster_index = cluster_tab.index

            for node in cluster_index:
                nodes_plus_neigh = nodes_plus_neigh.union(set(list(graph.neighbors(node))))
                nodes.append(node)

            asdb_dom_matches = [amp_gbk_dict[x] for x in list(nodes_plus_neigh) if x in amp_gbk_dict.keys()]
            amps_with_matches = [x for x in list(nodes_plus_neigh) if x in amp_gbk_dict.keys()]
            asdb_list = (list(set([item for sublist in asdb_dom_matches for item in sublist])))
            cluster_dict = {}

            for dom in asdb_list:
                for amp1 in range(len(amps_with_matches)):
                    for amp2 in range(amp1+1, len(amps_with_matches)):
                        if dom in amp_gbk_dict[amps_with_matches[amp1]] and dom in amp_gbk_dict[amps_with_matches[amp2]]:
                            bgc = '_'.join(dom.split('_')[:-1])
                            if bgc in cluster_dict.keys():
                                if dom in cluster_dict[bgc]:
                                    cluster_dict[bgc][dom] = cluster_dict[bgc][dom].union(set([amps_with_matches[amp1], amps_with_matches[amp2]]))
                                else:
                                    cluster_dict[bgc][dom] = set([amps_with_matches[amp1], amps_with_matches[amp2]])
                            else:
                                cluster_dict[bgc] = {dom : set([amps_with_matches[amp1], amps_with_matches[amp2]])}

            for bgc in cluster_dict.keys():
                domains = []
                for dom in cluster_dict[bgc].keys():
                    domains.extend(list(cluster_dict[bgc][dom]))
                appendix_dom += '{}\t{}\n'.format(bgc, '-'.join(list(set(domains))))
                taxas = (node_annot.loc[domains, ['genus']])
                taxa = (list(set([item for sublist in taxas.values for item in sublist])))
                out_txt += '{}\t{}\t{}\t{}\n'.format(bgc, len(cluster_dict[bgc]), len(set(domains)), '-'.join(taxa))

                if len(cluster_dict[bgc])>1:
                    if '-'.join(sorted(list(set(domains)))) not in multi_appendix_dom:
                        multi_out_txt += '{}\t{}\t{}\t{}\n'.format(bgc, len(cluster_dict[bgc]), len(set(domains)), ','.join(taxa))
                        multi_appendix_dom += '{}\t{}\n'.format(bgc, '-'.join(sorted(list(set(domains)))))
            outfile = open('{}subnetwork_{}_putative_clusters.txt'.format(outputFolder, cluster), 'w')
            outfile.write(header+out_txt+separator+appendix_dom)
            outfile.close()

            outfile = open('{}subnetwork_{}_multiple_domains_putative_clusters.txt'.format(outputFolder, cluster), 'w')
            outfile.write(header+multi_out_txt+separator+multi_appendix_dom)
            outfile.close()

if __name__ == "__main__":
    argOptions = commandLineParser()
    if argOptions['verbose'] == True:
        print ('loading feature_table')
    feature_table = pd.read_csv(argOptions['featureTable'], sep='\t', index_col=0, header=1)
    print (feature_table)
    print (feature_table.columns)

    #Load amplicons from nAMPs and in silico amps from reference dbs

    if argOptions['verbose'] == True:
        print ('loading amplicons')
    soil_amplicons_list, mibig_amplicons_list, antismashdb_amplicons_list, metagenome_amplicons_list = load_amplicons(feature_table.index, argOptions)    #with all the amplicons

    #sets output folders and creates the directories if they don't exist already
    outputFolder = argOptions['outputDir']
    mibigOutputFolder = argOptions['MiBIGOutputDir']
    metagenomeFolder = argOptions['metagenomeOutputDir']
    antismashdbOutputFolder = argOptions['antismashOutputDir']


    if not os.path.isdir(outputFolder):
        if argOptions['verbose'] == True:
            print ('Creating output folders')
        os.makedirs(outputFolder)

    if not os.path.isdir(mibigOutputFolder):
        os.makedirs(mibigOutputFolder)

    if not os.path.isdir(metagenomeFolder):
        os.makedirs(metagenomeFolder)

    if not os.path.isdir(antismashdbOutputFolder):
        os.makedirs(antismashdbOutputFolder)


    #groups amplicons that share the same sequence in amplicons clusters

    if argOptions['verbose'] == True:
        print ('Grouping similar amplicons')
    total_ampli = len(soil_amplicons_list)
    if total_ampli > 4500:
        print("Too many amplicons detected, the program will divide them by group of 4500") # good performance for 4500 amplicons
        n = 4500
        if not os.path.isdir(f"{argOptions['outputDir']}amplicon_processing"):
            os.makedirs(f"{argOptions['outputDir']}amplicon_processing")
        if not os.path.isdir(f"{argOptions['outputDir']}amplicon_processing/sub_amplicon"):
            os.makedirs(f"{argOptions['outputDir']}amplicon_processing/sub_amplicon")
        if not os.path.isdir(f"{argOptions['outputDir']}amplicon_processing/output/"):
            os.makedirs(f"{argOptions['outputDir']}amplicon_processing/output/")
        # creation of texte file with sub amplicons
        if not os.path.isfile(f"{argOptions['outputDir']}amplicon_processing/sub_amplicon/amplicon_1.faa"):
            print(f"{argOptions['outputDir']}amplicon_processing/sub_amplicon/amplicon_1.faa")
            divide_amplicon_list(f"{argOptions['amplicons']}", total_ampli // n)
        New_amplicon_list = [soil_amplicons_list[i * n:(i + 1) * n] for i in range((total_ampli + n - 1) // n )] # minus 1 to finish the file for the last group
        sub_cluster(New_amplicon_list, feature_table) # clustering for sub amplicons
        print("Regrouping the sub clusters in a dictionnary...")
        if not os.path.isfile(f"{argOptions['outputDir']}amplicon_processing/output/dico_cluster.p"):
            dict_cluster = merge_cluster(f"{argOptions['outputDir']}amplicon_processing/output/", len(New_amplicon_list)) # merge all the clusterin one dictionnary
        else:
            dict_cluster = pickle.load(open(f"{argOptions['outputDir']}amplicon_processing/output/dico_cluster.p", 'rb'))       
        if not os.path.isfile(f"{argOptions['outputDir']}grouped_amps.tsv"):
            merge_sequence(dict_cluster) # creation of the final sequence table
        print("Final sequence table created")
        print("Creation of amplicon_list.p file...")
        if not os.path.isfile(f"{argOptions['outputDir']}amplicon_list.p"):
            merge_amplicon(dict_cluster)
        print("amplicon_list.p file created")
        print("Creation of the grouped_feature_table.tsv file...")
        if not os.path.isfile(f"{argOptions['outputDir']}grouped-feature-table.tsv"):
            merge_feature(dict_cluster, len(New_amplicon_list))
        print("The feature table has been generated.\nThe normal program can now be executed.")
        
        feature_table = pd.read_csv(f"{argOptions['outputDir']}grouped-feature-table.tsv", sep = '\t', index_col = 0, header = 1)
        soil_amplicons_list = pickle.load(open(f"{argOptions['outputDir']}amplicon_list.p", 'rb'))
    else:
        feature_table, soil_amplicons_list = group_similar_amplicons(feature_table, soil_amplicons_list, argOptions)
    #rarefies the table

    if argOptions['verbose'] == True:
        print ('Creating rarefaction table')
    rare_table = make_rare_table(feature_table, argOptions['rarefactionThreshold'])
    print (rare_table)
    #groups MIBIG domains to BGC
    if argOptions['verbose'] == True:
        print ('Grouping MiBIG in silico domains belonging to the same BGC')
    if not os.path.isfile('{}mibig_bgc_dict.p'.format(mibigOutputFolder)):
        mibig_BGC_dict = group_insilico_amplicons(mibig_amplicons_list)
        pickle.dump(mibig_BGC_dict, open('{}mibig_bgc_dict.p'.format(mibigOutputFolder), 'wb'))
    else:
        mibig_BGC_dict = pickle.load(open('{}mibig_bgc_dict.p'.format(mibigOutputFolder), 'rb'))

    #creates all the annotation for MiBIG using pairwise identity
    dbtype = 'mibig'
    if not os.path.isfile('{}{}_soil_annotation.p'.format(mibigOutputFolder, dbtype)):
        if argOptions['verbose'] == True:
            print ('Matching MiBIG in silico amplicons to nAMPs')
        mibig_domains_annotation, mibig_BGC_matches, mibig_isamplicon_matches = calculate_similarity(mibig_BGC_dict, soil_amplicons_list, mibig_amplicons_list, mibigOutputFolder, 'mibig')
    else:
        if argOptions['verbose'] == True:
            print ('Found pairwise matching of MiBIG in silico amplicons and nAMPs, loading results')
        mibig_domains_annotation = pickle.load(open('{}{}_soil_annotation.p'.format(mibigOutputFolder, dbtype), 'rb'))
        mibig_BGC_matches = pickle.load(open('{}{}_BGC_to_soil.p'.format(mibigOutputFolder, dbtype), 'rb'))
        mibig_isamplicon_matches = pickle.load(open('{}{}_bgc_extended_dict.p'.format(mibigOutputFolder, dbtype), 'rb'))

    if argOptions['metagenomeAmplicons']:
        # groups shotgun metagenome domains to BGC
        if argOptions['verbose'] == True:
            print ('Grouping metagenome derived in silico domains belonging to the same BGC')
        BGC_metagenome_dict = group_metagenome_amplicons(metagenome_amplicons_list)

        # creates all the annotation for the shotgun metagenome using pairwise identity
        dbtype = 'metagenome'
        if not os.path.isfile('{}{}_soil_annotation.p'.format(metagenomeFolder, dbtype)):
            if argOptions['verbose'] == True:
                print ('Matching metagenome derived in silico amplicons to nAMPs')
            metagenome_domains_annotation, metagenome_BGC_matches, metagenome_isamplicon_matches = calculate_similarity(BGC_metagenome_dict, soil_amplicons_list, metagenome_amplicons_list, metagenomeFolder, 'meta')
        else:
            if argOptions['verbose'] == True:
                print ('Found pairwise matching of metagenome derived in silico amplicons and nAMPs, loading results')
            metagenome_domains_annotation = pickle.load(open('{}{}_soil_annotation.p'.format(metagenomeFolder, dbtype), 'rb'))
            metagenome_BGC_matches = pickle.load(open('{}{}_BGC_to_soil.p'.format(metagenomeFolder, dbtype), 'rb'))
            metagenome_isamplicon_matches = pickle.load(open('{}{}_bgc_extended_dict.p'.format(metagenomeFolder, dbtype), 'rb'))

    # groups antismashdb domains to BGC
    if argOptions['verbose'] == True:
        print ('Grouping antismash-database derived in silico domains belonging to the same BGC')
    if not os.path.isfile('{}antismashdb_bgc_dict2.p'.format(antismashdbOutputFolder)):
        BGC_antismashdb_dict = group_insilico_amplicons(antismashdb_amplicons_list)
        pickle.dump(BGC_antismashdb_dict, open('{}antismashdb_bgc_dict2.p'.format(antismashdbOutputFolder), 'wb'))
    else:
        BGC_antismashdb_dict = pickle.load(open('{}antismashdb_bgc_dict2.p'.format(antismashdbOutputFolder), 'rb'))
    # creates all the annotation for antismashdb in silico amplicons using pairwise identity
    dbtype = 'antismashdb'
    if not os.path.isfile('{}{}_soil_annotation.p'.format(antismashdbOutputFolder, dbtype)):
        if argOptions['verbose'] == True:
            print ('Matching antismash-database derived in silico amplicons to nAMPs')
        domains_annotation, BGC_matches, isamplicon_matches = calculate_similarity(BGC_antismashdb_dict, soil_amplicons_list, antismashdb_amplicons_list, antismashdbOutputFolder, 'antismashdb')
    else:
        if argOptions['verbose'] == True:
            print ('Found pairwise matching of antimsash-database in silico amplicons and nAMPs, loading results')
        domains_annotation = pickle.load(open('{}{}_soil_annotation.p'.format(antismashdbOutputFolder, dbtype), 'rb'))
        BGC_matches = pickle.load(open('{}{}_BGC_to_soil.p'.format(antismashdbOutputFolder, dbtype), 'rb'))
        isamplicon_matches = pickle.load(open('{}{}_bgc_extended_dict.p'.format(antismashdbOutputFolder, dbtype), 'rb'))

    # assignment of the product annotation [gene cluster type]
    if not os.path.isfile('{}antismashdb_product_dict.p'.format(outputFolder)):
        if argOptions['verbose'] == True:
            print ('Antismash-database product dictionary not found. Generating [this may take some time]')
        product_dict = parse_antismashdb_bgc_family(argOptions)
    else:
        if argOptions['verbose'] == True:
            print ('Antismash-database product dictionary found')
        product_dict = pickle.load(open('{}antismashdb_product_dict.p'.format(outputFolder), 'rb'))

    #assignment of the taxonomy annotation
    if not os.path.isfile('{}antismashdb2_taxa_dict.p'.format(argOptions['antismashOutputDir'])):
        if argOptions['verbose'] == True:
            print ('Antismash-database taxonomy annotation file found. Generating [this will take some time]')
        amplicon_taxonomy = assign_taxonomy(BGC_antismashdb_dict, argOptions)
    else:
        if argOptions['verbose'] == True:
            print ('Antismash-database taxonomy dictionary found')
        amplicon_taxonomy = pickle.load(open('{}antismashdb2_taxa_dict.p'.format(argOptions['antismashOutputDir']), 'rb'))
        
    # part to assign taxonomy to amplicons based on LCA
    if not os.path.isfile('{}antismashdb_taxa_df.p'.format(outputFolder)):
        if argOptions['verbose'] == True:
            print ('Annotating nAMPs with taxonomy')
        taxa_df = parse_taxonomy(amplicon_taxonomy, outputFolder)
    else:
        if argOptions['verbose'] == True:
            print ('Found nAMPs taxonomy annotation')
        taxa_df = pickle.load(open('{}antismashdb_taxa_df.p'.format(outputFolder), 'rb'))

    #assigns the taxonomy to the amplicons
    if not os.path.isfile('{}antismashdb_amplicon_taxonomy_dict.p'.format(outputFolder)):
        if argOptions['verbose'] == True:
            print ('Annotating nAMPs with their product family')
        amplicon_taxonomy_dict = assign_taxonomy_to_amplicons(isamplicon_matches, taxa_df, outputFolder)
    else:
        if argOptions['verbose'] == True:
            print ('Found nAMPs product annotation')
        amplicon_taxonomy_dict = pickle.load(open('{}antismashdb_amplicon_taxonomy_dict.p'.format(outputFolder),'rb'))

    #uses ncbi to reconstruct the LCA of the hits
    if not os.path.isfile('{}antismashdb_amp_lca.p'.format(outputFolder)):
        if argOptions['verbose'] == True:
            print ('Annotating nAMPs taxonomy with LCA')
        amp_lca = assign_lca_taxonomy(amplicon_taxonomy_dict, outputFolder)
    else:
        if argOptions['verbose'] == True:
            print ('Found LCA taxonomy annotation of nAMPs taxonomy')
        amp_lca = pickle.load(open('{}antismashdb_amp_lca.p'.format(outputFolder), 'rb'))
    
    #groups amplicons in the different product groups [ex. pure NRPS, t1pks-nrps, terpene-nrps etc..]

    if argOptions['verbose'] == True:
        print ('Assigning nAMPs to product families')
    product_amplicons, gbk_amplicons = make_upset_product_subgroups(product_dict, isamplicon_matches, argOptions)

    #writes amplicon annotation to file
    if not os.path.isfile('{}amplicon_annotation.csv'.format(outputFolder)):
        if argOptions['verbose'] == True:
            print ('Collecting nAMPs annotation and writing to output files')
        write_annotation_to_csv(feature_table, product_amplicons, amp_lca)

    #creates alpha and beta diversity plots
    # if argOptions['verbose'] == True:
    #     print ('Calculating alpha diversity [richness] measurements and creating figures')
    # alpha_diversity_plot(rare_table, argOptions)
    # if argOptions['verbose'] == True:
    #     print ('Calculating beta diversity [composition] measurements and creating figures')
    # beta_diversity_plot(feature_table, argOptions)

    # creates co-occurrence network
    if not os.path.isfile('{}{}_correlation_abuntant_amps.csv'.format(outputFolder, argOptions['minReplicates'])):
        if argOptions['verbose'] == True:
            print ('Calculating pairwise occurrence of nAMPs. This may take some time')
        corr_df = create_cooccurrence_network(rare_table, outputFolder, argOptions)
    else:
        if argOptions['verbose'] == True:
            print ('Loading dataframe with pairwise occurrence of nAMPs')
        corr_df = pd.read_csv('{}{}_correlation_abuntant_amps.csv'.format(outputFolder, argOptions['minReplicates']), sep='\t', index_col=0, header=0).astype(float)
    print (corr_df)
    #correlation network and annotation
    if not os.path.isfile('{}{}_corr_network_annot_table.csv'.format(outputFolder, argOptions['networkThreshold'])) or not os.path.isfile('{}{}_node_annotation.csv'.format(outputFolder, argOptions['networkThreshold'])):
        if argOptions['verbose'] == True:
            print ('Creating network from correlation table and transfering the annotation to the nodes')
        node_annotation = filter_correlation_table_and_annotate(corr_df, product_amplicons, outputFolder, amp_lca, mibig_BGC_matches, gbk_amplicons, argOptions, argOptions['networkThreshold'])
        corr_network_annot = pd.read_csv('{}{}_corr_network_annot_table.csv'.format(outputFolder, argOptions['networkThreshold']), sep='\t', index_col=[0,1])
    else:
        if argOptions['verbose'] == True:
            print ('Network and annotation files found')
        corr_network_annot = pd.read_csv('{}{}_corr_network_annot_table.csv'.format(outputFolder, argOptions['networkThreshold']), sep='\t', index_col=[0,1])
        node_annotation = pd.read_csv('{}{}_node_annotation.csv'.format(outputFolder,  argOptions['networkThreshold']), sep='\t', index_col=0, header=0)

    #if not os.path.isfile('{}amplicons_matches.txt'.format(outputFolder)):
    #    write_amplicon_matches(argOptions, product_dict, amp_lca, mibig_BGC_matches, gbk_amplicons)

    if not os.path.isfile('{}network_clustering.csv'.format(outputFolder)):
        clustering_result = correlation_table_to_network(corr_network_annot, outputFolder, argOptions)
    else:
        clustering_result = pd.read_csv('{}network_clustering.csv'.format(outputFolder), sep='\t', index_col=0)

    network_output_tables(corr_network_annot, node_annotation, clustering_result, gbk_amplicons, outputFolder)

    if argOptions['verbose'] == True:
        print ('Analysis completed. Richness/community estimates can be visualized in the different plots. \nInterpretation of the correlation network can be done in cytoscape, load the annotation for additional information.')

