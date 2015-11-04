
"""
Script to select n motifs of the TF and m co-regulatory motifs
Rule for selection of co-regulatory motifs:
1) Motif is not included if: in the tf annotation at least one TF is in the list of "covered TFs"
2) If motif of the same cluster is in the list of "covered TFs"

"""

from os import listdir
from os.path import isfile, join
import optparse
import sys

parser = optparse.OptionParser()
parser.add_option("-a", "--motif_annotation_path", action="store", type="string", dest="motif_annotation_path", help='Name of the TF to process')
parser.add_option("-s", "--stat_tbl_path", action="store", type="string", dest="stat_tbl_path", help='Name of the TF to process')
parser.add_option("-m", "--motif_path", action="store", type="string", dest="motif_path", help='Name of the TF to process')
parser.add_option("-t", "--tf_name", action="store", type="string", dest="tf_name", help='Name of the TF to process')
(options, args) = parser.parse_args()



PATH_TO_ANNOTATION_TABLE = options.motif_annotation_path
PATH_TO_STAT_TBL_RESULTS = options.stat_tbl_path
PATH_TO_SINGLETONS = options.motif_path
TFNAME = options.tf_name

db_name='hg19-regions-220330-9species.extracted'
NES_thershold = 2.5
num_motifs_per_tf = 10
num_motifs_for_cofactors = 10


def Tf2Motif(path_to_annot_table):
    tf_motif_dict = dict()
    with open(path_to_annot_table,'r') as annotation_file:
        for line in annotation_file:
            line = line.split("\t")
            TF_name = line[2]
            motif_name = line[0]
            if TF_name not in tf_motif_dict:
                tf_motif_dict[TF_name] = (motif_name)
            else:
                cur_values = tf_motif_dict[TF_name]
                cur_values = ",".join((cur_values, motif_name))
                tf_motif_dict[TF_name] = cur_values



def Motif2Tf(path_to_annot_table):
    motif2tf_dict = dict()
    with open(path_to_annot_table, 'r') as annotation_file:
        for line in annotation_file:
            line = line.split("\t")
            TF_name = line[2]
            motif_name = line[0]
            if motif_name not in motif2tf_dict:
                motif2tf_dict[motif_name] = (TF_name)
            else:
                cur_values = motif2tf_dict[motif_name]
                cur_values = ",".join((cur_values, TF_name))
                motif2tf_dict[motif_name] = cur_values
    return motif2tf_dict



def read_motf_collection_folder(path_to_folder):
    existing_motifs_in_collection = dict()
    onlyfiles = [ f for f in listdir(path_to_folder + "/") if isfile(join(path_to_folder,f)) ]
    for elem in onlyfiles:
        existing_motifs_in_collection[elem] = elem
    return existing_motifs_in_collection


def read_motf_collection_file_list(path_to_file):
    existing_motifs_in_collection = dict()
    with open(path_to_file,'r') as my_file:
        for line in my_file:
            line = line.split()
            existing_motifs_in_collection[line[0]] = line[0]
    return existing_motifs_in_collection


def read_cluster_tbl(path_to_cluster_tbl):
    cluster_motif_name = dict()
    with open(path_to_cluster_tbl, 'r') as my_file:
        for line in my_file:
            line = line.split()
            motif_name = line[2]
            cluster_label = line[0]
            cluster_motif_name[motif_name] = cluster_label
    return cluster_motif_name




def FindEnrichedMotifs_Linked_with_TF( path_to_stat_tbl, path_to_annot_table, path_to_singleton_folder, tf_name, extension = '.cb' ):

    """
    full_path_to_stat_tbl - path to statistics.tbl file from cluster buster
    all_motifs_list - motifs from singleton folder
    tf_name - name of transcription factor which we are processing now
    """
    motif2tf_dict = Motif2Tf( path_to_annot_table )
    motif_list = read_motf_collection_file_list(path_to_singleton_folder)
    singletone_list = []
    with open(path_to_stat_tbl, 'r') as stat_tbl_file:
        next(stat_tbl_file)
        for line in stat_tbl_file:
            line = line.split("\t")
            nes_score = float(line[7])
            # if nes_score>= NES_thershold:
            motif_name = line[2]
            if line[5] == db_name:
                if motif_name + extension in motif_list:
                    if motif_name in motif2tf_dict:
                        TFs = motif2tf_dict[motif_name]
                        # TFs = line[4]
                        TFs = TFs.split(",")
                        if tf_name in TFs:
                            singletone_list.append(motif_name + '.cb')
    return singletone_list




def FindEnrichedMotifs(path_to_stat_tbl, path_to_annot_table, path_to_singleton_folder, tf_name, extension = '.cb'):
    """
    full_path_to_stat_tbl - path to statistics.tbl file from cluster buster
    all_motifs_list - motifs from singleton folder
    tf_name - name of transcription factor which we are processing now
    """

    motif2tf_dict = Motif2Tf( path_to_annot_table )
    motifname_NES_score_dict = dict()
    motif_list = read_motf_collection_file_list(path_to_singleton_folder)
    enriched_cb_list = []
    with open(path_to_stat_tbl, 'r') as stat_tbl_file:
        for line in stat_tbl_file:
            line = line.split("\t")
            motif_name = line[2]
            nes_score = line[7]
            new_motif_id = motif_name + ".cb"
            motifname_NES_score_dict[new_motif_id] = nes_score
            if motif_name + extension in motif_list:
                if motif_name in motif2tf_dict:
                    TFs = motif2tf_dict[motif_name]
                    TFs = TFs.split(",")
                    if tf_name not in TFs:
                        enriched_cb_list.append(motif_name + '.cb')
    return enriched_cb_list





"""
Selects motifs per TF (co-regulatory motifs)
"""


def select_EnrichedMotifs_per_TF_uniq(path_to_annot_table, enriched_cb_list):
    motif2tf_dict = Motif2Tf( path_to_annot_table )
    TF_cluster_list = []
    selected_motifs=[]
    for motif in enriched_cb_list:
        if motif[:-3] in motif2tf_dict:
            tf_for_motif = motif2tf_dict[motif[:-3]]
            tf_for_motif = tf_for_motif.split(',')
            not_in_tf_list = [i for i in tf_for_motif if i not in TF_cluster_list]
            for tf in tf_for_motif:
                if tf not in TF_cluster_list:
                    TF_cluster_list.append(tf)
                    if motif not in selected_motifs:
                        selected_motifs.append(motif)
    return selected_motifs








def AnnotEnriched_merge_singletons(singletone_list, enriched_cb_list):
    motifs_to_use = singletone_list[:num_motifs_per_tf] + enriched_cb_list[:num_motifs_for_cofactors]
    return motifs_to_use



def find_co_regulatory_tfs(path_to_singleton_folder, path_to_annot_table, path_to_stat_tbl, path_to_cluster_tbl, forbiden_tfs_as_cofactors, forbiden_clusters, extension = '.cb'):

    """
    full_path_to_stat_tbl - path to statistics.tbl file from cluster buster
    all_motifs_list - motifs from singleton folder
    tf_name - name of transcription factor which we are processing now
    """

    motif2tf_dict = Motif2Tf( path_to_annot_table )
    motifname_NES_score_dict = dict()
    motif_list = read_motf_collection_file_list(path_to_singleton_folder)
    enriched_cb_list = []
    cluster2TF = read_cluster_tbl(path_to_cluster_tbl)
    with open(path_to_stat_tbl, 'r') as stat_tbl_file:
        for line in stat_tbl_file:
            line = line.split("\t")
            motif_name = line[2]
            nes_score = line[7]
            if nes_score>=NES_thershold:
                new_motif_id = motif_name + ".cb"
                motifname_NES_score_dict[new_motif_id] = nes_score
                if motif_name in  cluster2TF:
                    cluster_number = cluster2TF[motif_name]
                    if motif_name + extension in motif_list:
                        if motif_name in motif2tf_dict:
                            TFs = motif2tf_dict[motif_name]
                            # TFs = line[4]
                            to_select=False
                            TFs = TFs.split(",")
                            for i in range(0, len(TFs)):
                                tf = TFs[i]
                                if tf not in forbiden_tfs_as_cofactors:
                                    if cluster_number not in forbiden_clusters:
                                        if i==len(TFs)-1:
                                            to_select=True
                            if to_select:
                                enriched_cb_list.append(motif_name + '.cb')
                                ###---Add all TFs to forbiden set---###
                                for tf in TFs:
                                    forbiden_tfs_as_cofactors[tf]=tf
                                ###---Add all clusters to list of forbiden clusters---###
                                forbiden_clusters[cluster_number] = cluster_number
    return enriched_cb_list


"""
"""


tf_motifs_list = FindEnrichedMotifs_Linked_with_TF( PATH_TO_STAT_TBL_RESULTS, PATH_TO_ANNOTATION_TABLE, PATH_TO_SINGLETONS, TFNAME )

enriched_cb_list = FindEnrichedMotifs(PATH_TO_STAT_TBL_RESULTS, PATH_TO_ANNOTATION_TABLE, PATH_TO_SINGLETONS, TFNAME, extension = '.cb')

selected_motifs = select_EnrichedMotifs_per_TF_uniq(PATH_TO_ANNOTATION_TABLE, enriched_cb_list)


if len(tf_motifs_list)>num_motifs_per_tf:
    pass
else:
    num_motifs_per_tf=len(tf_motifs_list)



if len(selected_motifs)>num_motifs_for_cofactors:
    pass
else:
    num_motifs_for_cofactors=len(selected_motifs)

for i in range(0, num_motifs_per_tf):
    sys.stdout.write(tf_motifs_list[i] + "\t" + "tf_motifs" + "\n")
for i in range(0,num_motifs_for_cofactors):
    sys.stdout.write(selected_motifs[i] + "\t" + "co_regulatory_tfs" + "\n")


