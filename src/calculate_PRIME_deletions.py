"""
Input file example:
chr10	100028287	G	GCA
chr10	101292635	G	GC
chr10	102027666	T	TCTAAC
chr10	102027666	T	TCTAACCTAAC
chr10	102027666	T	TCTAACCTACC
"""

import sys
import subprocess
import bx.seq.twobit
from array import *
from collections import OrderedDict
from collections import defaultdict
from sklearn.externals import joblib
import numpy as np
import optparse


def select_fasta_for_region(chrID_fasta_dict, region_bed_format_list):
    chr_name = region_bed_format_list[0]
    start = int(region_bed_format_list[1])
    end = int(region_bed_format_list[2])
    fasta_sequence=chrID_fasta_dict[chr_name].get(start,end).upper()
    fasta_sequence = array('c',fasta_sequence)
    return fasta_sequence



def k_sliding_windows2fasta_deletion( chrom, start_indel, del_descr, sw_size ):

    regionID_fasta_dict = OrderedDict()
    regionID_mut_fasta_dict = OrderedDict()
    """
    regionID_fasta_dict - key: region ID (chr:start-end; value: fasta sequence)
    """
    ###----define start of the sliding window----###
    mut_start_chr = chrom
    mut_start = int(start_indel)
    ###----open file to save results-----###
    shift = int(sw_size/10)
    ##---iterate n_tmes=sw size----###
    for i in range( 0, sw_size+1, shift ):
        if i!=0:
            start_sw = mut_start - sw_size + i
            end_sw = start_sw + sw_size + 1
            chr_start_end_list= [mut_start_chr, start_sw, end_sw]
            ###----select fasta sequence for the region------###
            ###---find mutation to introduce----###
            ref_nuc = del_descr.split(":")[0][0]
            deleted_seq =  del_descr.split(":")[0][1:]
            mut_start_local_coord = mut_start - start_sw
            deletion_description = del_descr.split(":")[0]
            ###---Make predictions for reference sequence---###
            regionID = mut_start_chr + ":" + str(start_sw) + "-" + str(end_sw) + ":" + deletion_description + ":" + ref_nuc
            fasta_sequence = select_fasta_for_region(chrID_fasta_dict, chr_start_end_list)
            ###---Mutate fasta---###
            if fasta_sequence[mut_start_local_coord] == ref_nuc:
                fasta_sequence_new = fasta_sequence.tostring()
                if regionID not in regionID_fasta_dict:
                    regionID_fasta_dict[regionID] = fasta_sequence_new
            else:
                print "Reference nucleotide differs from hg19 version", mut_start_chr, mut_start ,ref_nuc, deleted_seq
                fasta_sequence_new = fasta_sequence.tostring()
                if regionID not in regionID_fasta_dict:
                    regionID_fasta_dict[regionID] = fasta_sequence_new

            ###----Make insertion to the fasta-----###
            if fasta_sequence[mut_start_local_coord] == ref_nuc:
                fasta_sequence_new = fasta_sequence.tostring()
                ###---make deletion---###
                insertion_start = mut_start_local_coord + 1
                end_of_delet_local_coord = insertion_start + len(deleted_seq)
                fasta_sequence_insertion = fasta_sequence_new[:insertion_start] + fasta_sequence_new[end_of_delet_local_coord:]
                if regionID not in regionID_mut_fasta_dict:
                    regionID_mut_fasta_dict[regionID] = fasta_sequence_insertion
            else:
                print "Reference nucleotide differs from hg19 version", mut_start_chr, mut_start ,ref_nuc, deleted_seq
                insertion_start = mut_start_local_coord+1
                fasta_sequence_insertion = fasta_sequence_new[:insertion_start] + fasta_sequence_new[end_of_delet_local_coord:]
                if regionID not in regionID_mut_fasta_dict:
                    regionID_mut_fasta_dict[regionID] = fasta_sequence_insertion
    return regionID_fasta_dict, regionID_mut_fasta_dict



def check_allowed_symbols(string):
    allowed_symbols = 'ACTG'
    for symbol in string:
        if symbol not in allowed_symbols:
            return False
    return True


def read_file_with_indels(path_to_file):
    allowed_symbols = 'ACTG'
    coord_muttype_dict = OrderedDict()
    with open(path_to_file) as my_file:
        for line in my_file:
            line = line.split()
            chrom = line[0]
            start = line[1]
            ###----check that all symbols are from ACTG alphabet-----###
            if not check_allowed_symbols(line[2]):
                continue
            if not check_allowed_symbols(line[3]):
                continue
            mut_type = line[2] + ":" + line[3]
            try:
                region_mutID = chrom + "_" + start + "_" + mut_type + "_" + line[4]
            except IndexError:
                print "No columne (#4) with ID"
                sys.exit(1)

            if region_mutID in  coord_muttype_dict:
                print region_mutID, " not unique"
            coord_muttype_dict[region_mutID] = mut_type
    return coord_muttype_dict


def saveFastaFile( path_to_save_fasta_file, seqID_fasta_dict ):
    file_to_save=open( path_to_save_fasta_file, "w" )
    for key in seqID_fasta_dict:
        file_to_save.write(">" + key + "\n")
        file_to_save.write(seqID_fasta_dict[key] + "\n")
    file_to_save.close()



def merge2oneFileSingletonScanWG(Text, oldText):
    """
    Output is in format:
    chr1    247279147    247282336    transfac_pro-M00801    136.0
    """
    Merged_file_RAM = ''
    ####Parsing of stdout####
    Text = Text.split("\n")
    motif_name = ''
    for line in Text:
        chrom = ''
        start = ''
        end = ''
        crm_score = ''
        if line.startswith("#"):
            tmp_line = line.split()
            motif_name = tmp_line[5]
        line = line.split()
        if len(line) == 5 and str(line[3]).startswith('chr'):
            chrom = str(line[3])   ###chrom:start-end
            start = str(line[1])
            end = str(line[2])
            crm_score = str(float(line[0]))
            ###----modify chrom output----###
            region_start_end = chrom.replace(":","\t")
            region_start_end = region_start_end.replace("-","\t")
            region_start_end = region_start_end.split("\t")
            chr_start = int(region_start_end[1]) + int(start)
            chr_end = int(region_start_end[1]) + int(end)
            chr_start_end = region_start_end[0] + "\t" + str(chr_start) + "\t" + str(chr_end)
        if chrom != '' and str(start) != '' and str(end) != '' and motif_name != '' and str(crm_score) != '':
            Merged_file_RAM += chr_start_end + "\t" + motif_name + "\t" + str(crm_score) + "\n"
            Merged_file_RAM += chrom + "\t" + str(chr_start)  + "\t" + str(chr_end) + "\t" + motif_name + "\t" + str(crm_score) + "\n"
    ###---Merge perocessed results with existing---###
    oldText = oldText + Merged_file_RAM
    return oldText


def run_cmd(cmd, stdin_str=None):
    """
    Run the program with the provided arguments.

      - When the program succeeds, return stdout and stderr.
      - When the program fails, print stdout and stderr messages.
    """

    try:
        pid = subprocess.Popen(args=cmd, stdin=subprocess.PIPE , stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        stdoutdata, stderrdata = pid.communicate(input=stdin_str)
    except OSError, msg:
        print >> sys.stderr, "\nERROR during execution of '" + ' '.join(cmd) + "': " + str(msg)
        sys.exit(1)

    if (pid.returncode != 0):
        print >> sys.stderr, "\nERROR during execution of '" + ' '.join(cmd) + "':"
        print >> sys.stderr, "\nStandard output:\n----------------"
        print >> sys.stderr, stdoutdata
        print >> sys.stderr, "\nStandard error:\n---------------"
        print >> sys.stderr, stderrdata
        sys.exit(1)

    return [stdoutdata, stderrdata]

def run_cbust_scan_pipe( singletone_list, modified_fasta_seq ):
    ###---Variable to save results of scanning---###
    initScanningData=''
    for motif in singletone_list:
        stdoutdata, stderrdata = run_cmd([PATH_TO_CBUST, '-c', '0', '-m', '0', '-f', '3', PATH_TO_FOLDER_WITH_SINGLETONS + motif], modified_fasta_seq)
        initScanningData = merge2oneFileSingletonScanWG( stdoutdata, initScanningData )
    return initScanningData



###----functions----###
def readFetureOrder(path_to_file):
    FeatureOrder_list=[]
    with open(path_to_file) as my_file:
        for line in my_file:
            motif_name = line.split()[0]
            FeatureOrder_list.append(motif_name)
    return FeatureOrder_list


def CRMmaxScore_Homo_Cluster(Merged_file_RAM, motif_order, regionOrder):
    """
    input file should have the format:
    chr16:262550-262750    262678    262697    transfac_pro-M01652    0.747
    chr16:262550-262750    262678    262697    transfac_pro-M01652    0.747
    chr16:262550-262750    262678    262697    transfac_pro-M01652    0.747
    """
    MaxScoredict = defaultdict(dict)
    ###---put regions to the dictionary:
    reg_idx=dict()
    i=0
    for region_id in regionOrder:
        reg_idx[region_id] = i
        i=i+1
    motif_idx=dict()
    i=0
    for motif in motif_order:
        ###----remove cb extension----###
        motif_idx[motif[:-3]]=i
        i=i+1
    ###-------------------------------###
    ####-----Read merged results for cbust scoring------####
    Merged_file_RAM = Merged_file_RAM.split("\n")
    for line in Merged_file_RAM:
        line=line.split()
        if len(line) ==5:
            regionID=line[0]
            motifID=line[3]
            score=float(line[4])
            if regionID not in MaxScoredict:
                MaxScoredict[regionID][motifID]=score
            else:
                if motifID not in MaxScoredict[regionID]:
                    MaxScoredict[regionID][motifID]=score
                else:
                    if MaxScoredict[regionID][motifID]<score:
                        MaxScoredict[regionID][motifID]=score
    ###---initialize the matrix with values
    featureMatrix_POS = np.zeros((len(regionOrder), len(motif_order)))
    row_idx=0
    col_idx=0
    for region in reg_idx:
        row_idx=reg_idx[region]
        for motif in motif_idx:
            col_idx=motif_idx[motif]
            if region in MaxScoredict and motif in MaxScoredict[region]:
                score=MaxScoredict[region][motif]
                featureMatrix_POS[row_idx,col_idx] = score
    return featureMatrix_POS



def convert_fasta_to_string( id_sequence_dict ):
    fasta_sequence_string = ''
    for key in id_sequence_dict:
        fasta_sequence_string += ">" + key + "\n"
        fasta_sequence_string += id_sequence_dict[key] + "\n"
    return fasta_sequence_string


def save_results_old(regionID_list, prob_posclass_ref, prob_posclass_mut, start, indel_descr):
    for i in range(0, len(regionID_list)):
        chr_start_end = regionID_list[i]
        # chr_start_end = chr_start_end.replace(":","\t")
        # chr_start_end = chr_start_end.replace("-","\t")
        FILE_TO_SAVE_HANDLE.write(chr_start_end + "\t" + str(start) + "\t" + indel_descr + "\t" + str(prob_posclass_ref[i]) + "\t" + str(prob_posclass_mut[i])  + "\n")
    FILE_TO_SAVE_HANDLE.flush()


def save_results(regionID_list, prob_posclass_ref, prob_posclass_mut, start, indel_descr, indel_additiona_information):
    ref_max_score_dict = dict()
    mut_max_score_dict = dict()
    for i in range(0, len(regionID_list)):
        insertion_desciption = regionID_list[i].split(":")[0] + "\t" + str(start) + "\t" + indel_descr
        if insertion_desciption not in ref_max_score_dict:
            ref_max_score_dict[insertion_desciption] = float(prob_posclass_ref[i])
        else:
            if ref_max_score_dict[insertion_desciption] < float(prob_posclass_ref[i]):
                ref_max_score_dict[insertion_desciption] = float(prob_posclass_ref[i])

        if insertion_desciption not in mut_max_score_dict:
            mut_max_score_dict[insertion_desciption] = float(prob_posclass_mut[i])
        else:
            if mut_max_score_dict[insertion_desciption] < float(prob_posclass_mut[i]):
                mut_max_score_dict[insertion_desciption] = float(prob_posclass_mut[i])

    for region in ref_max_score_dict:
        if region in mut_max_score_dict:
            string_to_write = region + "\t" + indel_additiona_information + "\t" + str(ref_max_score_dict[region]) + "\t" + str(mut_max_score_dict[region]) + "\n"
            FILE_TO_SAVE_HANDLE.write(string_to_write)
            FILE_TO_SAVE_HANDLE.flush()
    FILE_TO_SAVE_HANDLE.flush()



def LoadRFClassifier(path_to_model):
    RFClassifier = joblib.load(path_to_model)
    return RFClassifier


def check_if_to_continue_scoring(path_to_scored_results):
    """
    chr6	161092882	G:GA	0.05353200883	0.05353200883
    """
    mutID_dict=dict()
    with open(path_to_scored_results,'r') as my_file:
        for line in my_file:
            line = line.split()
            chrom = line[0]
            indel_data = line[2]
            mut_start = line[1]
            indel_information = line[3]
            mut_ID = chrom + "_" + mut_start + "_" + indel_data + "_" + indel_information
            mutID_dict[mut_ID] = mut_ID
    return mutID_dict


def score_ref_mut_region(regionID_fasta_dict_ref, regionID_fasta_dict_mut,start, indel_descr, indel_additiona_information):
    ###---Run Cbust---###
    region_order_ref = regionID_fasta_dict_ref.keys()
    fasta_string_ref = convert_fasta_to_string(regionID_fasta_dict_ref)
    cbust_allmotifs_merged_results_ref = run_cbust_scan_pipe(FeatureOrder_list, fasta_string_ref)
    featureMatrix_region_ref = CRMmaxScore_Homo_Cluster(cbust_allmotifs_merged_results_ref, FeatureOrder_list, region_order_ref)
    ###----Make classification----###
    probas_ref = Classifier.predict_proba(featureMatrix_region_ref)
    probas_ref =  probas_ref[:,1]
    ###---Run Cbust---###
    region_order_mut = regionID_fasta_dict_mut.keys()
    fasta_string_mut = convert_fasta_to_string(regionID_fasta_dict_mut)
    cbust_allmotifs_merged_results_mut = run_cbust_scan_pipe(FeatureOrder_list, fasta_string_mut)
    featureMatrix_region_mut = CRMmaxScore_Homo_Cluster(cbust_allmotifs_merged_results_mut, FeatureOrder_list, region_order_mut)
    ###----Make classification----###
    probas_mut = Classifier.predict_proba(featureMatrix_region_mut)
    probas_mut =  probas_mut[:,1]
    ###---chek that region order for mutant and reference is the same---###
    IS_ORDER_OK = True
    for i in range(0, len(region_order_mut)):
        if region_order_mut[i]!=region_order_ref[i]:
            IS_ORDER_OK = False
            print "ERROR!!! Order of regions for ref and mut sequence is not the same"
    if IS_ORDER_OK:
        save_results( region_order_mut, probas_ref ,probas_mut, start, indel_descr, indel_additiona_information)

"""
Main block
"""
if __name__ == '__main__':
    fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1)
    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage = usage, version = "%prog v1.0", formatter = fmt)
    parser = optparse.OptionParser()
    parser.add_option("-o", "--feature_order", action = "store", type = "string", dest = "feat_order_path", help = 'Name of the TF to process')
    parser.add_option("-b", "--bed_file", action = "store", type = "string", dest = "bed_file", help = 'Path to .bed file with mutations')
    parser.add_option("-s", "--save_path", action = "store", type = "string", dest = "path_to_save_results", help = 'Path to save results')
    parser.add_option("-f", "--fasta_path", action = "store", type = "string", dest = "path_to_fasta", help = 'Path to fasta file')
    parser.add_option("-m", "--pcl_path", action = "store", type = "string", dest = "pcl_path", help = 'Path to .pkl file with model')
    parser.add_option("-w", "--window_size", action = "store", type = "int", dest = "window_size", help = 'Size of the sliding window')
    parser.add_option("-p", "--pwms_path", action = "store", type = "string", dest = "pwms_path", help = 'Path to the folder with PWMs')
    parser.add_option("-c", "--cbust_path", action = "store", type = "string", dest = "cbust_path", help = 'Path to Cluster-buster tool')
    (options, args) = parser.parse_args()
    # Check if we have an expression matrix filea FASTA or twobit file is given as input.
    if ( (options.feat_order_path is None) or (options.bed_file is None) or (options.path_to_save_results is None)  \
                 or (options.path_to_fasta is None) or (options.pcl_path is None) or (options.window_size is None) \
                 or (options.pwms_path is None) or (options.cbust_path is None)):
        parser.print_help()
        print >> sys.stderr, '\nERROR: minimum required options not satisfied:\n'
        sys.exit(1)

    PATH_TO_FEATURE_ORDER = options.feat_order_path
    PATH_TO_FASTA_FILE = options.path_to_fasta
    PATH_TO_INDELS = options.bed_file
    PATH_TO_SAVE_RESULTS = options.path_to_save_results
    PATH_TO_PKL_FILE = options.pcl_path
    SLIDING_WINDOW_SIZE = options.window_size
    PATH_TO_CBUST = options.cbust_path
    PATH_TO_FOLDER_WITH_SINGLETONS = options.pwms_path

    chrID_fasta_dict = bx.seq.twobit.TwoBitFile( open(PATH_TO_FASTA_FILE, 'rb') )
    FeatureOrder_list = readFetureOrder( PATH_TO_FEATURE_ORDER )
    Classifier = LoadRFClassifier(PATH_TO_PKL_FILE)
    coord_muttype_dict = read_file_with_indels(PATH_TO_INDELS)
    FILE_TO_SAVE_HANDLE = open(PATH_TO_SAVE_RESULTS,'a')

    scored_already_insertions = check_if_to_continue_scoring(PATH_TO_SAVE_RESULTS)

    for indel in coord_muttype_dict:
        if indel not in scored_already_insertions:
            indel_specific_info = "_".join(indel.split("_")[3:])
            indel_data =  indel.split( "_" )
            chrom = indel_data[0]
            start = int(indel_data[1])
            indel_descr = coord_muttype_dict[indel]
            regionID_fasta_dict_ref, regionID_mut_fasta_dict = k_sliding_windows2fasta_deletion( chrom, start, indel_descr, SLIDING_WINDOW_SIZE )
            score_ref_mut_region( regionID_fasta_dict_ref, regionID_mut_fasta_dict,start, indel_descr, indel_specific_info)
    FILE_TO_SAVE_HANDLE.close()