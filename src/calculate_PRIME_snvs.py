import sys
import subprocess
import bx.seq.twobit
from array import *
import os
from collections import defaultdict
import numpy as np
import random
from sklearn.externals import joblib
from collections import OrderedDict
import optparse


###----functions----###
def readFetureOrder(path_to_file):
    FeatureOrder_list=[]
    with open(path_to_file) as my_file:
        for line in my_file:
            motif_name = line.split()[0]
            FeatureOrder_list.append(motif_name)
    return FeatureOrder_list


def read_bed_file(path_to_bed_file):
    bed_file_regionsID_coord_dict=dict()
    with open(path_to_bed_file) as my_file:
        for line in my_file:
            line=line.split()
            chr=line[0]
            start=int(line[1])
            end=int(line[2])
            regionID=line[3]
            if regionID in bed_file_regionsID_coord_dict:
                print regionID, " is not unique"
                continue
            bed_file_regionsID_coord_dict[regionID] = [chr,start,end]
    return bed_file_regionsID_coord_dict



def select_fasta_for_region(chrID_fasta_dict, region_bed_format_list):
    chr_name = region_bed_format_list[0]
    start = int(region_bed_format_list[1])
    end = int(region_bed_format_list[2])
    fasta_sequence=chrID_fasta_dict[chr_name].get(start,end).upper()
    fasta_sequence = array('c',fasta_sequence)
    return fasta_sequence



###-----read filr with type of mutation and position of the mutation-----###
def read_mut_file(path_to_file):
    coord_muttype_dict = OrderedDict()
    with open(path_to_file) as my_file:
        for line in my_file:
            line = line.split()
            chrom = line[0]
            start = line[1]
            end = line[2]
            mut_type = line[3] + ":" + line[4]

            mut_additional_info = line[5]
            save_to_dict = mut_type + ":" + mut_additional_info
            mut_type = line[3].split(",")[0] + ":" + line[4].split(",")[0]
            region_mutID = chrom + "_" + start + "_" + end + "_" + mut_type + "_" + mut_additional_info
            if region_mutID in  coord_muttype_dict:
                print region_mutID, " not unique"
            coord_muttype_dict[region_mutID] = mut_type
    return coord_muttype_dict


def read_mut_file_recurrent_only(path_to_file):
    coord_muttype_dict = dict()
    coord_muttype_dict_recurent = dict()
    with open(path_to_file) as my_file:
        for line in my_file:
            line = line.split()
            chrom = line[0]
            start = line[1]
            end = line[2]
            mut_type = line[3] + ":" + line[4]
            region_mutID=chrom + "_" + start + "_" + end
            if region_mutID in  coord_muttype_dict:
                print region_mutID, " not unique"
                coord_muttype_dict_recurent[region_mutID] = mut_type
            coord_muttype_dict[region_mutID] = mut_type
    return coord_muttype_dict_recurent



###---read fasta file----###
def read_fastafile(filename):
    """
    Fasta file containing whole genome must
    have chomosomes ID as fastaID names
    """
    id = ''
    ids = []
    seqs = []
    seq = []
    with open(filename, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                regionID=line[1:].rstrip('\n')
                regionID = regionID.split(":")[0]
                ids.append(regionID)
                if seq != []: seqs.append("".join(seq))
                seq = []
            else:
                seq.append(line.rstrip('\n').upper())
        if seq != []:
            seqs.append("".join(seq))
        seqID_fasta = dict(zip(ids,seqs))
        ###----convert strings to python array----###
        for key in seqID_fasta:
            seqID_fasta[key] = array('c',seqID_fasta[key])
    return seqID_fasta





def saveFastaFile(path_to_save_fasta_file, seqID_fasta_dict):
    file_to_save=open(path_to_save_fasta_file,"w")
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
            chrom = str(line[3])   ###chr:start-end
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



def run_cbust_scan_pipe( singletone_list, modified_fasta_seq ):
    ###---Variable to save results of scanning---###
    initScanningData=''
    for motif in singletone_list:
        stdoutdata, stderrdata = run_cmd([PATH_TO_CBUST, '-c', '0', '-m', '0', '-f', '3', PATH_TO_FOLDER_WITH_SINGLETONS + motif], modified_fasta_seq)
        initScanningData = merge2oneFileSingletonScanWG( stdoutdata, initScanningData )
    return initScanningData


def convert_fasta_to_string( id_sequence_dict ):
    fasta_sequence_string = ''
    for key in id_sequence_dict:
        fasta_sequence_string += ">" + key + "\n"
        fasta_sequence_string += id_sequence_dict[key] + "\n"
    return fasta_sequence_string



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



def generate_tmp_file_name(path_to_candidate_name):
    while os.path.exists(path_to_candidate_name):
        rnd_int = random.randint(10,1000)
        path_to_candidate_name = path_to_candidate_name[:-4] + str(rnd_int) + '.tmp'
        if not os.path.exists(path_to_candidate_name):
            return path_to_candidate_name
    else:
        return path_to_candidate_name


def LoadRFClassifier(path_to_model):
    RFClassifier = joblib.load(path_to_model)
    return RFClassifier


def sliding_windows2fasta(mut_coord, mut_type, sw_size, regionID_fasta_dict):
    """
    regionID_fasta_dict - key: region ID(chr:start-end; value: fasta sequence)
    """
    ###----define start of the sliding window----###
    mut_start_chr = mut_coord[0]
    mut_start = int(mut_coord[1])

    ###----open file to save results-----###


    ##---iterate n_tmes=sw size----###
    for i in range(0, sw_size+1):
        start_sw = mut_start - sw_size + i
        end_sw = start_sw + sw_size+1
        chr_start_end_list= [mut_start_chr, start_sw, end_sw]
        ###----select fasta sequence for the region------###
        fasta_sequence_ref = select_fasta_for_region(chrID_fasta_dict, chr_start_end_list)
        ###---find mutation to introduce----###
        ref_nuc = mut_type[0]
        mut_nuc = mut_type[1]
        mut_start_local_coord = mut_start - start_sw

        ###---Make predictions for reference sequence---###
        regionID = mut_start_chr + ":" + str(start_sw) + "-" + str(end_sw) + ":" + ref_nuc + ":" + mut_nuc

        fasta_sequence = select_fasta_for_region(chrID_fasta_dict, chr_start_end_list)
        ###---Mutate fasta---##
        if fasta_sequence[mut_start_local_coord] == ref_nuc:
            fasta_sequence = fasta_sequence.tostring()
            if regionID not in regionID_fasta_dict:
                regionID_fasta_dict[regionID] = fasta_sequence
        else:
            # print "Wrong reference nucleotide reference (sliding_windows2fasta function)", mut_start_chr, mut_start ,ref_nuc, mut_nuc
            break
        ###---Save fasta file corresponding to sequence---###


def k_sliding_windows2fasta( mut_coord, mut_type, sw_size, regionID_fasta_dict ):
    """
    regionID_fasta_dict - key: region ID(chr:start-end; value: fasta sequence)
    """
    ###----define start of the sliding window----###
    mut_start_chr = mut_coord[0]
    mut_start = int(mut_coord[1])

    ###----open file to save results-----###

    shift = int(sw_size/10)
    ##---iterate n_tmes=sw size----###
    for i in range(0, sw_size+1,shift ):
        start_sw = mut_start - sw_size + i
        end_sw = start_sw + sw_size + 1
        chr_start_end_list= [mut_start_chr,start_sw,end_sw]
        ###----select fasta sequence for the region------###
        ###---find mutation to introduce----###
        ref_nuc = mut_type[0]
        mut_nuc = mut_type[1]
        mut_start_local_coord = mut_start - start_sw

        ###---Make predictions for reference sequence---###
        regionID = mut_start_chr + ":" + str(start_sw) + "-" + str(end_sw) + ":" + ref_nuc + ":" + mut_nuc

        fasta_sequence = select_fasta_for_region(chrID_fasta_dict, chr_start_end_list)
        ###---Mutate fasta---##
        if fasta_sequence[mut_start_local_coord] == ref_nuc:
            fasta_sequence = fasta_sequence.tostring()
            if regionID not in regionID_fasta_dict:
                regionID_fasta_dict[regionID] = fasta_sequence
        else:
            print "Reference nucleotide differs from hg19 version", mut_start_chr, mut_start ,ref_nuc, mut_nuc
            fasta_sequence = fasta_sequence.tostring()
            if regionID not in regionID_fasta_dict:
                regionID_fasta_dict[regionID] = fasta_sequence
        ###---Save fasta file corresponding to sequence---###

def k_sliding_windows2fasta_mut( mut_coord, mut_type, sw_size, regionID_fasta_dict ):
    """
    regionID_fasta_dict - key: region ID(chr:start-end; value: fasta sequence)
    """
    ###----define start of the sliding window----###
    mut_start_chr = mut_coord[0]
    mut_start = int(mut_coord[1])

    ###----open file to save results-----###

    shift = int(sw_size/10)
    num_windows = int(sw_size/shift)
    ##---iterate n_tmes=sw size----###
    for i in range(0, sw_size+1,shift ):
        start_sw = mut_start - sw_size + i
        end_sw = start_sw + sw_size + 1
        chr_start_end_list= [mut_start_chr, start_sw, end_sw]
        ###----select fasta sequence for the region------###
        ###---find mutation to introduce----###
        ref_nuc = mut_type[0]
        mut_nuc = mut_type[1]
        mut_start_local_coord = mut_start - start_sw

        ###---Make predictions for reference sequence---###

        regionID = mut_start_chr + ":" + str(start_sw) + "-" + str(end_sw) + ":" + ref_nuc + ":" + mut_nuc
        fasta_sequence = select_fasta_for_region(chrID_fasta_dict, chr_start_end_list)
        ###---Mutate fasta---##
        if fasta_sequence[mut_start_local_coord] == ref_nuc:
            fasta_sequence[mut_start_local_coord] = mut_nuc
            fasta_sequence = fasta_sequence.tostring()
            if regionID not in regionID_fasta_dict:
                regionID_fasta_dict[regionID] = fasta_sequence
        else:
            print "Reference nucleotide in the input file differs .2bit file", mut_start_chr, mut_start ,ref_nuc, mut_nuc

            fasta_sequence[mut_start_local_coord] = mut_nuc
            fasta_sequence = fasta_sequence.tostring()
            if regionID not in regionID_fasta_dict:
                regionID_fasta_dict[regionID] = fasta_sequence

def check_if_to_continue_scoring(path_to_scored_results):
    """
    chr1	224548373	224548374	C	T	CNIH4:TCGA-A7-A426	0.0694260485651	0.0627404604226
    chr1	222895653	222895654	C	A	C1orf58:TCGA-A7-A426	0.0524282560706	0.0524282560706
    chr1	220213428	220213429	G	T	EPRS:TCGA-A7-A426	0.0419426048565	0.0440397350993
    chr1	207143259	207143260	G	T	FCAMR:TCGA-A7-A4266	0.297114474929	0.30307473983
    """
    mutID_dict=dict()
    with open(path_to_scored_results) as my_file:
        for line in my_file:
            line = line.split()
            chrom = line[0]
            mut_type_ref = line[3]
            mut_type_mut = line[4]
            mut_start = line[1]
            mut_end = line[2]
            regionID = line[5]
            mut_ID = chrom + "_" + mut_start + "_" + mut_end + "_" + mut_type_ref + ":" + mut_type_mut + "_" + regionID
            mutID_dict[mut_ID] = mut_ID
    return mutID_dict



###----Save results of the predictions-----###
def save_results(regionID_list, prob_posclass_ref, prob_posclass_mut, mut_coord):
    regionID_max_score_ref = dict()
    regionID_max_score_mut = dict()
    for i in range(0,len(regionID_list)):
        mut_region_ID =mut_coord[0] + "\t" + str(mut_coord[1]) + "\t" + str(mut_coord[2]) + "\t" \
                              + str(mut_coord[3].split(":")[0]) + "\t" \
                              + str(mut_coord[3].split(":")[1]) + "\t" + str("_".join(mut_coord[4:]))
        score_ref = prob_posclass_ref[i]
        score_mut = prob_posclass_mut[i]
        if mut_region_ID not in regionID_max_score_ref:
            regionID_max_score_ref[mut_region_ID] = score_ref
        else:
            if regionID_max_score_ref[mut_region_ID] < score_ref:
                regionID_max_score_ref[mut_region_ID] = score_ref

        if mut_region_ID not in regionID_max_score_mut:
            regionID_max_score_mut[mut_region_ID] = score_mut
        else:
            if regionID_max_score_mut[mut_region_ID] < score_mut:
                regionID_max_score_mut[mut_region_ID] = score_mut

    for region in regionID_max_score_ref:
        if region in regionID_max_score_mut:
            string_to_write = region + "\t" + str(regionID_max_score_ref[region]) + "\t" + str(regionID_max_score_mut[region]) + "\n"
            FILE_TO_SAVE_HANDLE.write(string_to_write)
            FILE_TO_SAVE_HANDLE.flush()




def k_sw_mut_scoring(region_mut):
    regionID_fasta_dict_ref = OrderedDict()
    regionID_fasta_dict_mut = OrderedDict()
    ###---For each mutation in run scoring with sliding window---###
    mut_coord = region_mut.split("_")
    mut_type = mut_coordID_dict[region_mut].split(":")
    sw_size = SLIDING_WINDOW_SIZE
    ###----Select fasta regions around mutation----###
    k_sliding_windows2fasta( mut_coord, mut_type, sw_size, regionID_fasta_dict_ref )
    ###---Run Cbust---###
    region_order_ref = regionID_fasta_dict_ref.keys()
    fasta_string_ref = convert_fasta_to_string(regionID_fasta_dict_ref)
    cbust_allmotifs_merged_results_ref = run_cbust_scan_pipe(FeatureOrder_list, fasta_string_ref)
    featureMatrix_region_ref = CRMmaxScore_Homo_Cluster(cbust_allmotifs_merged_results_ref, FeatureOrder_list, region_order_ref)
    ###----Make classification----###
    probas_ref = Classifier.predict_proba(featureMatrix_region_ref)
    probas_ref =  probas_ref[:,1]
    ###----Select fasta regions around mutation for each mutation----###
    k_sliding_windows2fasta_mut( mut_coord, mut_type, sw_size, regionID_fasta_dict_mut)
    ###---Run Cbust---###
    region_order_mut = regionID_fasta_dict_mut.keys()
    fasta_string_mut = convert_fasta_to_string(regionID_fasta_dict_mut)
    cbust_allmotifs_merged_results_mut = run_cbust_scan_pipe(FeatureOrder_list, fasta_string_mut)
    featureMatrix_region_mut = CRMmaxScore_Homo_Cluster(cbust_allmotifs_merged_results_mut, FeatureOrder_list, region_order_mut)
    ###----Make classification----###
    probas_mut = Classifier.predict_proba(featureMatrix_region_mut)
    probas_mut =  probas_mut[:,1]
    ###---chek that region order for mutant and reference is the same---###
    IS_ORDER_OK=True
    for i in range(0,len(region_order_mut)):
        if region_order_mut[i]!=region_order_ref[i]:
            IS_ORDER_OK=False
            print "ERROR!!! Order of regions for ref and mut sequence is not the same"
    # save results if order is OK
    if IS_ORDER_OK:
        save_results(region_order_mut, probas_ref, probas_mut, mut_coord)



"""
Main block
"""


if __name__ == '__main__':

    fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1)
    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage = usage, version = "%prog v1.0", formatter = fmt)
    parser.add_option("-o", "--feature_order", action = "store", type = "string", dest = "feat_order_path", help = 'File with feature names')
    parser.add_option("-b", "--bed_file", action = "store", type = "string", dest = "bed_file", help = 'file with mutations (SNVs)')
    parser.add_option("-s", "--save_path"   , action = "store", type = "string", dest = "path_to_save_results", help = 'Path to save results')
    parser.add_option("-f", "--fasta_path", action = "store", type = "string", dest = "path_to_fasta", help = 'Path to .2bit file')
    parser.add_option("-m", "--pcl_path", action = "store", type = "string", dest = "pcl_path", help = 'Path to .pcl file with model')
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
    PATH_TO_BED_FILE_MUST = options.bed_file
    PATH_TO_SAVE_RESULTS = options.path_to_save_results
    PATH_TO_PKL_FILE = options.pcl_path
    SLIDING_WINDOW_SIZE = options.window_size
    PATH_TO_CBUST = options.cbust_path
    PATH_TO_FOLDER_WITH_SINGLETONS = options.pwms_path
    ###---Control parameters-----###
    ###---Read fasta---###
    print 'Read fasta'
    chrID_fasta_dict = bx.seq.twobit.TwoBitFile( open(PATH_TO_FASTA_FILE, 'rb') )
    ###----Read file with mutations----###
    mut_coordID_dict = read_mut_file( PATH_TO_BED_FILE_MUST )
    mut_coordID_dict_keys = mut_coordID_dict.keys()
    ###---Read feature order---###
    FeatureOrder_list = readFetureOrder( PATH_TO_FEATURE_ORDER )
    ###----Open file descriptor----###
    num_regions_to_score = len(mut_coordID_dict.keys())
    print num_regions_to_score, " mutations to score"
    ###---run sliding window scoring for each mutation----###
    ###---Load RF model(.pkl file)---###
    Classifier = LoadRFClassifier(PATH_TO_PKL_FILE)
    ###---Read what mutations are already scored----##
    FILE_TO_SAVE_HANDLE = open(PATH_TO_SAVE_RESULTS,'a')
    scored_already_mutations_dict = check_if_to_continue_scoring(PATH_TO_SAVE_RESULTS)
    for key_mut in mut_coordID_dict_keys:
        chr_name_mut = key_mut.split("_")[0]
        if key_mut not in scored_already_mutations_dict:
            if chr_name_mut in chrID_fasta_dict:
                k_sw_mut_scoring(key_mut)

    FILE_TO_SAVE_HANDLE.close()