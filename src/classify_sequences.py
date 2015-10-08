__author__ = 'dsvet'

"""
This scripts makes WG prediction using FT
In pandas format: with row names (sliding window regionID ) and colnames: Motif names ID
"""

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt
import os
import optparse
from sklearn.externals import joblib
import pandas as pd
import sys




###----add command line arguments----###
fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1)
usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage = usage, version = "%prog v1.0", formatter = fmt)

parser.add_option("-d", "--data_path", action="store", type="string", dest="data_path", help='Path to feature-table file')
parser.add_option("-m", "--model_file", action="store", type="string", dest="path_to_model", help='Path to .pcl file with trained classifier')
parser.add_option("-s", "--save_results", action="store", type="string", dest="path_to_save_results", help='Path to save results')
(options, args) = parser.parse_args()

# Check if we have an expression matrix filea FASTA or twobit file is given as input.
if ( (options.data_path is None) or (options.path_to_model is None) or (options.path_to_save_results is None) ):
    parser.print_help()
    print >> sys.stderr, '\nERROR: minimum required options not satisfied:\n'
    sys.exit(1)


###----global variables-----###

PATH_TO_PKL_FILE = options.path_to_model
PATH_TO_SCORED_CHROMOSOMES = options.data_path
PATH_TO_SAVE_PREDICTIONS_per_CHR = options.path_to_save_results

###---Load RF classifier---###
def LoadRFClassifier(path_to_pkl_file):
    path_to_learned_clf = path_to_pkl_file
    RFClassifier = joblib.load(path_to_learned_clf)
    return RFClassifier


def read_region_ids(path_to_regionID_file):
    regionID_list=[]
    with open(path_to_regionID_file,'r') as regionid_file:
        for line in regionid_file:
            line = line.split()
            regionID_list.append(line[3])
    return regionID_list

def save_results(regionID_list, prob_posclass, path_to_save):
    file_to_save=open(path_to_save,'w')
    for i in range(0,len(regionID_list)):
        chr_start_end = regionID_list[i]
        file_to_save.write(chr_start_end + "\t" + str(prob_posclass[i])+"\n")



"""
Main block
"""


###---Load trained classifier---###
RFClassifier = LoadRFClassifier(PATH_TO_PKL_FILE)

###----Read FT----###
data_to_predict_df = pd.read_table(PATH_TO_SCORED_CHROMOSOMES,sep="\t", header=0, index_col=0)
regionID_list = data_to_predict_df.index
###---conver DF to numpy array---###
data_to_predict = pd.np.array(data_to_predict_df)
probas_ = RFClassifier.predict_proba(data_to_predict)
probas_ =  probas_[:,1]
print "Prediction was finished"
print "Scored number of regions...", len(regionID_list)
print "Saving results"
save_results(regionID_list, probas_, PATH_TO_SAVE_PREDICTIONS_per_CHR)


























