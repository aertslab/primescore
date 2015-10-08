###-----import-------###
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import os
import optparse
import sys
from sklearn.externals import joblib
import pandas as pd

# Build command line options.
fmt = optparse.IndentedHelpFormatter(indent_increment = 2, max_help_position = 9, width = 79, short_first = 1)

usage = "Usage: %prog [options]"

parser = optparse.OptionParser(usage=usage, version="%prog v1.0", formatter=fmt)

###----add command line arguments----###
#parser = optparse.OptionParser()
parser.add_option("-p", "--path_to_pos_sample",action="store", type="string",dest="path_to_pos_sample", help='Path to FT for positives')
parser.add_option("-n", "--path_to_neg_sample",action="store", type="string",dest="path_to_neg_sample", help='Path to FT for negatives')
parser.add_option("-s", "--path_save_model",action="store", type="string",dest="path_save_model", help='Path to save model with trained classifier')
parser.add_option("-t", "--tree_number",action="store", type="int",dest="n_tries",default = 151, help='Path to save model with trained classifier')
parser.add_option("-f", "--max_Features",action="store", type="string",dest="max_Features",default = "sqrt", help='Path to save model with trained classifier')
(options, args) = parser.parse_args()

# Check if we have an expression matrix filea FASTA or twobit file is given as input.
if ( (options.path_to_pos_sample is None) or (options.path_to_neg_sample is None) or (options.path_save_model is None) ):
    parser.print_help()
    print >> sys.stderr, '\nERROR: minimum required options not satisfied:\n         - Path to FT for positives\n         - Path to FT for negatives\n         - Path to save model with trained classifier\n'
    sys.exit(1)



###----global variables-----###
PATH_TO_POS_FILE = options.path_to_pos_sample
PATH_TO_NEG_FILE = options.path_to_neg_sample



###----Classifier parameters----###
n_trees = options.n_tries
max_Features = options.max_Features

Min_samples_split = 2
n_jobs = 1
Min_samples_leaf=2
max_depth = None
random_state=None
###################################

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        print "Created folder", f


class TrainClassifier():

    RFClassifier=None
    def __init__(self,path_to_train_data_pos_= PATH_TO_POS_FILE,
                 path_to_train_data_neg_=PATH_TO_NEG_FILE,
                 n_trees=n_trees,
                 max_Features=max_Features,
                 Min_samples_split=Min_samples_split,
                 n_jobs=n_jobs):
        self.path_to_train_data_pos = path_to_train_data_pos_
        self.path_to_train_data_neg = path_to_train_data_neg_
        self.n_trees = n_trees
        self.max_Features=max_Features
        self.Min_samples_split=Min_samples_split
        self.min_samples_leaf=Min_samples_leaf
        self.n_jobs=n_jobs

    def TrainRFClassifier(self):
        data_pos = pd.read_table(self.path_to_train_data_pos,sep="\t", header=0, index_col=0)
        data_pos = pd.np.array(data_pos)
        data_neg = pd.read_table(self.path_to_train_data_neg,sep="\t", header=0, index_col=0)
        data_neg = pd.np.array(data_neg)


        data_label = np.array([1]*data_pos.shape[0] + [0]*data_neg.shape[0])
        data_all = np.append(data_pos, data_neg, axis = 0)

        self.RFClassifier=RandomForestClassifier(max_depth=max_depth,
                                                 n_estimators = self.n_trees,
                                                 max_features = self.max_Features,
                                                 random_state = random_state,
                                                   min_samples_split = self.Min_samples_split,
                                                 min_samples_leaf=Min_samples_leaf,
                                                 n_jobs = self.n_jobs)


        self.RFClassifier.fit(data_all,data_label)
        ###---Save trained model to the file----###
        path_to_save_trained_model = options.path_save_model
        ###---open file in binary mode to pickle the classifier----##
        joblib.dump(self.RFClassifier,path_to_save_trained_model,compress=1)
        print "Classifier was trained and saved to the file..", path_to_save_trained_model



"""
Main block
"""

###---Make object----###
clf=TrainClassifier()
###----Train classifier----###
clf.TrainRFClassifier()