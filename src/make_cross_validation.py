
"""
This script makes cross validation procedure
"""

import optparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import precision_recall_curve
import pandas as pd
import numpy as np



###----Classifier parameters----###
n_trees = 151
max_Features = "sqrt"
Min_samples_split = 2
n_jobs = 1
Min_samples_leaf = 2
max_depth = None
random_state = None
n_folds=10
###################################





ROC_Cbust = []
PR_Cbust = []

ROC_MEAN_RF = []
PR_MEAN_RF = []


parser = optparse.OptionParser()
parser.add_option("-p", "--pos_path", action="store", type="string", dest="pos_path", help='Name of the TF to process')
parser.add_option("-n", "--neg_path", action = "store", type = "string", dest = "neg_path", help = 'Path to ctx-rcc results')
parser.add_option("-s", "--save_path", action = "store", type = "string", dest = "save_path", help = 'Path to ctx-rcc results')
parser.add_option("-c", "--save_cross_validation", action = "store", type = "string", dest = "save_cross_validation", help = 'Path to ctx-rcc results')
(options, args) = parser.parse_args()



PATH_TO_POS_FILE = options.pos_path
PATH_TO_NEG_FILE = options.neg_path
PATH_TO_SAVE = options.save_path
PATH_TO_SAVE_CROSS_VALIDATION_SCOES = options.save_cross_validation

file_to_save = open(PATH_TO_SAVE, 'w')
###---BEGIN


###----Read data using pandas FT style-----###

data_table_seq_pos_pandas = pd.read_table(PATH_TO_POS_FILE,sep="\t", header=0, index_col=0)
data_table_seq_neg_pandas = pd.read_table(PATH_TO_NEG_FILE,sep="\t", header=0, index_col=0)
###---convert data frame to numpy array---###

###---convert data to
data_table_pos = pd.np.array(data_table_seq_pos_pandas)
data_table_neg = pd.np.array(data_table_seq_neg_pandas)
class_labels_pos = [1]*data_table_pos.shape[0]
class_labels_neg = [0]*data_table_neg.shape[0]



DataTable_all = np.append(data_table_pos, data_table_neg, axis = 0)
class_labels_all = class_labels_pos + class_labels_neg




###---------------------------------###
classifier = RandomForestClassifier(max_depth=max_depth, n_estimators = n_trees, max_features = max_Features, random_state = random_state, min_samples_split = Min_samples_split, min_samples_leaf=Min_samples_leaf)

###--classifier.fit(DataTable_all, np.array(class_labels_all))
X = DataTable_all
y = np.array(class_labels_all)


###------------------------###


path_to_save_region_rfsocres=PATH_TO_SAVE[:-4] + 'rf_scores.txt'

cv = StratifiedKFold(np.array(y), n_folds = n_folds)
mean_tpr = 0.0
mean_fpr = np.linspace(0, 1, 200)
all_tpr = []
###----ROC and PR curve with cross-validation-----###
auc_roc_list = []
aupr_mean = []
file_to_save_scoring=open(PATH_TO_SAVE_CROSS_VALIDATION_SCOES, 'w')
for i, (train, test) in enumerate(cv):
    probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
    ###----Save cross_validation scores----###
    prob_idx=0
    for idx in test:
        class_label_tmp = y[idx]
        score_tmp = probas_[prob_idx,1]
        prob_idx = prob_idx + 1
        file_to_save_scoring.write(str(class_label_tmp) + "\t" + str(score_tmp) + "\n")


    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
    mean_tpr += interp(mean_fpr, fpr, tpr)
    mean_tpr[0] = 0.0
    roc_auc = auc(fpr, tpr)
    auc_roc_list.append(roc_auc)
    ###---Compute PR for
    precision, recall, thresholds = precision_recall_curve(y[test], probas_[:, 1])
    area = auc(recall, precision)
    aupr_mean.append(area)
    print "fold", i, "ROC=", roc_auc, "PR=", area
PR_MEAN_RF.append(mean(aupr_mean))
ROC_MEAN_RF.append(mean(auc_roc_list))

file_to_save.write(str(mean(aupr_mean)) + "\t" + str(mean(auc_roc_list)) + "\n")
###----RF end-----###
file_to_save.close()
