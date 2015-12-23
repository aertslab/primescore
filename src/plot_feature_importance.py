import numpy as np
import matplotlib.pyplot as plt
from scipy import interp
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from pylab import *
import optparse


PATH_TO_DATA_POS = '/home/dsvet/mydata/45TFs_new/TRAINING_SET/ft_pos_neg_1pertf_all45_tfs_with_fmf/STAT5A.pos.FT.ftf'
PATH_TO_DATA_NEG = '/home/dsvet/mydata/45TFs_new/TRAINING_SET/ft_pos_neg_1pertf_all45_tfs_with_fmf/STAT5A.neg.FT.ftf'


###----Calculate feature importance----###

def calculate_feature_importance(data_for_training_numpy, list_all_features, class_labels_bin,filename_to_save_fig):
    n_trees=151
    ###---Make cross validation and plot curve using all features---###
    y = np.array(class_labels_bin)
    X = data_for_training_numpy
    feature_importance_per_trie=np.zeros((n_trees, len(list_all_features)))
    forest = RandomForestClassifier(n_estimators = n_trees)
    forest.fit(X, y)
    importances = forest.feature_importances_
    indices = np.argsort(importances)[::-1]
    row_idx=0
    for tree in forest.estimators_:
        feature_importance_per_trie[row_idx,:] = tree.feature_importances_
        row_idx=row_idx+1
    feature_importance_std = np.std(feature_importance_per_trie, axis=0)
    importance_std_mean = feature_importance_std/np.sqrt(len(feature_importance_std))

    ###---Make cross validation and plot curve using all features---###

    ###---Make cross validation and plot curve using all features---###
    color_PWM = '#80E800'
    num_features = len(list_all_features)
    # plt.figure(figsize=(6, 8), dpi=300)
    plt.plot(np.arange(0, num_features), importances, marker='o',linestyle='--',color='#8667D7')
    plt.fill_between(np.arange(0,num_features), importances[:num_features], importances[:num_features] + importance_std_mean[:num_features], color = color_PWM, alpha=0.5)
    plt.fill_between(np.arange(0,num_features), importances[:num_features], importances[:num_features] - importance_std_mean[:num_features], color = color_PWM, alpha=0.5)
    plt.xticks(np.arange(0, num_features), list_all_features, rotation='vertical',fontsize = 10)
    plt.tight_layout()
    plt.savefig(filename_to_save_fig, dpi=80)


    ###---Print to stdout feature importance value---###

    # plt.show()


if __name__ == '__main__':

    fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1)
    usage = "Usage: %prog [options]"
    parser = optparse.OptionParser(usage = usage, version = "%prog v1.0", formatter = fmt)
    parser.add_option("-p", "--pos_data", action = "store", type = "string", dest = "pos_data_path", help = 'File with feature names')
    parser.add_option("-n", "--neg_data", action = "store", type = "string", dest = "neg_data_path", help = 'file with mutations (SNVs)')
    parser.add_option("-s", "--save_path"   , action = "store", type = "string", dest = "path_to_save_results", help = 'Path to save results')
    (options, args) = parser.parse_args()

    # Check if we have an expression matrix filea FASTA or twobit file is given as input.
    if ( (options.pos_data_path is None) or (options.neg_data_path is None) or (options.path_to_save_results is None)):
        parser.print_help()
        print >> sys.stderr, '\nERROR: minimum required options not satisfied:\n'
        sys.exit(1)

    PATH_TO_DATA_POS = options.pos_data_path
    PATH_TO_DATA_NEG = options.neg_data_path
    PATH_TO_SAVE = options.path_to_save_results
    ###----Read file with posotives----###
    data_for_training_pos = pd.read_table( PATH_TO_DATA_POS, sep = "\t", header = 0, index_col = 0 )
    data_for_training_neg = pd.read_table( PATH_TO_DATA_NEG, sep = "\t", header = 0, index_col = 0 )

    ###---Get feature order---##
    feature_order_pos = data_for_training_pos.columns.values.tolist()
    feature_order_neg = data_for_training_neg.columns.values.tolist()

    if feature_order_pos == feature_order_neg:
        filename_to_save_fig = PATH_TO_SAVE
        ###----Append data frames----###
        class_labels = [1]*data_for_training_pos.shape[0] + [0]*data_for_training_neg.shape[0]
        data_for_training_pos_neg = data_for_training_pos.append(data_for_training_neg)
        class_labels = np.array(class_labels)
        data_for_training_pos_neg = pd.np.array(data_for_training_pos_neg)
        calculate_feature_importance(data_for_training_pos_neg, feature_order_pos, class_labels,filename_to_save_fig)




















