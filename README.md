# PRIME: Predicted Regulatory Impact of a Mutation in an Enhancer



## Requirements:

1. Python 2.7
2. NumPy (>= 1.6.1)
3. SciPy (>= 0.9)
4. scikit-learn (>=14.0)
5. bx.python(>=0.7)
6. Pandas
7. Matplotlib

Istalling from PyPI:


```
pip install -U numpy scipy scikit-learn bx-python pandas matplotlib
```

Then you can clone the code:

```
git clone https://github.com/aertslab/primescore.git
```

## Description of the folders:

1. ./scr: folder with scripts
2. ./cbust_folder: folder with Cluster-buster tool
3. ./pwms_folder: folder with PWMs (Position Weight Matrix)
4. ./example: folder with example data
5. ./results: folder with results


## Step 1. Feature-vector representation of the positives and negatives

Make table with feature values of the data.
At this step you represent your fasta sequences as feature-vectors. 
For each sequence in the fasta file will be assigned set of scores and number of the features correspond to the number of used PWMs (PATH_TO_FEATURE_ORDER). 
```bash
# for positives
# file with names of the PWMs (correspond to the name of the PWM in the folder ./pwms_folder)
PATH_TO_FEATURE_ORDER=./example/feature_order.txt 
#.fasta file
PATH_TO_FASTA_FILE=./example/pos.fa
# path to save results
PATH_TO_SAVE=./results/pos.FT.txt
# python script
PATH_TO_SCRIPT=./src/make_feature_table.py
# folder where files with PWMs are
PATH_TO_SINGLEOTNS=./pwms
# path to compiled Cluster-buster tool
PATH_TO_CBUST=./cbust_folder/cbust
# command to run the script
python ${PATH_TO_SCRIPT} -f ${PATH_TO_FASTA_FILE} -o ${PATH_TO_FEATURE_ORDER} -m ${PATH_TO_SINGLEOTNS}/ -s ${PATH_TO_SAVE} -c ${PATH_TO_CBUST}

# for negatives
# file with names of the PWMs (correspond to the name of the PWM in the folder ./pwms_folder)
PATH_TO_FEATURE_ORDER=./example/feature_order.txt
#.fasta file
PATH_TO_FASTA_FILE=./example/neg.fa
# path to save results
PATH_TO_SAVE=./results/neg.FT.txt
# python script
PATH_TO_SCRIPT=./src/make_feature_table.py
# folder where files with PWMs are
PATH_TO_SINGLEOTNS=./pwms
# path to compiled Cluster-buster tool
PATH_TO_CBUST=./cbust_folder/cbust
# command to run the script
python ${PATH_TO_SCRIPT} -f ${PATH_TO_FASTA_FILE} -o ${PATH_TO_FEATURE_ORDER} -m ${PATH_TO_SINGLEOTNS}/ -s ${PATH_TO_SAVE} -c ${PATH_TO_CBUST}
```

## Step 2. Train classifier


At this step will be generated **.pcl** (pickle) file with trained classifier.
Models generated with RF classifiers are always different due to internal training procedure and in order to calculate PRIME score the same classifier shoulde be applied to reference and mutated sequences.


```bash
# name of the file where model will be saved
PATH_TO_SAVE_MODEL=./results/model.pcl
# python script
PATH_TO_SCRIPT=./src/train_model.py
# file for table with feature values for positives
PATH_TO_POSITIVES=./results/pos.FT.txt
# file for table with feature values for negatives
PATH_TO_NEGATIVES=./results/neg.FT.txt
# command to run the script
python ${PATH_TO_SCRIPT} -p ${PATH_TO_POSITIVES} -n ${PATH_TO_NEGATIVES} -s ${PATH_TO_SAVE_MODEL}
```

### Make cross-validation
You can estimate quality of the classifier using cross-validation.

```bash

# file to save average across folds are ander precision-recall and roc corve values
PATH_TO_SAVE_SCORES_PER_SAMPLE=./results/aupr_auroc.txt
# file to save RF score for each region in the test set
PATH_TO_SAVE_AVERAGE_SCORES=./results/rf_score_per_region.txt
# python script
PATH_TO_SCRIPT=./src/make_cross_validation.py
# # file for table with feature values for positives
PATH_TO_POSITIVES=./results/pos.FT.txt
# file for table with feature values for negatives
PATH_TO_NEGATIVES=./results/neg.FT.txt
# command to run
python ${PATH_TO_SCRIPT} -p ${PATH_TO_POSITIVES} -n ${PATH_TO_NEGATIVES} -c ${PATH_TO_SAVE_AVERAGE_SCORES} -s ${PATH_TO_SAVE_SCORES_PER_SAMPLE}
```

## Step 3: Score regions with classifier
To score new sequences with trained classifier use **classify_sequences.py** script

```bash
PATH_TO_TRAINED_CLASSIFIER=./results/model.pcl
# Path to feature-table file
PATH_TO_DATA=./results/pos.FT.txt
# path to save Rf prediction results
PATH_TO_SAVE=./results/pos.scored_with_model.txt
# path to script
PATH_TO_SCRIPT=./src/classify_sequences.py
# command
python ${PATH_TO_SCRIPT} -d ${PATH_TO_DATA} -m ${PATH_TO_TRAINED_CLASSIFIER} -s ${PATH_TO_SAVE}
```


## Step 4: Calculate PRIME score
To claculate PRIME scores the following data should be provided

1. **.2bit** file ( here we Yassume it is in the folder **./example/**
2. File with SNVs

For each mutation (SNV, insertion or deletion) script checks 10 regions (windows) around mutations. 

Output:
```
chr start end ref mut	ID	ref_score	mut_score
```

Difference between classification score in mutant versus reference sequence

**PRIME=mut_score - ref_score**

### SNVs

Structure of the file:

```
chr	start	end	ref_nucleotide	mut_nucleotide	ID
```
```bash
# python script to calculate PRIME score for each SNV
PATH_TO_SCRIPT=./src/calculate_PRIME_snvs.py
# file with named of the PWMs (correspond to the name of the PWM in the folder ./pwms_folder)
# use the same file (with the same PWMs and order of the PWMs inside the file) as for making of feature-tables
PATH_TO_FEATURE_ORDER=./example/feature_order.txt
# path to .2bit file 
PATH_TO_2BIT_FILE=./example/hg19.2bit
# path to file with SNVs
PATH_TO_SNVs=./example/snvs.txt
# path to file where you want to save results
PATH_TO_SAVE=./results/snvs_prime_test_2.txt
# path to trained classifier file
PATH_TO_MODEL=./results/model.pcl
# path to folder with PWMs
PATH_TO_SINGLEOTNS=./pwms
# path to compiled Cluster-buster tool
PATH_TO_CBUST=./cbust_folder/cbust
# command to run
python ${PATH_TO_SCRIPT} -f ${PATH_TO_2BIT_FILE} -o ${PATH_TO_FEATURE_ORDER} -b ${PATH_TO_SNVs} -s ${PATH_TO_SAVE} -m ${PATH_TO_MODEL} -w 800 -p ${PATH_TO_SINGLEOTNS}/ -c ${PATH_TO_CBUST}
```


### Insertions
Structure of the file for insertions:

```
chr	start	ref_nucleotide ref_nucleotide_insertion	ID
chr6	161094665	C	CA	LPA
161094665: zero-based coordinate of the nucleotids before the insertion
```

```bash
# python script to calculate PRIME score for each insertion
PATH_TO_SCRIPT=./src/calculate_PRIME_insertions.py
# file with named of the PWMs (correspond to the name of the PWM in the folder ./pwms_folder)
# use the same file (with the same PWMs and order of the PWMs inside the file) as for making of feature-tables
PATH_TO_FEATURE_ORDER=./example/feature_order.txt
# path to .2bit file 
PATH_TO_2BIT_FILE=./example/hg19.2bit
# path to file with insertions
PATH_TO_INSERTIONS=./example/insertions.txt
# path to file where you want to save results
PATH_TO_SAVE=./results/insertions_prime.txt
# path to trained classifier file
PATH_TO_MODEL=./results/model.pcl
# path to folder with PWMs
PATH_TO_SINGLEOTNS=./pwms
# path to compiled Cluster-buster tool
PATH_TO_CBUST=./cbust_folder/cbust
python ${PATH_TO_SCRIPT} -f ${PATH_TO_2BIT_FILE} -o ${PATH_TO_FEATURE_ORDER} -b ${PATH_TO_INSERTIONS} -s ${PATH_TO_SAVE} -m ${PATH_TO_MODEL} -w 800 -p ${PATH_TO_SINGLEOTNS}/ -c ${PATH_TO_CBUST}
```



### Deletions
Structure of the file for deletions
chr6	start	ref_nucleotide_deleted_sequence	ref_nucleotide	ID

```
chr10	100164477	GA	G ID
100164477: zero-based coordinate of the nucleotids before the deletion
```

```bash
# python script to calculate PRIME score for each deletions
PATH_TO_SCRIPT=./src/calculate_PRIME_deletions.py
# file with named of the PWMs (correspond to the name of the PWM in the folder ./pwms_folder)
# use the same file (with the same PWMs and order of the PWMs inside the file) as for making of feature-tables
PATH_TO_FEATURE_ORDER=./example/feature_order.txt
# path to .2bit file 
PATH_TO_2BIT_FILE=./example/hg19.2bit
# path to file with insertions
PATH_TO_INSERTIONS=./example/deletions.txt
# path to file where you want to save results
PATH_TO_SAVE=./results/deletions_prime.txt
# path to trained classifier file
PATH_TO_MODEL=./results/model.pcl
# path to folder with PWMs
PATH_TO_SINGLEOTNS=./pwms
# path to compiled Cluster-buster tool
PATH_TO_CBUST=./cbust_folder/cbust
python ${PATH_TO_SCRIPT} -f ${PATH_TO_2BIT_FILE} -o ${PATH_TO_FEATURE_ORDER} -b ${PATH_TO_INSERTIONS} -s ${PATH_TO_SAVE} -m ${PATH_TO_MODEL} -w 800 -p ${PATH_TO_SINGLEOTNS}/ -c ${PATH_TO_CBUST}
```


## Additional scripts:

### Select features (PWMs) for model 


If you are using i-cis-Target tool for motif enrichment anaysis we provide additional script to select PWMs (features) for your RF model.

To select motifs the following data are necessery:

1. File with motif enrichment results (from i-cis-Target tool): **statistics.tbl**
2. File with list of all singletones (PWMs) you have (all **.cb** files in the **./pwms_folder**). This files can be found in the folder with examples.



```bash
# python script making selection of PWMs based on i-cis-Target enrichment analysis
PATH_TO_SCRIPT=./src/select_tf_model_features.py
# i-cis-Target enrichment analysis results
PATH_TO_ENRICHMENT_RESULTS=./example/ctx_rcc_enrichment_results.txt
# path to save results: selected PWMs for TF specific RF model
PATH_TO_SAVE=./example/my_selected_pwms.txt
# file with annotation of the motifs (PWM->TF)
PWM2TF_ANNOTATION=./example/pwm2tf.tbl
# text file with list of all avialable PWMs for further scoring
PATH_TO_SINGLEOTNSLIST=./example/pwm_list.txt
# path to save selected features
PATH_TO_SAVE=./example/selected_model_features.txt
# command to run
python ${PATH_TO_SCRIPT}  -s ${PATH_TO_ENRICHMENT_RESULTS} -t ATF2 -a ${PWM2TF_ANNOTATION} -m ${PATH_TO_SINGLEOTNSLIST} > ${PATH_TO_SAVE}
```

## Scripts to make .fasta files with mutations
```bash
# PAth to script
PATH_TO_SCRIPT=./src/mut2fasta.py
# Path to save fasta with mutations
PATH_TO_SAVE_REF_FASTA=./results/ref.fa
# Path to save refeence fasta sequence
PATH_TO_SAVE_MUT_FASTA=./results/mut.fa

PATH_TO_2BIT_FILE=./example/hg19.2bit
# file with mutations (can me mixed insertions, deletions and SNVs in the same file)
PATH_TO_MUTATIONS=./example/mutations.txt
python ${PATH_TO_SCRIPT} -r ${PATH_TO_SAVE_REF_FASTA} -m ${PATH_TO_SAVE_MUT_FASTA} -f ${PATH_TO_2BIT_FILE} -v ${PATH_TO_MUTATIONS}
```

This **fasta** files coulde be used to make feature-tables (apply **make_feature_table.py** ) and then be scored with **classify_sequences.py**. To calculate **PRIME** score you can substract mutant and reference classifier scores. 
This approach is faster comparing to usage of **calculate_PRIME<>.py** scripts because only one region is scored (mutation will be in the center of the window). But when using **calculate_PRIME<>.py**  **10** windows around each mutation are checked. 


## References



 * [i-cis-Target tool](https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/)
 * [Cluster-Buster](http://zlab.bu.edu/cluster-buster/download.html)
 
1. Imrichová,H., Hulselmans,G., Kalender Atak,Z., Potier,D. and Aerts,S. (2015) i-cisTarget 2015 update: generalized cis-regulatory enrichment analysis in human, mouse and fly. Nucleic Acids Res. doi: 10.1093/nar/gkv395

2. Herrmann,C., Van de Sande,B., Potier,D. and Aerts,S. (2012) i-cisTarget: an integrative genomics method for the prediction of regulatory features and cis-regulatory modules. Nucleic Acids Res. doi: 10.1093/nar/gks543

3. Frith MC1, Li MC, Weng Z. Cluster-Buster: Finding dense clusters of motifs in DNA sequences. Nucleic Acids Res. 2003 Jul 1;31(13):3666-8. 
