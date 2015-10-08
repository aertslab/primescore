#!/usr/bin/env bash

# Author: Bram Van de Sande

########################################################################################################################
# DEPENDENCIES
########################################################################################################################
SCRIPT_FOLDERNAME=$(dirname $0)
CBUST_COMMAND=$(which cbust)
if [ -z "${CBUST_COMMAND}" ]; then
  echo "The executable for cbust is not present in the PATH." > /dev/stderr
  exit 1
fi
TWOBITTOFA_COMMAND=$(which twoBitToFa)
if [ -z "${TWOBITTOFA_COMMAND}" ]; then
  echo "The executable for twoBitToFa is not present in the PATH." > /dev/stderr
  exit 1
fi
if [ ! -e ${SCRIPT_FOLDERNAME}/parseCbust.pl ]; then
  echo "The script depends on parseCbust.pl, which must be in the same directory as this script." > /dev/stderr
  exit 1
fi

########################################################################################################################
# FUNCTIONS
########################################################################################################################

# Create FASTA file from BED file ...
function bedToFasta() {
    local BED_FILENAME=$1
    local GENOME_2BIT_FILENAME=$2
	${TWOBITTOFA_COMMAND} -bed=${BED_FILENAME} ${GENOME_2BIT_FILENAME} stdout
}

# Score via CLusterBuster ...
function crm_locations() {
    local FASTA_FILENAME=$1
    local BED_FILENAME=$2
    local CRM_FILENAME=$3
    local INCLUDE_MOTIFS=$4 # "Y" or "N" ...
    local EXTENSION=$5      # in base pairs ...
	local MOTIF_SCORE_THRESHOLD=$6
	local CLUSTER_SCORE_THRESHOLD=$7

    local TF_NAME=$(basename ${CRM_FILENAME%.cb})
    # Caveat:
    # 1. The locations of the regions are given in 0-based half-open intervals (end position is not included in the
    # interval). In contrast, the CRM locations provided by ClusterBuster are provided as 1-based closed intervals.
    # 2. For sequences located on the negative strand, twoBitToFa (from the UCSC Kent Toolkit) returns the reverse
    # complemented nucleotide sequence. The CRM locations provided by ClusterBuster on these negative strands not
    # be converted.
    # 3. The strand of a CRM or motif is derived based on the following rules:
    # + CRM/motif & + region => + CRM/motif
    # - CRM/motif & + region => - CRM/motif
    # + CRM/motif & - region => - CRM/motif
    # - CRM/motif & - region => + CRM/motif
    if [ ${INCLUDE_MOTIFS} == "Y" ]; then
        ${CBUST_COMMAND} -c${CLUSTER_SCORE_THRESHOLD} -m${MOTIF_SCORE_THRESHOLD} ${CRM_FILENAME} ${FASTA_FILENAME} 2>/dev/null | perl ${SCRIPT_FOLDERNAME}/parseCbust.pl | sort -k1 | \
		    join -o 2.1,2.2,2.3,1.1,1.3,1.4,1.5,1.6,2.6,2.4,1.7 -1 1 -2 4 - ${BED_FILENAME} | \
		    awk -v name=${TF_NAME} -v ext=${EXTENSION} \
		    '$9 == "+" {
		        print $1 "\t" ($2+$6-1)-ext "\t" ($2+$7)+ext "\t" $5 "[" name "@" $10 "]\t" $8 "\t" $11
		     }
		     $9 == "-" {
		        l=($3-$2);
		        if ($11 == "+") s="-"; else s="+";
		        print $1 "\t" ($2+(l-$7))-ext "\t" ($2+(l-$6)+1)+ext "\t" $5 "[" name "@" $10 "]\t" $8 "\t" s
		     }' | \
		    sed "s!misc_feature!motif!"
	else
        ${CBUST_COMMAND} -c${CLUSTER_SCORE_THRESHOLD} -m${MOTIF_SCORE_THRESHOLD} ${CRM_FILENAME} ${FASTA_FILENAME} 2>/dev/null | perl ${SCRIPT_FOLDERNAME}/parseCbust.pl | sort -k1 | \
		    join -o 2.1,2.2,2.3,1.1,1.3,1.4,1.5,1.6,2.6,2.4,1.7 -1 1 -2 4 - ${BED_FILENAME} | \
		    awk -v name=${TF_NAME} -v ext=${EXTENSION} \
		    '$5 == "CRM" && $9 == "+" {
		        print $1 "\t" ($2+$6-1)-ext "\t" ($2+$7)+ext "\t" $5 "[" name "@" $10 "]\t" $8 "\t" $11
		     }
		     $5 == "CRM" && $9 == "-" {
		        l=($3-$2);
		        if ($11 == "+") s="-"; else s="+";
		        print $1 "\t" ($2+(l-$7))-ext "\t" ($2+(l-$6)+1)+ext "\t" $5 "[" name "@" $10 "]\t" $8 "\t" s
		     }' | \
		    sed "s!misc_feature!motif!"
	fi
}

function display_usage() {
	echo
	echo "Usage: cbust.sh [<options>] <bed_filename> <2bit_genome_filename> <cb_filename>"
	echo "Options: -c <thr> : Cluster score threshold (default = 5)"
	echo "         -m <thr> : Motif score threshold (default = 6) "
	echo "         -u       : Normalize scores between range of 0 to 1000 (useful for viewing in UCSC genome browser) "
}

########################################################################################################################
# INPUT PARAMETERS
########################################################################################################################

MOTIF_SCORE_THRESHOLD=6
CLUSTER_SCORE_THRESHOLD=5
NORMALIZE_SCORES=
# -c <thr> Cluster score threshold (5)
# -m <thr> Motif score threshold (6)
# -n Normalize scores to range of 0-1000
while getopts 'hc:m:u' OPTION; do
     case ${OPTION} in
         h)
             display_usage
             exit 1
             ;;
         c)
             CLUSTER_SCORE_THRESHOLD=${OPTARG}
             ;;
         m)
             CLUSTER_SCORE_THRESHOLD=${OPTARG}
             ;;
         u)
             NORMALIZE_SCORES="Yes"
             ;;
         ?)
             display_usage
             exit 2
             ;;
     esac
done

shift $((OPTIND-1))
if [ $# -ne 3 ]; then
	echo "Wrong number of input arguments."
	display_usage
	exit 2
fi
BED_FILENAME=$1
GENOME_2BIT_FILENAME=$2
CRM_FILENAME=$3

########################################################################################################################
# ACTUAL CALCULATIONS
########################################################################################################################

# Adjust BED file: make sure 6 columns are present and file is sorted based on IDs ...
TMP_BED_FILENAME=tmp.$$.input.bed
awk 'NF == 6 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t0\t" $6} NF == 4 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t0\t+"}' ${BED_FILENAME} | sort -k4 > ${TMP_BED_FILENAME}
# Create FASTA file ...
TMP_FASTA_FILENAME=tmp.$$.input.fa
bedToFasta ${TMP_BED_FILENAME} ${GENOME_2BIT_FILENAME} > ${TMP_FASTA_FILENAME}
# Find locations of CRM ...
if [ ! -n "${NORMALIZE_SCORES}" ]; then
	crm_locations ${TMP_FASTA_FILENAME} ${TMP_BED_FILENAME} ${CRM_FILENAME} "Y" 0 ${MOTIF_SCORE_THRESHOLD} ${CLUSTER_SCORE_THRESHOLD}
else
	RESULT_BED_FILENAME=tmp.$$.cbust.bed
	crm_locations ${TMP_FASTA_FILENAME} ${TMP_BED_FILENAME} ${CRM_FILENAME} "Y" 0 ${MOTIF_SCORE_THRESHOLD} ${CLUSTER_SCORE_THRESHOLD} > ${RESULT_BED_FILENAME}
	echo "track name=CRMs description=\"Predicted CRMs and motifs\" useScore=1"
	MIN_CRM_SCORE=$(cat ${RESULT_BED_FILENAME} | grep 'CRM\[' | cut -f5 | sort -gu | head -1)
	MIN_MOTIF_SCORE=$(cat ${RESULT_BED_FILENAME} | grep 'motif\[' | cut -f5 | sort -gu | head -1)

	MAX_CRM_SCORE=$(cat ${RESULT_BED_FILENAME} | grep 'CRM\['| cut -f5 | sort -gu | tail -1)
	MAX_MOTIF_SCORE=$(cat ${RESULT_BED_FILENAME} | grep 'motif\[' | cut -f5 | sort -gu | tail -1)

	awk -v minCRMScore=${MIN_CRM_SCORE} -v maxCRMScore=${MAX_CRM_SCORE} -v minMotifScore=${MIN_MOTIF_SCORE} -v maxMotifScore=${MAX_MOTIF_SCORE} \
				'$4 ~ /CRM*/ { print $1 "\t" $2 "\t" $3 "\t" $4 "\t" (1000.0 * ($5-minCRMScore))/(maxCRMScore-minCRMScore) "\t" $6 }
				 $4 ~ /motif*/ { print $1 "\t" $2 "\t" $3 "\t" $4 "\t" (1000.0 * ($5-minMotifScore))/(maxMotifScore-minMotifScore) "\t" $6}' \
				${RESULT_BED_FILENAME}
fi

# Clean up ...
rm -f ${TMP_BED_FILENAME} ${RESULT_BED_FILENAME} ${TMP_FASTA_FILENAME}

