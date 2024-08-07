#!/bin/bash

#Usage: runTransWarp.sh runEventLoopMC.root warped.root varName

VARIABLE=$3 #pTmu
MIGRATION_FILE=$1
TRUE_HIST=${VARIABLE}_efficiency_numerator
WARPED_FILE=$2
RECO_HIST=${VARIABLE}_selected_signal_reco
STAT_SCALE=$4
UNC_FACTOR=$5

OUTFILE_NAME=$(basename $2)

TransWarpExtraction --stat_scale $STAT_SCALE --corr_factor $UNC_FACTOR --output_file Warping_${VARIABLE}_StatScale_${STAT_SCALE}_UncFactor_${UNC_FACTOR}_$OUTFILE_NAME --data $RECO_HIST --data_file $WARPED_FILE --data_truth $TRUE_HIST --data_truth_file $WARPED_FILE --migration ${VARIABLE}_migration --migration_file $MIGRATION_FILE --reco $RECO_HIST --reco_file $MIGRATION_FILE --truth $TRUE_HIST --truth_file $MIGRATION_FILE --num_iter 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,40,50,60,70,80,90,100 --num_uni 100 -C 2
