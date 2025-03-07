#! /usr/bin/env bash

script=/home/users/nus/e1083772/fhr/pipeline/06_fit_plsda.R

module load gcc
module load r/4.2.0

shuffle=$((PBS_ARRAY_INDEX / 5 + 1))
fold=$((PBS_ARRAY_INDEX % 5 + 1))
echo "shuffle: $shuffle, fold: $fold"
Rscript $script $shuffle $fold