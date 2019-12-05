#!/bin/bash

RANDOM_SEED=$1
OUTPUT_PATH=$2

matlab -nodisplay -nosplash -r "example1($RANDOM_SEED, '$OUTPUT_PATH'); exit"