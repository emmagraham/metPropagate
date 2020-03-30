#!/bin/bash

#name of patients you want to analyze. Separate each with space. 
names='sampleA'

qsub_output="output/$names"
qsub_error="error/$names"

#qsub -l nodes=6,mem=50gb,walltime=20:00:00 -o "$qsub_output" -e "$qsub_error" -v patient=$names do_met_processing.pbs
sbatch --mem=50gb --time=20:00:00 --output="/zfs3/users/grahamemma9/grahamemma9/March_2018_metabolomic_data/metPropagate/METABOLOMIC_PROCESSING_PIPELINE/output/$names" --error="/zfs3/users/grahamemma9/grahamemma9/March_2018_metabolomic_data/metPropagate/METABOLOMIC_PROCESSING_PIPELINE/error/$names" --export=patient=$names /zfs3/users/grahamemma9/grahamemma9/March_2018_metabolomic_data/metPropagate/METABOLOMIC_PROCESSING_PIPELINE/do_met_processing.pbs
echo 'end script'

	