#!/bin/bash

names='sampleA'

#qsub_output="output/$patient"
#qsub_error="error/$patient"
for patient in $names
do
	#qsub -l nodes=1,mem=15gb,walltime=0:30:00 -o "$qsub_output" -e "$qsub_error" -v patient=$patient ../integration/run_label_prop.pbs
	sbatch --mem=15gb --time=00:30:00 --output="/zfs3/users/grahamemma9/grahamemma9/March_2018_metabolomic_data/metPropagate/label_propagation/output/$names" --error="/zfs3/users/grahamemma9/grahamemma9/March_2018_metabolomic_data/metPropagate/label_propagation/error/$names" --export=patient=$names /zfs3/users/grahamemma9/grahamemma9/March_2018_metabolomic_data/metPropagate/integration/run_label_prop.pbs

done

