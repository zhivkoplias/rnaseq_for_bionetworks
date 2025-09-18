#!/bin/bash

# Job allocation
SBATCH -A 2020-32 # This is your allocation number, you get it from https://supr.snic.se/ and you need to be in project.

# The name of the script is myjob
SBATCH -J SalmonControls

# Only 24 hours (1 days) wall-clock time will be given to this job
SBATCH -t 128:00:00

# Number of nodes (computers)
SBATCH --nodes=1

# Number of processes per node (24 is recommended for most cases)
# 48 is the default to allow the possibility of hyperthreading
SBATCH --ntasks-per-node=48

SBATCH -e error_file_controls.e
SBATCH -o output_file_controls.o

# Go to your home directory on klemmming (no space limit, but no backups)
#cd /cfs/klemming/nobackup/m/mateuszk/InParanoid
# Load anaconda
module load anaconda
# Load environment 
conda activate salmon_env
# Adding something to PATH must be done after activating new environment
#export PATH=/cfs/klemming/nobackup/m/mateuszk/InParanoid/blast-2.2.18/bin/:$PATH
# Start your script, when this ends, node is released
start=`date +%s`
for i in /cfs/klemming/nobackup/e/erikzhi/ENCODE/controls/*.bam;
do
    echo "Processing sample $i"
    sample=`basename ${i}`
    out_dir=$'/cfs/klemming/nobackup/e/erikzhi/ENCODE/salmon_output/controls/'
    STARTTIME=$(date +%s)
    bamToFastq -i $i -fq ${i}.fq1.fq -fq2 ${i}.fq2.fq
    echo "Salmon part..."
    salmon quant -i /cfs/klemming/nobackup/e/erikzhi/ENCODE/human_index -l A -1 ${i}.fq1.fq -2 ${i}.fq2.fq --validateMappings -o $out_dir$sample
    echo "DONE"
    rm ${i}.fq1.fq
    rm ${i}.fq2.fq
    ENDTIME=$(date +%s)
    echo "It took $((($ENDTIME - $STARTTIME)/60)) minutes to complete"
done
echo "all done"
