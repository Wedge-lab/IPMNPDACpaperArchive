#!/bin/bash --login
#SBATCH -J bbtmpaca
#SBATCH -p multicore
#SBATCH -n 2
#SBATCH -t 2-0
#SBATCH -a 1-2
#SBATCH -o /yourFolder/logs/%A_%a.out
#SBATCH -e /yourFolder/logs/%A_%a.err

module purge
conda activate xtpy311

# Get the Nth line from a file based on array task ID
sampleID=`sed -n "${SLURM_ARRAY_TASK_ID}p" <(awk '{ print $1 }' sampleList.txt)`
ssDPCIfolderpath='/Data/dpcInputs'
dpcoutfolderpath='/Data/dpcOutputs/'
outpath4in='/Data/timingOutputs/'

python3 ./bbDPCtiming.py ${sampleID} ${ssDPCIfolderpath} ${dpcoutfolderpath} ${outpath4in}


















