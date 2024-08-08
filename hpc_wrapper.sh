#!/bin/bash

# Array of sd values
sd_array=("2.5" "3" "1000")
file_array=("no_effect_simulation.R" "small_effect_simulation.R" "large_effect_simulation.R")
runs=100
participants=200
trials_per_condition=80


for ((i=0; i<3; i++)); do

	current_sd=${sd_array[i]}

  for ((j=0; j<3; j++)); do
	  current_file=${file_array[j]}
	  sbatch << EOF
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --partition=hpc2019
#SBATCH --job-name=cse_simul_sd_${current_sd}_${current_file}
#SBATCH --mem=32000
#SBATCH --output=log_cse_simul_sd_${current_sd}_${current_file}.log
#SBATCH --error=err_cse_simul_sd_${current_sd}_${current_file}.err

cd /users/pmxtny/cse_simul
source /users/pmxtny/.bashrc
export PATH=/users/pmxtny/local/bin:$PATH
export LD_LIBRARY_PATH=/users/pmxtny/local/lib:$LD_LIBRARY_PATH

Rscript ${current_file} ${runs} ${participants} ${trials_per_condition} ${current_sd}
EOF
  done
done
