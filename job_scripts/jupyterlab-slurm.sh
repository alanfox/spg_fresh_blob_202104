#!/usr/bin/env bash

#
# Login node usage
#
# $ ./jupyterlab-slurm.sh
#
# Compute node usage
#
# $ sbatch jupyterlab-slurm.sh  # request Jupyter server on a compute node
# $ ./jupyterlab-slurm.sh <jobid>  # fetch Jupyter server address
# $ scancel <jobid>  # free resources, after you are done!
#
# Compute node usage options
#
# $ sbatch --time=<HH:MM:SS> jupyterlab-slurm.sh
# $ sbatch --cpus-per-task=<number-of-CPUs> jupyterlab-slurm.sh
# $ sbatch --mem=<memory-size> jupyterlab-slurm.sh
# $ sbatch --partition=<batch-class> jupyterlab-slurm.sh
# $ sbatch --account=<compute-project> jupyterlab-slurm.sh
#

#
#SBATCH --cpus-per-task=32
#SBATCH --mem=1500GB
#SBATCH --time=24:00:00
##SBATCH --partition=
##SBATCH --account=
#
#SBATCH --nodes=1 --ntasks-per-node=1  # change these only if you exactly know what you are doing!
#SBATCH --job-name=jupyterlab
#SBATCH --output=slurm-jupyterlab-%j.log
#SBATCH --error=slurm-jupyterlab-%j.log
#

set -e  # for triggering SLURM job states correctly

if [ -n "${1}" ]; then

  scontrol show job "${1}" &> /dev/null || { echo 'Job does not exist...'; exit; }
  STATUS=$(scontrol show job "${1}" | grep 'JobState' | awk '{print $1}')

  if [ "${STATUS}" == 'JobState=PENDING' ] || [ "${STATUS}" == 'JobState=CONFIGURING' ]; then

    scontrol show job "${1}" | grep 'StartTime' | awk '{print $1}'

  elif [ "${STATUS}" == 'JobState=RUNNING' ]; then

    END_TIME=$(scontrol show job "${1}" | grep 'EndTime' | awk '{print $2}')
    END_TIME_SEC=$(date -d "${END_TIME##*=}" +"%s")
    REMAINING_TIME_SEC=$((END_TIME_SEC-$(date +"%s")))
    TZ=Berlin printf RemainingTime="%d-%(%H:%M:%S)T " $((REMAINING_TIME_SEC/86400)) ${REMAINING_TIME_SEC}; echo "${END_TIME}"
    JOB_STD_OUT=$(scontrol show job "${1}" | grep 'StdOut' | awk '{print $1}')
    JOB_STD_OUT=${JOB_STD_OUT##*=}
    REMOTE_NODE_LOCATION=$(cat "${JOB_STD_OUT}" | grep node_network_location | cut -d ' ' -f 2)
    cat "${JOB_STD_OUT}" | grep -m 1 -e "://${REMOTE_NODE_LOCATION}:" || echo 'JupyterLab not yet properly running... please wait.'

  elif [ "${STATUS}" == 'JobState=COMPLETING' ] || [ "${STATUS}" == 'JobState=COMPLETED' ]; then

    echo 'JupyterLab has already been shutdown...'

  else

    echo "JupyterLab server not running... ${STATUS}"

 fi

else
	
  ##### Specify base directory for serving Jupyter notebooks. #####

  # cd "${HOME}"
  echo "Jupyter notebooks are served from: ${PWD}"

  ##### Specify a node network location that is visible from the login nodes. #####

  VISIBLE_NETWORK_LOCATION=$(hostname) # NESH, HLRN-G, ...
  #VISIBLE_NETWORK_LOCATION=$(ifconfig ib0 | grep "inet\ " | awk '{print $2}') # JUWELS, ...

  echo "node_network_location ${VISIBLE_NETWORK_LOCATION}"

  ##### Launch the Jupyter server. #####

  # This is necessary to properly initialize conda inside non-interactive Bash sessions.
  # Works only for conda environments >=v4.6 that were set-up with a `conda init bash` command.
  # For details see: https://github.com/conda/conda/issues/7980

  # eval "$(conda shell.bash hook)"

  # If conda was not setup with `conda init bash` use e.g. this:
  #source "${HOME}"/miniconda3/etc/profile.d/conda.sh

  # Activate Jupyter environment.
  # conda activate base

  # Launch Jupyter server.
  module load singularity/3.5.2
  singularity run -B /sfs -B /gxfs_work1 -B $PWD:/work --pwd /work parcels-container_2021.09.29-09ab0ce.sif jupyter lab --ip="${VISIBLE_NETWORK_LOCATION}" --no-browser

  ##### Report usage of resources to guide adapting the job script header. #####

  jobinfo || echo '' # NESH-only

fi

