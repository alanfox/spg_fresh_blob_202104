#!/usr/bin/env bash
#SBATCH --job-name=037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=12:00:00
#SBATCH --partition=cluster

module load singularity/3.5.2 

srun --ntasks=1 --exclusive singularity run -B /sfs -B /gxfs_work1 -B $PWD:/work --pwd /work \
    parcels-container_2021.03.02-7088d90.sif bash -c \
        ". /opt/conda/etc/profile.d/conda.sh && conda activate base && papermill 037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill.ipynb 037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill.executed_job.ipynb -k python -p path_name \"/gxfs_work1/geomar/smomw355/model_data/ocean-only/VIKING20X.L46-KKG36107B/nemo/output/\" -p mask_path_name \"/gxfs_work1/geomar/smomw355/model_data/ocean-only/VIKING20X.L46-KKG36107B/nemo/suppl/\" -p w_name_extension \"\" -p mesh_mask_filename \"1_mesh_mask.nc\""

jobinfo
