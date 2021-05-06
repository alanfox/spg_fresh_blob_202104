#!/usr/bin/env bash
#SBATCH --job-name=037_afox_1990-2019_rerun_2nd_time
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=35G
#SBATCH --time=30:00:00
#SBATCH --partition=cluster
#SBATCH --mail-user=chschmidt@geomar.de
#SBATCH --mail-type=ALL

module load singularity/3.5.2 

for date in '2003-08-21' '2000-03-09'; do
    sleep 0.3  # spread out a little to avoid a bug resulting from a race
               # condition for an unused port
               # see https://github.com/jupyter/jupyter_client/issues/487

    srun --ntasks=1 --exclusive singularity run -B /sfs -B /gxfs_work1 -B $PWD:/work --pwd /work \
        containers/parcels-container_2021.03.17-6c459b7.sif bash -c \
            "cd notebooks; . /opt/conda/etc/profile.d/conda.sh && conda activate base && papermill 037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill.ipynb executed/037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill_executed_${date}.ipynb -k python -p path_name /gxfs_work1/geomar/smomw355/model_data/ocean-only/VIKING20X.L46-KKG36107B/nemo/output/ -p data_resolution 5d -p w_name_extension '' -p mask_path_name /gxfs_work1/geomar/smomw355/model_data/ocean-only/VIKING20X.L46-KKG36107B/nemo/suppl/ -p mesh_mask_filename 1_mesh_mask.nc -p year_prefix '' -p runtime_in_days 3650 -p create_number_particles 4000000 -p use_number_particles 4000000 -p max_release_depth 1000 -p max_current 2.0 -p t_0_str '1980-01-03T12:00:00' -p t_start_str '${date}T12:00:00' -p use_dask_chunks False" &
done

wait

jobinfo
