#!/usr/bin/env bash
#SBATCH --job-name=037_afox_1990-1997
#SBATCH --ntasks=584
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=35G
#SBATCH --time=30:00:00
#SBATCH --partition=cluster
#SBATCH --mail-user=chschmidt@geomar.de
#SBATCH --mail-type=ALL

module load singularity/3.5.2 

for year in {1990..1997}; do
    for month_day in '01-03' '01-08' '01-13' '01-18' '01-23' '01-28' '02-02' '02-07' '02-12' '02-17' '02-22' '02-28' '03-04' '03-09' '03-14' '03-19' '03-24' '03-29' '04-03' '04-08' '04-13' '04-18' '04-23' '04-28' '05-03' '05-08' '05-13' '05-18' '05-23' '05-28' '06-02' '06-07' '06-12' '06-17' '06-22' '06-27' '07-02' '07-07' '07-12' '07-17' '07-22' '07-27' '08-01' '08-06' '08-11' '08-16' '08-21' '08-26' '08-31' '09-05' '09-10' '09-15' '09-20' '09-25' '09-30' '10-05' '10-10' '10-15' '10-20' '10-25' '10-30' '11-04' '11-09' '11-14' '11-19' '11-24' '11-29' '12-04' '12-09' '12-14' '12-19' '12-24' '12-29'; do
        sleep 0.3  # spread out a little to avoid a bug resulting from a race
                   # condition for an unused port
                   # see https://github.com/jupyter/jupyter_client/issues/487

        srun --ntasks=1 --exclusive singularity run -B /sfs -B /gxfs_work1 -B $PWD:/work --pwd /work \
            containers/parcels-container_2021.03.17-6c459b7.sif bash -c \
                "cd notebooks; . /opt/conda/etc/profile.d/conda.sh && conda activate base && papermill 037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill.ipynb executed/037_afox_RunParcels_TS_MXL_Multiline_Randomvel_Papermill_executed_${year}-${month_day}.ipynb -k python -p path_name /gxfs_work1/geomar/smomw355/model_data/ocean-only/VIKING20X.L46-KKG36107B/nemo/output/ -p data_resolution 5d -p w_name_extension '' -p mask_path_name /gxfs_work1/geomar/smomw355/model_data/ocean-only/VIKING20X.L46-KKG36107B/nemo/suppl/ -p mesh_mask_filename 1_mesh_mask.nc -p year_prefix '' -p runtime_in_days 3650 -p create_number_particles 4000000 -p use_number_particles 4000000 -p max_release_depth 1000 -p max_current 2.0 -p t_0_str '1980-01-03T12:00:00' -p t_start_str '${year}-${month_day}T12:00:00' -p use_dask_chunks False" &
    done
done

wait

jobinfo
