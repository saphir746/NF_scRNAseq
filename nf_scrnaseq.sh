#!/bin/sh

module load Nextflow/22.04.0
module load Singularity/3.6.4

WORK_DIR=/camp/stp/babs/scratch/schneid/SC21137_CAF_HPV_G2M/

PROJECT_dir=/camp/stp/babs/working/schneid/projects/sahaie/giovanni.giangreco/Characterisation_of_CAF_in_HPV_cancer_scrnaseq/
OUT_DIR=${PROJECT_dir}Data_interim_files_2
OUT_SCV=${PROJECT_dir}scviewer_files_G2M

nextflow run main_scrnaseq.nf -resume \
                        --WD ${PROJECT_dir} \
			--outdir ${OUT_DIR} \
			--outscv ${OUT_SCV} \
			-profile cluster \  ## local \
                        -work-dir $WORK_DIR
