#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//// SPLIT + REGRESS +COMBINE  ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process seurat_split {

        module params.MD_ANACONDA
        conda params.CONDA_ENV

        cpus 4
        time "1h"
        memory "10G"

        input:
                path(RDSobj)

        output:
                path("*.RDS")
        script:
        """
        split_raw_indv_samples.R $RDSobj
        """

}

process G2MREGRESS {

        executor "slurm"
        cpus 4
        time "12h"
        memory "100G"
        module params.MD_ANACONDA
        conda params.CONDA_ENV

        tag { sample }

	label 'save_output'

        input:
                tuple val(sample), path(rds)

        output:
                tuple val(sample), path("SC21137_Umapped_*.RDS")

        script:
        """
        Preprocess_G2Mregress_indv.R $rds
        """
}

process filter_combine {

        module params.MD_ANACONDA
        conda params.CONDA_ENV

	label 'save_output'

        cpus 4
        time "1h"
        memory "100G"

        input:
                path(rdss)

        output:
                path("SC21137_Umapped_Filteredseurat_object.RDS")
        script:
        """
        filter_combine_indv_samples.R
        """
}


///////////////////////////////////////////////////////////////////////////////

workflow CELLCYCLE_REGRESS {
  take: rds

  main:
    seurat_split(rds)   
    seurat_split.out
                .flatten()
                .map{[
                it.toString().replaceAll("(.*)/SC21137_raw_scrubletted_(.*).RDS", "\$2"),
                it]}
                .set{RAWrds}

   G2MREGRESS(RAWrds)
   G2MREGRESS.out
            .map{[ it[1] ]}
            .collect()
            .set{UMAPed}

  filter_combine(UMAPed)

  emit:
    filter_combine.out
}


