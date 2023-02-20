#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////
  
MD_ANACONDA = "Anaconda3/2020.07"

////////////////////////////////////////
  
include { SCRUBLET } from './workflows/scrublet'
include { CELLCYCLE_REGRESS } from './workflows/cell_cycle_regress'
include { LEAD_MARKERS } from './workflows/lead_cluster_markers'
include { CELLS_SPLIT } from './workflows/cells_separate'
//
include { integration } from './modules/integration'
include { make_scv_file } from './modules/make_scv_files'
include { assign_identities } from './modules/assign_cellIDs'

Channel
  .from( "'0.3'", "'0.9'" )
  .set{ RESOL }


///////////////////////////////////////////////////////////////////////////////
//// PROCESSES ////////////////////////////////////////////////////////////////
        
process unpack_CellRanger {
        
	label 'save_output'

        container "/camp/stp/babs/working/schneid/Singularity/sceasy/sceasy_nourd.sif"
        
        cpus 4
        time "1h"
        memory "10G"

        input:
                path(loc)

        output:
                path("*.RDS")
        script:
        """
        Unpack_CellRanger.R $loc
        """
}

//// Norm + PCA +Umap ////////////////////////////////////////////////////////

process process_1 {
        
        label 'save_output'

	module params.MD_ANACONDA
        conda params.CONDA_ENV 
        
	cpus 4
        time "6h"
        memory "100G"
        
        input:
                path(Rds)

        output:
                path("*_Umapped_Filteredseurat_object.RDS"), emit: main
		path("*_Umapped_seurat_object.RDS"), emit: sub
        script:
                """
                Preprocess_G2Mregress.R $Rds
                """
}

process integration_indv {

        label 'save_output'

        module params.MD_ANACONDA
        conda params.CONDA_ENV

        cpus 4
        time "6h"
        memory "100G"

        input:
                path(Rds)

        output:
                path("SC21137_Integrated_Filtered.RDS")
        
	script:
        """
        Integration_from_indv.R $Rds
        """


}

////////////////////////////////////////

workflow {

  unpack_CellRanger(params.WD)
  
  SCRUBLET(unpack_CellRanger.out)
  
 // process_1(SCRUBLET.out)

  CELLCYCLE_REGRESS(SCRUBLET.out)
  integration_indv(CELLCYCLE_REGRESS.out) | assign_identities

  make_scv_file(assign_identities.out)
  
  LEAD_MARKERS(
      RESOL,
      assign_identities.out
  )
  
  CELLS_SPLIT(
      RESOL,
      assign_identities.out
  )
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  THAT'S IT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
