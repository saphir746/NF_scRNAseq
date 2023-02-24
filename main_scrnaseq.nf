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
include { assign_identities; transfer_identities } from './modules/assign_cellIDs'

Channel
  .from( "'0.3'", "'0.9'" )
  .set{ RESOL }

Channel
  .from( "'0.1'", "'0.5'","'0.7'" )
  .set{ RESOL_split }

Channel 
  .fromPath('./assets/cc_genes_mm10.rds')
  .set{ CC_GENES }

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
                tuple path(Rds), path(ccgenes)

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
  CELLCYCLE_REGRESS.out.combine(CC_GENES) | integration_indv | transfer_identities

  make_scv_file(transfer_identities.out)
  
//  LEAD_MARKERS(
//      RESOL,
//      transfer_identities.out
//  )
//  
//  CELLS_SPLIT(
//      RESOL_split,
//      transfer_identities.out
//  )
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
