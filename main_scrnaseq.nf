#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////
  
MD_ANACONDA = "Anaconda3/2020.07"

PY_DIR = params.WD
//OUT_DIR = Paths.get( params.WD , "Data_interim_files_2" ).toString()
//OUT_SCV = Paths.get( params.WD , "scviewer_files" ).toString()

PY_DES_SCRIPT = Paths.get(PY_DIR,"scrub_things.py").toString()
R0=Paths.get(PY_DIR,"Unpack_CellRanger.R").toString()
R3=Paths.get(PY_DIR,"Preprocess_G2Mregress.R").toString()

////////////////////////////////////////
  
include { SCRUBLET } from './workflows/scrublet'
include { LEAD_MARKERS } from './workflows/lead_cluster_markers'
//include { CELL_SPLIT } from './workflows/cells_split'
//
include { integration } from './modules/integration'
include { make_scv_file } from './modules/make_scv_files'



Channel
  .from( "'0.3'", "'0.9'" )
  .set{ RESOL }


///////////////////////////////////////////////////////////////////////////////
//// PROCESSES ////////////////////////////////////////////////////////////////
        
process unpack_CellRanger {
        
//        container "/camp/stp/babs/working/schneid/Singularity/sceasy/sceasy_nourd.sif"
        
	module "Anaconda3/2019.07"
        conda "/camp/stp/babs/working/schneid/conda/envs/R-4.2-Seurat"

        cpus 4
        time "1h"
        memory "10G"

        input:
                path(loc)

        output:
                path("*.RDS")
        script:
        """
        Rscript ${R0} $loc
        """
}

//// Norm + PCA +Umap ////////////////////////////////////////////////////////

process process_1 {
        
        module "Anaconda3/2019.07"
        conda "/camp/stp/babs/working/schneid/conda/envs/R-4.2-Seurat"
        
	      cpus 1
        time "2h"
        memory "100G"
        
        publishDir params.outdir,
        	mode: "copy",
        	overwrite: true
        	
        input:
                path(Rds)

        output:
                path("*_Umapped_Filteredseurat_object.RDS"), emit: main
		            path("*_Umapped_seurat_object.RDS"), emit: sub
        script:
                """
                Rscript ${R3} $Rds
                """
}

////////////////////////////////////////

workflow {

  unpack_CellRanger(params.WD)
  
  SCRUBLET(unpack_CellRanger.out)
  
  process_1(SCRUBLET.out)
  integration(process_1.out.main) 
  make_scv_file(integration.out)
  
  LEAD_MARKERS(
      RESOL,
      integration.out
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
