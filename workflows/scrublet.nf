#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////
  
MD_ANACONDA = "Anaconda3/2020.07"
CONDA_ENV = '/camp/stp/babs/working/schneid/conda/envs/Scrublet_py'

PY_DIR = workflow.projectDir.toString()
OUT_DIR = Paths.get( workflow.projectDir.toString() , "Data_interim_files" ).toString()
OUT_SCV = Paths.get( workflow.projectDir.toString() , "scviewer_files" ).toString()
PY_DES_SCRIPT = Paths.get(PY_DIR,"scrub_things.py").toString()
R1=Paths.get(PY_DIR,"seurat_to_AnnData.R").toString()
R2=Paths.get(PY_DIR,"AnnData_to_seurat.R").toString()

///////////////////////////////////////////////////////////////////////////////
//// ANNDATA + SCRUBLET ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process seurat_anndata {
        
        container "/camp/stp/babs/working/schneid/Singularity/sceasy/sceasy_nourd.sif"
        cpus 4
        time "1h"
        memory "10G"

        input:
                path(RDSobj)

        output:
                path("*_raw.h5ad")
        script:
                """
		          Rscript ${R1} $RDSobj
	           	"""

}


process scrublet {

        executor "slurm"
        cpus 4
        time "3h"
        memory "10G"
        module MD_ANACONDA
        conda CONDA_ENV

        tag { sample }

        input:
                tuple val(sample), path(h5d)

        output:
                tuple val(sample), path("${sample}*_scrubletted.h5ad")

        script:
                script = Paths.get(PY_DES_SCRIPT).toString()
                """
                python ${script} $h5d
                """
}

process anndata_seurat {
        
        container "/camp/stp/babs/working/schneid/Singularity/sceasy/sceasy_nourd.sif"
        
	cpus 1
        time "1h"
        memory "20G"
        
        input:
                path(h5ds)

        output:
                path("*_scrubletted.RDS")
        script:
                """
                Rscript ${R2}
                """
}



///////////////////////////////////////////////////////////////////////////////

workflow SCRUBLET {
  take: cellranger
  
  main:
    seurat_anndata(cellranger)
    seurat_anndata.out
		 .flatten()
   	        .map{[
                it.toString().replaceAll("(.*)/(.*)_raw.h5ad", "\$2"),
                it]}
		 .set{RAWh5ad}

   scrublet(RAWh5ad) 
   scrublet.out
	    .map{[ it[1] ]}
	    .collect()
	    .set{Scrubbedh5ad}

  anndata_seurat(Scrubbedh5ad)
    
  emit: 
    anndata_seurat.out
}
