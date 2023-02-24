import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////
  

///////////////////////////////////////////////////////////////////////////////
//// ANNDATA + SCRUBLET ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process seurat_anndata {
        
        container "/camp/stp/babs/working/schneid/Singularity/sceasy/sceasy_nourd.sif"
//        container "quay.io/biocontainers/r-sceasy:0.0.7--r42hdfd78af_1"

        cpus 4
        time "1h"
        memory "10G"

        input:
                path(RDSobj)

        output:
                path("*_raw.h5ad")
        script:
        """
	seurat_to_AnnData.R $RDSobj
	"""

}


process scrublet {

        executor "slurm"
        cpus 4
        time "3h"
        memory "10G"
 //       module params.MD_ANACONDA
 //       conda params.CONDA_ENV_PY
        
        container "docker://saphir746/scrublet:v1"

        tag { sample }

	label 'save_h5s'

        input:
                tuple val(sample), path(h5d)

        output:
                tuple val(sample), path("${sample}*_scrubletted.h5ad")

        script:
        """
        scrub_things.py $h5d
        """
}

process anndata_seurat {
        
//        container "quay.io/biocontainers/r-sceasy:0.0.7--r42hdfd78af_1"
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
        AnnData_to_seurat.R
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
