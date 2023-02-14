#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////


process lead_cluster_markers {
  
  module "Anaconda3/2019.07"
  conda "/camp/stp/babs/working/schneid/conda/envs/R-4.2-Seurat"
  
  cpus 1
  time "6h"
  memory "100G"
  
  input:
    tuple path(Rds), val(res)
  
  output:
    path("*.csv")
  
  script:
    """
    Rscript /camp/stp/babs/working/schneid/projects/sahaie/giovanni.giangreco/Characterisation_of_CAF_in_HPV_cancer_scrnaseq/Find_lead_markers_clusters.R $Rds $res
    """
}


///////////////////////////////////////

workflow LEAD_MARKERS {
  take:
    resol
    rds
  
  main:
    rds.combine(resol) | lead_cluster_markers
  
  emit:
    lead_cluster_markers.out
}


