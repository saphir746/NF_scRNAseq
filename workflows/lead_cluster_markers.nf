#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2

//////////////////////////////////////////////////////
//// Made on my own (schneid@crick.ac.uk)  ///////////
//////////////////////////////////////////////////////


process lead_cluster_markers {
  
  label 'save_output'

  module params.MD_ANACONDA
  conda params.CONDA_ENV

  cpus 1
  time "6h"
  memory "100G"
  
  input:
    tuple path(Rds), val(res)
  
  output:
    path("*.csv")
  
  script:
    """
    Find_lead_markers_clusters.R $Rds $res
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


