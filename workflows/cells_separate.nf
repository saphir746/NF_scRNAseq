#!/usr/bin/env nextflow

import java.nio.file.Paths

nextflow.enable.dsl=2


include { integration } from '../modules/integration'
include { make_scv_file } from '../modules/make_scv_files'
include { LEAD_MARKERS } from '../workflows/lead_cluster_markers'

Channel
  .from( "'0.3'", "'0.9'" )
  .set{ RESOL }

//////////////////// processes/////////////////////

process split_cells {
  
  label 'save_output'
  
  module params.MD_ANACONDA
  conda params.CONDA_ENV
  
  cpus 1
  time "2h"
  memory "100G"

  input:
      path(RDSobj)

  output:
      path("*.RDS")

  script:
  """
  Cells_split.R $RDSobj
  """

}


/////////////////////////////////////////////////

workflow CELLS_SPLIT {
  take:
    resol
    rds
    
  main:
    split_cells(rds)
    split_cells.out | flatten | integration | make_scv_file
    integration.out.toList().set{ split_rds }

    LEAD_MARKERS(
      RESOL,
      integration.out
    )
 
  //  integration.out.combine(resol) | scCatch 
    //integration.out.combine(lead_cluster_markers.out).toList().view() 
  
  emit:
    split_rds
    LEAD_MARKERS.out
    
}
      
