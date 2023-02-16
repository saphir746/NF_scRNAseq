process make_scv_file {
  
  module "Anaconda3/2019.07"
  conda "/camp/stp/babs/working/schneid/conda/envs/R-4.2-Seurat"
  
  cpus 1
  time "6h"
  memory "100G"
    
  input:
    path(Rds)
  
  output:
    path("*.scv")
  
  script:
    """
    scviewer_makeFiles.R $Rds
    """
}
