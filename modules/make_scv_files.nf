process make_scv_file {
  
  module params.MD_ANACONDA
  conda params.CONDA_ENV
  
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
