process integration {
 
  label 'save_output'

  module params.MD_ANACONDA
  conda params.CONDA_ENV

  cpus 1
  time "4h"
  memory "100G"
    
  input:
    path(Rds)
  
  output:
    path("*.RDS")
  
  script:
    """
    Integration_general.R $Rds
    """
}
