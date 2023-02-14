process integration {
  
  module "Anaconda3/2019.07"
  conda "/camp/stp/babs/working/schneid/conda/envs/R-4.2-Seurat"
  
  cpus 1
  time "4h"
  memory "100G"
    
  input:
    path(Rds)
  
  output:
    path("*.RDS")
  
  script:
    """
    Rscript /camp/stp/babs/working/schneid/projects/sahaie/giovanni.giangreco/Characterisation_of_CAF_in_HPV_cancer_scrnaseq/Integration_general.R $Rds
    """
}
