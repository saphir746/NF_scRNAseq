process assign_identities {

  label 'save_output'

  module params.MD_ANACONDA
  conda params.CONDA_ENV

  input:
        path(Rds)

        output:
                path("*_augmented.RDS")

        script:
        """
        Assign_cell_identities_custom.R $Rds    
        """
}

process transfer_identities {

  label 'save_output'

  module params.MD_ANACONDA
  conda params.CONDA_ENV

  input:
        path(Rds)

  output:
        path("*_augmented.RDS")

  script:
  """
  Transfer_cell_identities.R $Rds
  """
}

