process assign_identities {

  label 'save_output'

  module params.MD_ANACONDA
  conda params.CONDA_ENV

  input:
        	tuple path(Rds), path(annot)

        output:
                path("*_augmented.RDS")

        script:
        """
        Assign_cellID_custom.R $Rds $annot
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

