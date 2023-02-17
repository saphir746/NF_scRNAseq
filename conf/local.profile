/*
 * configuration for local execution
 */

profiles {

	local {

		includeConfig "singularity.config"
		
		process.executor = "local"
	}

	cluster {

		includeConfig "singularity.config"
		process.executor = "slurm"
	}
}
