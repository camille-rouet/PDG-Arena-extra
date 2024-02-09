You will found here inventories, a launch script, original data, logs and analysis script related to the simulations described in Rouet et al. (2024): PDG-Arena: An eco-physiological model for characterizing tree-tree interactions in heterogeneous stands (submitted to PCI Forest and Wood Sciences).


# Description of the files

- 1_inventories_files
	- 2023-12-15_GMAP_hauteurMoy_goodAlignment
	>>	Contains all inventories (.inv files) used for PDG-Arena and CASTANEA simulations

		- 0_plots_image
		>>	Contains images of the simulated inventories

- 2_simulation_logs
	- multiScripts_copy.txt
	>>	Bash script used to run PDG-Arena and CASTANEA using Capsis in script mode

	- simulations-meta-informations.txt
	>>	Meta-informations about the simulation (elapsed time, capsis version..)

- 3_analysis_script
	- PDG-Arena_analyse_publication.R
	>>	The main script of the publication

	- CRMethodsForSimulations.R
	>>	A set of methods used in the main script

	- 2023-12-15_goodHeightGoodAlignment.RData
	>>	Contains the simulation data that are necessary to run the main script