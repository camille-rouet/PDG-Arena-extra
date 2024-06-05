Scripts and data associated with the PDG-Arena forest growth model. PDG-Arena is based on ecophysiological processes and tree-tree interactions. It is part of Physio-Demo-Genetics, a model that is maintained on the [Capsis platform](https://capsis.cirad.fr/capsis/help_en/physiodemogenetics).

Here, you will found :
- scripts that serves at inventory generation and simulation analysis;
- simulation inventories related to publications;
- simulation results files and logs related to publications.

## How to use 
- Git clone this repository or download and extract the repository ZIP file
- Open PDG-Arena-extra.Rproj using RStudio

## Publications

### Publication 1
Rouet et al. (2024): PDG-Arena: An eco-physiological model for characterizing tree-tree interactions in heterogeneous stands (preprint submitted to PCI Forest and Wood Sciences). [https://doi.org/10.1101/2024.02.09.579667](https://doi.org/10.1101/2024.02.09.579667)

Here are the files used in this publication:

> inventories/2023-12-15_GMAP_hauteurMoy_NoSlope_goodAlignment/ 

— This folder contains the virtual inventories used in the simulation

> inventories/2023-12-15_GMAP_hauteurMoygoodAlignment/sites_description/LAI_GMAP.csv

— Data table describing the simulated plots and in particular their Leaf Area Index

> inventories/2023-12-15_GMAP_hauteurMoy_goodAlignment/sites_description/soil_CASTANEA_GMAP.csv

— Data table describing the soil texture of GMAP sites

> scripts/generation-of-inventories/GMAPtoPDG.R 

— The script used to generates inventories

> scripts/analysis-of-simulations/PDG-Arena_publication_simulation_analysis.R 

— The script used to analyse the simulations and plot the result figures

> simulation-data/paper-PDG-Arena/2023-12-15_goodHeightGoodAlignment.RData 

— A Rdata file containing the simulation results and the dataset used in the publication (see the accompanying readme)

> simulation-data/paper-PDG-Arena/multiScripts_copy.txt 

— The script used to launch the simulation

> simulation-data/paper-PDG-Arena/simulations-meta-informations.txt 

— It contains simulation meta-information such as software version

> simulation-data/paper-PDG-Arena/LAI/

— It contains the values of LAI used in the simulations



## License
This work © 2024 by Camille Rouet is licensed under [CC BY-NC 4.0](http://creativecommons.org/licenses/by-nc/4.0/).
