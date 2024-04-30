# This work © 2024 by Camille Rouet is licensed under [CC BY-NC 4.0](http://creativecommons.org/licenses/by-nc/4.0/).
# Genetics computation segments were written by Sylvie Oddou-Muratorio.

# INIT ----
rm(list = ls()) ; gc()

library(tidyverse)
library(openxlsx)
library(data.table)

# README ----
"
HOW TO USE
Source this script
"


source("scripts/define_folders.R")
source("scripts/generation-of-inventories/CR_methods_for_inventory_generation.R")
source("scripts/analysis-of-simulations/CRMethodsForSimulations.R")


# PARAMETERS ----


outputInConsole = FALSE # if TRUE, output is redirected in console instead of in files
# outputInConsole = TRUE # to comment

# folder for writing inventories
outputInventoryFolder = "inventories/2024-04-30_GMAP_publication/"

# standard = trees are located and defined based on the GMAP inventory, output is one (eventually mixed) plot
# regdemo_monosp = trees are all of the same (mean) age and (quadratic mean) dbh and they are located on a regular basis, two monospecific plot are made
# regdemo_plurisp = trees of the same species are all of the same (mean) age and (quadratic mean) dbh and trees are located on a regular (randomized ?) basis, one plurispecific plot
# irregdemo_monosp = trees age and dbh are conserved based on inventories, trees are located on a regular (randomized ?) basis, two monospecific plot are made 
# allbutE1B = all modes but irregdemoMonosp
# all = all modes are written
inventoryMode = "allbutE1B" 

# if TRUE, LAI will be adjusted for each species (based on forrester allometry table). # If FALSE, LAI will be the same in regularized monospecific inventories and demography-fine plurispecific inventoriy
useSpeciesLAIduringRegularization = TRUE

# Put stand to zero slope
forceNoSlope = FALSE

# nTree per stand for regdemo mode
nTree_regularization = 36

nbAltitudeClasses = 5

# selection of sites
sites = c("vl", "vtx", "bg")

# first year of simulation / year before the simulation begins (=the climatic year) / annee avant le debut des simulation (=année climatique) (cr-15.03.2024?)
initYear = 1995

# If TRUE, diameter at a given are computed based on core data
useTimeMachine = TRUE

thisScriptOptions = tibble(outputInConsole = outputInConsole,
                           outputInventoryFolder = outputInventoryFolder,
                           inventoryMode = inventoryMode, 
                           useSpeciesLAIduringRegularization = useSpeciesLAIduringRegularization, 
                           forceNoSlope = forceNoSlope,
                           nTree_regularization = nTree_regularization,
                           nbAltitudeClasses = nbAltitudeClasses,
                           initYear = initYear,
                           useTimeMachine = useTimeMachine)



# Create the output inventory folder
if(!dir.exists(outputInventoryFolder)){
  dir.create(outputInventoryFolder)
}


# IMPORT ----

# GMAP sites and trees
pathFolder = paste0(myFolder, "01_docs-these/O_Sites_GMAP/GMAP_donnees/arbres_indiv_horsprod/")
fileName = "Data_arbres_placettes_horsprod.csv"
path = paste0(pathFolder, fileName)
dGMAP_horsprod_all = read_delim(path, delim = "\t", show_col_types = FALSE)

# Clean dGMAP_horsprod_all (dead trees not properly marked)
dGMAP_horsprod_all[grepl(dGMAP_horsprod_all$num_arbre, pattern = "mort") | grepl(dGMAP_horsprod_all$num_tige, pattern = "mort"), ]$mort = 1


# GMAP dendro : age and growth data
pathFolder = paste0(myFolder, "01_docs-these/O_Sites_GMAP/GMAP_donnees/croissance_matt_dendro/")
fileName = "matt_dendro_densite_croiss_interdat_age_tot.csv"
path = paste0(pathFolder, fileName)
dGMAP_dendro = read_delim(path, delim = ";", show_col_types = FALSE)


# GMAP soil : ventoux and vercors
pathFolder = paste0(myFolder, "01_docs-these/O_Sites_GMAP/sol/")
fileName = "soil_CASTANEA_GMAP_depth1600.csv"
path = paste0(pathFolder, fileName)
dGMAP_soil = read_delim(path, delim = ";", show_col_types = FALSE)


# GMAP LIDAR
pathFolder = paste0(myFolder, "01_docs-these/N_LIDAR/__gitlab/donnees/")
fileName = "GMAP_sites.csv"
path = paste0(pathFolder, fileName)
dGMAP_LIDAR = read_delim(path, delim = ";", show_col_types = FALSE)


# # CASTANEA SPECIES FILE
CASTANEASpeciesFile = readCASTANEASpeciesFile(pathtofile = paste0(capsisPath, "data/castaneaonly/species/CastaneaSpecies_2022_05.txt"))
speciesParameters = CASTANEASpeciesFile$parameters
CASTANEASpeciesTable = CASTANEASpeciesFile$speciesTable





# *****************************************
# 1 = LOADING Parameters [you can change here their value] / chargement des parametres [vous pouvez changer leur valeur ici] ----

# *****************************************

# # Inventory option file
pdgOptionFileFolder = paste0(myFolder, "01_docs-these/O_Sites_GMAP/GMAP_donnees/")
pdgOptionFileName = "PDG_inventory_options_GMAP_simulations.txt"
pdgInventoryOptions = readInventoryOptionFile(filePath = paste0(pdgOptionFileFolder, pdgOptionFileName))

# Can use this :
# pdgInventoryOptions = readInventoryOptionFile("scripts/generation-of-inventories/PDG_inventory_options.txt")



# MODIFY options

# demographic parameters are considered equal between species
pdgInventoryOptions$demographic_parameters[2, c(is.na(pdgInventoryOptions$demographic_parameters[2, ])) ] = pdgInventoryOptions$demographic_parameters[1,  c(is.na(pdgInventoryOptions$demographic_parameters[2, ])) ]

pdgInventoryOptions$PDG_simulation_parameters$initYear = initYear


# IMPORT options

demographic_parameters = pdgInventoryOptions$demographic_parameters
PDG_simulation_parameters = pdgInventoryOptions$PDG_simulation_parameters
CASTANEA_simulation_parameters = pdgInventoryOptions$CASTANEA_simulation_parameters
genetics_parameters = pdgInventoryOptions$genetics_parameters


speciesFullNames = demographic_parameters$speciesFullNames
speciesFrenchNames = demographic_parameters$speciesFrenchNames
fmIdSp = demographic_parameters$fmIdSp # castanea species codes
idSp = demographic_parameters$idSp # species name for PDG inventory

meanPAI_LIDAR = round(mean(dGMAP_LIDAR$PAI, na.rm = T), digits = 2)



cellWidth = PDG_simulation_parameters$targetCellWidth





# LOADING GENETIC PARAMETERS / PARAMETRES GENETIQUES

# FOR NOW, GENETICS VARIATION IN MULTISPECIES IS NOT AVAILABLE
# (to implement is, one should make a table including genetics parameters per species + implement the modification during file importation in PDG in capsis)
# THUS, GENETIC PARAMETERS ARE EQUAL FOR ALL SPECIES
# If disableGenetics if TRUE, value of FCRITBB AND G1VAR will be found in CASTANEA species file
disableGenetics = PDG_simulation_parameters$disableGenetics

# genetic parameters uniques
gen_param_u = genetics_parameters[1,]


# Import genetics parameters
numberOfGeneticParameters = gen_param_u$numberOfGeneticParameters
numberOfTraits = gen_param_u$numberOfTraits

# genetics parameters FCRITBB
nbLocus_FCRITBB=	gen_param_u$nbLocus_FCRITBB
nbAllelePerLoc_FCRITBB= gen_param_u$nbAllelePerLoc_FCRITBB  ##### BE CAREFUL : can be only 2
interStepEnvirVar_FCRITBB= gen_param_u$interStepEnvirVar_FCRITBB ##### BE CAREFUL : totEnvirVar_FCRITBB can be only 0 
herit_FCRITBB=	gen_param_u$herit_FCRITBB ##### BE CAREFUL : herit_FCRITBB can be either 1 or 0; if herit_FCRITBB = 1, totEnvirVar_FCRITBB will be considered as null

target_mean_FCRITBB= gen_param_u$target_mean_FCRITBB 
target_CV_FCRITBB=	gen_param_u$target_CV_FCRITBB
alleleEffectMultiplCoeffFCRITBB= gen_param_u$alleleEffectMultiplCoeffFCRITBB 

# genetics parameters g1max
nbLocus_g1max= gen_param_u$nbLocus_g1max
nbAllelePerLoc_g1max= gen_param_u$nbAllelePerLoc_g1max ##### BE CAREFUL : can be only 2
interStepEnvirVar_g1max= gen_param_u$interStepEnvirVar_g1max  ##### BE CAREFUL : interStepEnvirVar_g1max can be only 0 
herit_g1max=	gen_param_u$herit_g1max  ##### BE CAREFUL : herit_g1max can be either 1 or 0; if herit_g1max = 1, totEnvirVar_g1max will be considered as null

target_mean_g1max= gen_param_u$target_mean_g1max     
target_CV_g1max=	gen_param_u$target_CV_g1max     
alleleEffectMultiplCoeffg1max= gen_param_u$alleleEffectMultiplCoeffg1max 

nbMsat= gen_param_u$nbMsat
nbAllPerMsat= gen_param_u$nbAllPerMsat
nbSNP= gen_param_u$nbSNP
nbAllPerSNP= gen_param_u$nbAllPerSNP
distanceClassBoundsForSGS = pdgInventoryOptions$distance_class_bounds_for_SGS


# Compute derivated genetics parameters
targetSd_FCRITBB=target_CV_FCRITBB*target_mean_FCRITBB #target sd in FCRITBB
targetVariance_FCRITBB=targetSd_FCRITBB^2

targetSd_g1max=target_CV_g1max*target_mean_g1max #target sd in g1max
targetVariance_g1max=targetSd_g1max^2


# end genetic parameters


# DATA TREATMENT ----

# SITE SELECTION
dGMAP_horsprod = subset(dGMAP_horsprod_all, site %in% sites)
dGMAP_dendro = subset(dGMAP_dendro, site %in% sites)
dGMAP_LIDAR = subset(dGMAP_LIDAR, site %in% sites)
dGMAP_soil = dGMAP_soil[vapply(dGMAP_soil$code_site, function(x) grepl(x, pattern = "vl_|vtx_|bg_"), FUN.VALUE = T), ]


# CLEAN GMAP TABLES
dGMAP_horsprod_save_mort = dGMAP_horsprod # save all tree, including dead ones
# Keep alive trees only
dGMAP_horsprod = subset(dGMAP_horsprod,  mort == 0) 
dGMAP_horsprod_save_vivant_allspecies = dGMAP_horsprod
# Keep species of interest only
dGMAP_horsprod = subset(dGMAP_horsprod,  essence %in% speciesFrenchNames)

# num arbre as factor
dGMAP_horsprod$num_arbre = as.numeric(dGMAP_horsprod$num_arbre)

dGMAP_dendro$num_index = getNumIndex(dGMAP_dendro$num_arbre, dGMAP_dendro$num_tige)
dGMAP_dendro$treeGlobalId = paste0(dGMAP_dendro$code_site, '_', dGMAP_dendro$num_index)
dGMAP_horsprod$num_index = getNumIndex(as.numeric(dGMAP_horsprod$num_arbre), dGMAP_horsprod$num_tige)
dGMAP_horsprod$treeGlobalId = paste0(dGMAP_horsprod$code_site, '_', dGMAP_horsprod$num_index)


# GMAP tree-year table
treeYearTable_GMAP = makeGMAP_TreeYearTable(dGMAP_dendro, dGMAP_horsprod, years = initYear:2013)



# *****************************************
# 2 = GMAP Data / Données GMAP ----
# *****************************************

all_code_site = unique(dGMAP_horsprod$code_site)

# create a table with all possible inventories
all_inventoryType = c("irregdemo_plurisp", "regdemo_monosp_hetre", "regdemo_monosp_sapin", "CASTANEA_hetre", "CASTANEA_sapin",
                      "regdemo_plurisp", "irregdemo_monosp_hetre", "irregdemo_monosp_sapin")
inventoryTable = tibble(code_site = rep(all_code_site, each = length(all_inventoryType)), inventoryType = rep(all_inventoryType, times = length(all_code_site)))
inventoryTable$bee_BAprop_init = -1
inventoryTable$fir_BAprop_init = -1
inventoryTable$comp = ""
inventoryTable$created = FALSE
inventoryTable$inventoryName = ""
inventoryTable$bee_nprop = -1
inventoryTable$fir_nprop = -1
inventoryTable$bee_BAprop = -1 # BA_beech / (BA_beech + BA_fir)
inventoryTable$fir_BAprop = -1
inventoryTable$GHA_init = -1 # BA / standArea (cm3 / m2)
inventoryTable$GHA = -1

inventoryTable$bee_dbh = -1 # dbh quadratique moyen
inventoryTable$fir_dbh = -1
inventoryTable$bee_dbh_init = -1
inventoryTable$fir_dbh_init = -1

code_site = all_code_site[11] # for line-per-line testing purpose, no impact

# 2 LOOP ON CODESITE ----
# for(code_site in c("vl_inter_m_2")){
for(code_site in all_code_site){
  
  
  if(grepl(code_site, pattern = "old")){
    inventoryTable = inventoryTable[inventoryTable$code_site != code_site, ]
    next
  }
  
  if( !code_site %in% dGMAP_LIDAR$code_site & code_site %in% dGMAP_LIDAR$code_site_2 ){
    stop(paste0( code_site, "is in dGMAP_LIDAR$code_site2"))
  }
  
  
  
  # 2-a SELECT SITE AND TREES ----
  cat("\nReading: ", code_site, " (", which(all_code_site == code_site), "/", length(all_code_site), ")\n", sep = "")
  
  # Selection du site
  
  # horsprod
  dGMAP_hp_site = dGMAP_horsprod[dGMAP_horsprod$code_site == code_site, ]
  dGMAP_hp_site_save_mort =  dGMAP_horsprod_save_mort[dGMAP_horsprod_save_mort$code_site == code_site, ]
  dGMAP_hp_site_save_vivant_allspecies = dGMAP_horsprod_save_vivant_allspecies[dGMAP_horsprod_save_vivant_allspecies$code_site == code_site, ]
  
  # dendro
  if(code_site %in% dGMAP_dendro$code_site){
    dGMAP_dendro_site = dGMAP_dendro[dGMAP_dendro$code_site == code_site, ]
  }else{
    code_site_bis = substr(code_site, 1, nchar(code_site)-4)
    if(code_site_bis %in% dGMAP_dendro$code_site){
      dGMAP_dendro_site = dGMAP_dendro[dGMAP_dendro$code_site == code_site_bis, ]
    }else{
      # redirect output to console
      while(sink.number() > 0){sink()}
      cat("No dendro data found\n")
      warning(paste0("Stand ", code_site, " dont have dendro data"))
      inventoryTable = inventoryTable[inventoryTable$code_site != code_site, ]
      next
    }
  }
  
  # LIDAR
  LAI_LIDAR = meanPAI_LIDAR
  if(code_site %in% dGMAP_LIDAR$code_site){
    dGMAP_LIDAR_site = dGMAP_LIDAR[dGMAP_LIDAR$code_site == code_site, ]
    if(dGMAP_LIDAR_site$hasLIDAR){
      LAI_LIDAR = dGMAP_LIDAR_site$PAI
    }else{
      warning(paste0("Site ", code_site, " was not LIDARized."))
    }
  }else{
    warning(paste0("Site ", code_site, " is not in LIDAR data."))
  }
  
  
  globalSite = unique(dGMAP_hp_site$site)
  subSite = unique(dGMAP_hp_site$sous_site)
  standType = unique(dGMAP_hp_site$peuplement)
  
  
  
  nTiges = dim(dGMAP_hp_site)[1]
  
  
  # 2-b SITE PARAMETERS ----
  
  referenceAltitude = unique(dGMAP_hp_site$altitude)   # Altitude at which x = y = 0 (or, simply, y = 0 when aspect is 0)
  longitude= unique(dGMAP_hp_site$position_e)     # longitude in degree for CASTANEA library / longitude en degré pour la librairie CASTANEA
  latitude= unique(dGMAP_hp_site$position_n)     # latitude in degree for CASTANEA library / latitude en degré pour la librairie CASTANEA
  
  # if slope_deg = 45, then 1m along y axis correpond to 1m a=of altitudinal variation
  slope_deg = unique(dGMAP_hp_site$pente) # Slope of the plot in degree / pente en degré
  
  aspect_deg = as.numeric(unique(dGMAP_hp_site$exposition)) # exposition of the slope in degree
  
  if(forceNoSlope){
    slope_deg = 0
    aspect_deg = 0
  }

  
  pdgInventoryOptions$plotLocationParameters = tibble(referenceAltitude = referenceAltitude, latitude = latitude, longitude = longitude,
                                                      slope_deg = slope_deg, aspect_deg = aspect_deg)
  
  # evaluate distance to center
  minDistCenterX = 1e9
  minDistCenterY = 1e9
  maxDistCenterX = -1e9
  maxDistCenterY = -1e9
  
  # find maximum and minimum coordinates
  for(i in 1:nTiges){
    dataTree = dGMAP_hp_site[i, ]
    
    # L'angle de position des arbres est donné dans le sens des aiguilles d'une montre (sens anti-trigonométrique)
    # et il est mesuré par rapport au Nord (axe y)
    # nous voulons l'angle dans le sens trigonométrique et partant de l'axe x :
    angle_deg = 90 - dataTree$angle
    angle_rad = angle_deg / 360 * 2 * pi
    
    # la distance au centre est vue d'en haut 
    distCenter = dataTree$distance_centre
    
    treeDistCenterX = distCenter * cos(angle_rad)
    treeDistCenterY = distCenter * sin(angle_rad)
    
    minDistCenterX = min(minDistCenterX, treeDistCenterX)
    maxDistCenterX = max(maxDistCenterX, treeDistCenterX)
    minDistCenterY = min(minDistCenterY, treeDistCenterY)
    maxDistCenterY = max(maxDistCenterY, treeDistCenterY)
  }
  rm(i)
  
  # min and maximum plot positions
  xorigin	= 0
  yorigin	= 0
  xMaxTree = max( abs(maxDistCenterX), abs(minDistCenterX) ) * 2 + 0.01
  yMaxTree = max( abs(maxDistCenterY), abs(minDistCenterY) ) * 2 + 0.01
  # xMaxTree = maxDistCenterX - minDistCenterX
  # yMaxTree = maxDistCenterY - minDistCenterY
  
  
  nbColumns	= ceiling(xMaxTree / cellWidth)       # number of cells along x-axis / nombre de cellules lel long de l'axe x
  nbLines	= ceiling(yMaxTree / cellWidth)        # number of cells along y-axis / nombre de cellules lel long de l'axe y
  
  nbCell = nbLines*nbColumns  # number of cells
  
  xMin = 0; xMax = xMin + cellWidth * nbColumns ; yMin = 0; yMax = yMin + cellWidth * nbLines
  
  # center of the study field position
  centerX = xMax / 2
  centerY = yMax / 2
  # centerX = - minDistCenterX
  # centerY = - minDistCenterY
  
  # # # altitude consideration 
  
  
  altitudeClasses = getAltittudeClases(xMin, xMax, yMin, yMax, slope_deg, aspect_deg, referenceAltitude, nbAltitudeClasses)
  
  
  
  # 2-b SOIL PARAMETERS ----
  code_site_soil = paste(globalSite, subSite, standType, sep = '_')
  
  dGMAP_soil_site = subset(dGMAP_soil, code_site == code_site_soil)
  
  # potentially variable among cells but  at this stage only for soilHeight
  soilHeight = dGMAP_soil_site$soilHeight         # soil Height in mm / hauteur de sol en mm
  VariabSoilHeith = 0
  
  stoneContent = dGMAP_soil_site$stoneContent      #percentage of stone in the soil / % de pierre dans le sol
  wfc = dGMAP_soil_site$wfc 
  wilt = dGMAP_soil_site$wilt
  propMacro = dGMAP_soil_site$propMacro
  propMacroDeep = dGMAP_soil_site$propMacroDeep
  bulk = dGMAP_soil_site$bulk
  SOLCLAYtop = dGMAP_soil_site$SOLCLAYtop   # proportion d'argile dans le sol (0 - 30 cm)
  SOLCLAYsol = dGMAP_soil_site$SOLCLAYsol	 # proportion d'argile dans le sol (30 - 100 cm)
  SOLFINtop = dGMAP_soil_site$SOLFINtop # proportion de particules fines (argiles + limons) (0 - 30 cm)
  SOLFINsol = dGMAP_soil_site$SOLFINsol # proportion de particules fines (argiles + limons) (30 - 100 cm)
  SOLSANDtop = dGMAP_soil_site$SOLSANDtop # proportion de sables (0 - 30 cm)
  SOLSANDsol = dGMAP_soil_site$SOLSANDsol # proportion de sables (30 - 100 cm)
  deepSoilDepth = dGMAP_soil_site$deepSoilDepth #depth of deep soil (especially for karstic zones)
  stoneContentDeep = dGMAP_soil_site$stoneContentDeep # rate of stone content in deep soil
  prac = dGMAP_soil_site$prac
  pracDeep = dGMAP_soil_site$pracDeep
  # end soil parameters
  
  
  
  # 2-b NA correction ----
  # dbh, height, hcb, age
  
  # Correction of dbh / circonference for NA
  circonferenceIsNa = which( is.na(dGMAP_hp_site$circonference) )
  
  if( length(circonferenceIsNa) > 0 ) {
    dGMAP_validCirconference = dGMAP_hp_site[-circonferenceIsNa, ]
    meanCirconference = mean(dGMAP_validCirconference$circonference)
    
    dGMAP_hp_site[circonferenceIsNa, ]$circonference = round( meanCirconference, 3)
  }
  
  
  # Correction of height for NA
  heightIsNa = is.na(dGMAP_hp_site$htot)
  
  if( length(which(heightIsNa)) > 0 ) {
    # linear model by species to estimate missing values
    dGMAP_hp_site = computeMissingVariableProportionnaly(dGMAP_hp_site, "dbh", "htot", "essence", heightIsNa)
    
    # old way
    # dGMAP_hp_site$dbh = dGMAP_hp_site$circonference / pi
    # dGMAP_validHeight = dGMAP_hp_site[-heightIsNa, ]
    # 
    # 
    # heightRatio_1 = computeBasalAreaRatio(dGMAP_validHeight, "dbh", "htot", species = "hetre", variableSpecies = "essence")
    # heightRatio_2 = computeBasalAreaRatio(dGMAP_validHeight, "dbh", "htot", species = "sapin", variableSpecies = "essence")
    # estimatedHeight = dGMAP_hp_site$dbh * ifelse(dGMAP_hp_site$essence == "hetre", heightRatio_1, heightRatio_2)
    # meanHeight =  mean(dGMAP_validHeight$htot)
    # 
    # dGMAP_hp_site[heightIsNa, ]$htot = round( meanHeight , 3)
  }
  
  
  # Correction of hcb for NA
  crownSpanIsTooLow = dGMAP_hp_site$htot - dGMAP_hp_site$h_base_houppier < 0.5
  hcbIsNa = is.na(dGMAP_hp_site$h_base_houppier)
  hcbToCorrect = crownSpanIsTooLow | hcbIsNa
  
  if( length(which(hcbToCorrect)) > 0 ) {
    # linear model by species to estimate missing values
    dGMAP_hp_site = computeMissingVariableProportionnaly(dGMAP_hp_site, "htot", "h_base_houppier", "essence", hcbToCorrect)
    
    # old way
    # dGMAP_validHcb = dGMAP_hp_site[-hcbIsNa, ]
    # 
    # hcbRatios = dGMAP_validHcb$h_base_houppier / dGMAP_validHcb$htot
    # meanHcbRatio = mean( hcbRatios )
    # 
    # dGMAP_hp_site[hcbIsNa, ]$h_base_houppier = round( dGMAP_hp_site$htot[hcbIsNa] * meanHcbRatio , 3)
  }
  
  # Correction of age for NA (HERE, FOR TREES THAT ARE IN GMAP_dendro /!\)
  ageIsNa = which(is.na(dGMAP_dendro_site$age))
  
  if(length(ageIsNa) > 0){
    
    # 1. (a) extract the lines of Matt_dendro for the case where several lines are represented for the same stem
    # (b) If so, age can be NA for a line and defined for another.
    # (c) Then, every lines are filled with the same age
    for(index in ageIsNa){ # for each NA as age
      # (a)
      num_tige_index = dGMAP_dendro_site[index, ]$num_tige
      lines = subset(dGMAP_dendro_site, num_tige == num_tige_index)
      if(dim(lines)[1] > 1){
        ageNotNaInLines = !is.na(lines$age)
        if(sum(ageNotNaInLines) > 0){ # (b)
          dGMAP_dendro_site[index, ]$age = unique(lines[ageNotNaInLines,]$age) # (c)
        }
      }
    } # end for loop
  }
  
  ageIsNa_table = is.na(dGMAP_dendro_site$age)
  ageIsNa = which(ageIsNa_table)
  
  if(length(ageIsNa) > 0){
    # 2. Now, the trees whose age could not be corrected by several lines in matt dendro.
    
    # linear model by species to estimate missing values
    dGMAP_dendro_site = computeMissingVariableProportionnaly(dGMAP_dendro_site, "perimetre", "age", "essence", ageIsNa_table)
    dGMAP_dendro_site$age = round(dGMAP_dendro_site$age)
    
    if(TRUE %in% (dGMAP_dendro_site$age <= 1 + 2013 - initYear)){
      # arbitrary : trees can not be under 2 at initYear
      dGMAP_dendro_site$age = ifelse(dGMAP_dendro_site$age < 2013 - initYear + 2, 2013 - initYear + 2, dGMAP_dendro_site$age)
    }
    
    
    # old way
    # Their age is the mean age of the plot
    # dGMAP_validAge =  dGMAP_dendro_site[ -ageIsNa, ]
    # meanAge = mean(dGMAP_validAge$age)
    # dGMAP_dendro_site[ageIsNa, ]$age = round(meanAge)
    

    # # # alternative : Their age is computed with a dbh-based ratio
    # # meanAgeRatio = mean(dGMAP_validAge$age/dGMAP_validAge$perimetre)
    # # dGMAP_dendro_site[ageIsNa, ]$age = round(dGMAP_dendro_site[ageIsNa, ]$perimetre * meanAgeRatio)
  }
  
  # 2-c TREES ATTRIBUTES ----
  
  # new table
  trees = tibble("numTigeChr" = dGMAP_hp_site$num_tige)
  
  trees$idTree = -1
  trees$pop = -1
  trees$sp = -1
  trees$cell_IDs = -1
  trees$x_abs = -1
  trees$y_abs = -1
  trees$xnorm = -1
  trees$ynorm = -1
  trees$dbh = -1
  trees$dbh2013 = -1
  trees$height = -1
  trees$height2013 = -1
  trees$hcb = -1
  trees$hcb2013 = -1
  trees$age = -1
  trees$age2013 = -1
  trees$LAIinit = -1
  
  nbIndividuals = c(0,0)
  
  
  # LOOP ON STEMS
  for(numTigeChr in trees$numTigeChr){
    
    whichTree = which(trees$numTigeChr == numTigeChr)
    tree = trees[whichTree, ]
    
    # load data for this stem
    dataTree = dGMAP_hp_site[dGMAP_hp_site$num_tige == numTigeChr, ]
    
    treeYearTable_tree = subset(treeYearTable_GMAP, treeGlobalId == dataTree$treeGlobalId)
    treeYearTable_tree_initYear = subset(treeYearTable_tree, year == initYear)
    
    # stem index
    numArbre = dataTree$num_arbre
    numTige = dataTree$num_tige
    
    # load dendro data for this stem
    hasDendro = numTige %in% dGMAP_dendro_site$num_tige
    if(hasDendro){
      dataTree_dendro = dGMAP_dendro_site[dGMAP_dendro_site$num_tige == numTige, ]
    }
    
    # -> as.numeric("11e") gives numeric 11 ..
    # Then, isCoppice = suppressWarnings( is.na(as.numeric(numTige)) ) # with coppice, numTige as a letter
    # is replaced by : 
    numbers_only <- function(x) !grepl("\\D", x)
    isCoppice = !numbers_only(numTige)
    
    # define a numerical index for each stem
    if(isCoppice){
      tige_letter = substr(numTige, start = nchar(numArbre) + 1, stop = nchar(numTige))
      tige_letter_as_index = utf8ToInt(tige_letter) - utf8ToInt('a') + 1 # from letter to integer : index from 1 to XX
      num_index = numArbre + tige_letter_as_index * 1000
    }else{
      num_index = numArbre
    }
    
    trees$idTree[whichTree] = num_index
    
    # species and counting
    if(dataTree$essence == "hetre") {
      trees$sp[whichTree] = idSp[1]
      nbIndividuals[1] = nbIndividuals[1] + 1
    }else if(dataTree$essence == "sapin"){
      trees$sp[whichTree] = idSp[2]
      nbIndividuals[2] = nbIndividuals[2] + 1
    }
    
    trees$pop[whichTree] = 1
    
    # diameter (of 2013 and eventually older)
    dbh2013 = dataTree$circonference / pi
    height2013 = dataTree$htot
    hcb2013 = dataTree$h_base_houppier
    
    trees$dbh2013[whichTree] = dbh2013
    trees$height2013[whichTree] = height2013
    trees$hcb2013[whichTree] = hcb2013
    
    
    # cr-2024.04.30, this is replaced by treeYearTable for harmonizing between inventory generation and analysis
    # # use dbh from initYear computed using core data
    # timeMachineDBHCorrection = 1 # proportion dbhInitYear / dbh2013 
    # if (useTimeMachine) {
    #   
    #   if(hasDendro){
    #   
    #     coreIndexes = which(colnames(dataTree_dendro) == paste0("X", initYear)):which(colnames(dataTree_dendro) == "X2013")
    #     
    #     # if there is several lines in dataTree_dendro (it happens..), choose the line with max value of X2010 (arbitrary)
    #     max2010GrowthLine = which.max(dataTree_dendro[["X2010"]])
    #     dataTree_dendro1 = dataTree_dendro[max2010GrowthLine,]
    #     
    #     BAI = sum(dataTree_dendro1[, coreIndexes]) # cr-2024.04.26 should total BAI ignore the 2013 year, or the 1996 year ?
    #     dbhInitYear = sqrt(dbh2013 ** 2 - BAI * 4 / pi) # calcul de dbh init en cm a partir de dbh 2013 et du BAI en cm2 
    #     
    #     
    #     trees$dbh[whichTree] = round(dbhInitYear, 3)
    #     
    #     timeMachineDBHCorrection = dbhInitYear / dbh2013 # for later use with height
    #     timeMachineBasalAreaCorrection = dbhInitYear**2 / dbh2013**2 # for later use with height
    #     
    #   }else{ # no dendrometric data
    #     trees$dbh[whichTree] = -11
    #   }
    #   
    # } else{
    #   trees$dbh[whichTree] = round(dbh2013, 3) # $circonference pour Data_arbre_placettes*.csv / $perimetre pour matt_dendro*.csv
    # }
    # 
    # 
    # # height alteration by timeMachine (only if useTimeMachine == TRUE)
    # trees$height[whichTree] = round(dataTree$htot * timeMachineDBHCorrection, 3) # $htot pour Data_arbre_placettes*.csv / $hverticale pour matt_dendro*.csv
    # trees$hcb[whichTree] = round(dataTree$h_base_houppier * timeMachineDBHCorrection, 3) # $h_base_houppier pour Data_arbre_placettes*.csv / $hhouppier pour matt_dendro*.csv
    
    if (useTimeMachine) {
      trees$dbh[whichTree] = treeYearTable_tree_initYear$dbhRetro
      trees$height[whichTree] = treeYearTable_tree_initYear$hauteurRetro
      trees$hcb[whichTree] = treeYearTable_tree_initYear$hauteurRetro * (hcb2013 / height2013) # hcb is kept proportionnal to height
    }else{
      trees$dbh[whichTree] = dbh2013
      trees$height[whichTree] = height2013
      trees$hcb[whichTree] = hcb2013
    }
    
    
    
    # handle age in matt_dendro file (there can be two line for one stem..)
    if(hasDendro){
      anAge = dataTree_dendro$age
      
      # Some num_tige appears several times in the matt_dendro file
      # Then, we search value that is not NA
      nValuesDendro = length(dataTree_dendro$age)
      
      if(nValuesDendro > 1){ # if they are several lines in dendro data, find the not-NA one
        for(value in dataTree_dendro$age){
          anAge = value
          if(!is.na(anAge)){
            break
          }
        } # end loop
      }
      
      trees$age2013[whichTree] = anAge
      trees$age[whichTree] = anAge
      
      # tree age at the year initYear with time machine
      if (useTimeMachine) {
        trees$age[whichTree] = trees$age[whichTree] - (2013 - initYear)
        trees$age[whichTree] = max(trees$age[whichTree], 1)
      }
      
      
    }else{ # if no dendro data
      trees$age[whichTree] = -1
    }
    
    
    trees$LAIinit[whichTree] = 0
    
    
    # 2-d TREES POSITION ----
    
    # L'angle de position des arbres est donné dans le sens des aiguilles d'une montre (sens anti-trigonométrique)
    # et il est mesuré par rapport au Nord (axe y)
    # nous voulons l'angle dans le sens trigonométrique et partant de l'axe x :
    angle_deg = 90 - dataTree$angle
    angle_rad = angle_deg / 360 * 2 * pi
    
    # la distance au centre est vue d'en haut 
    distCenter = dataTree$distance_centre
    
    treeDistCenterX = distCenter * cos(angle_rad)
    treeDistCenterY = distCenter * sin(angle_rad)
    
    treeX = centerX + treeDistCenterX
    treeY = centerY + treeDistCenterY
    
    treeCellColumn = floor(treeX / cellWidth)
    treeCellLineUpWard = floor(treeY / cellWidth)
    treeCellLine = (nbLines - 1) - treeCellLineUpWard # y axe and cell vertical index are opposite 
    
    
    # coordinate relative to the tree cell origin
    treeRelativeX = treeX - treeCellColumn * cellWidth
    treeRelativeY = treeY - treeCellLineUpWard * cellWidth
    
    treeCell_ID = treeCellLine * nbColumns + treeCellColumn + 1
    
    trees$x_abs[whichTree] = treeX
    trees$y_abs[whichTree] = treeY
    trees$cell_IDs[whichTree] = treeCell_ID
    trees$xnorm[whichTree] = round(treeRelativeX, 3)
    trees$ynorm[whichTree] = round(treeRelativeY, 3)
    
    if(treeCell_ID < 0){
      browser()
    }
    
  } # end loop on trees
  
  
  
  # NON-CORED TREES, age, dbh.. ----
  # FOR TREES THAT ARE NOT IN GMAP_dendro /!\
  
  # Correction of age
  trees$ageInDendro = trees$age >= 0
  ageNotInDendro = which(!trees$ageInDendro)
  if(length(ageNotInDendro) > 0){
    
    # linear model
    # "pour les arbres dont on connait l'age, pour un dbh de X, on a un age de Y"
    trees = computeMissingVariableProportionnaly(trees, "dbh2013", "age2013", "sp", variableYisMissing = !trees$ageInDendro)
    trees$age2013 = round(trees$age2013)
    trees$age = trees$age2013 - (2013 - initYear)
    
    if(TRUE %in% (trees$age <= 1)){
      trees$age = ifelse(trees$age <= 2, 2, trees$age)
    }
    
    # old way
    # mean age
    # trees_validAge = trees[- ageNotInDendro, ]
    # meanAge = mean(trees_validAge$age)
    # trees[ageNotInDendro, ]$age = round(meanAge)
    
    # oldway, alternative dbh-based proportionnal age
    # trees_validAge = trees[- ageNotInDendro, ]
    # meanAgeRatio = mean(trees_validAge$age2013 /  trees_validAge$dbh2013)
    # trees[ageNotInDendro, ]$age2013 = round(meanAgeRatio * trees[ageNotInDendro, ]$dbh2013)
    # trees[ageNotInDendro, ]$age = trees[ageNotInDendro, ]$age2013 - (2013- initYear)
  }
 
  # check for unexpected value
  if(TRUE %in% trees$dbh == -1){
    stop("Some tree has -1 as dbh. Which means that they have no dbh attributed..")
  }
  
  # is there unexpected negative or null value ?
  if(TRUE %in% (trees$dbh <= 0)){
    stop("Some unexpected negative of null value of dbh has pop up..")
  }
  
  # After these two checks, only positive or -11 value of dbh can be found
  
  # cr-2024.04.30, this is replaced by using treeYearTable for harmonizing between inventory generation and analysis
  # # Correction of passed dbh, for trees whose basalAreaIncrement is unknown
  # if(useTimeMachine){
  #   
  #   dbhIsMissingList = trees$dbh == -11
  #   
  #   if(TRUE %in% dbhIsMissingList){
  #     
  #     # # restart debug
  #     # trees$dbh[trees$missingGrowthData] = -11
  #     # trees$height[trees$missingGrowthData] = trees$height2013[trees$missingGrowthData]
  #     # trees$hcb[trees$missingGrowthData] = trees$hcb2013[trees$missingGrowthData]
  #     # # end restart debug
  #     # 
  #     # # debug
  #     # trees = trees[trees$sp == 1, ]
  #     # trees$missingGrowthData = trees$dbh == -11
  #     # valid_trees = trees[!trees$missingGrowthData, ]
  #     # 
  #     # # proportionnal variation dbh
  #     # trees = computeMissingVariableProportionnaly(trees, "dbh2013", "dbh", "sp", variableYisMissing = trees$missingGrowthData, method = "proportionnal_mean")
  #     # trees = computeMissingVariableProportionnaly(trees, "height2013", "height", "sp", variableYisMissing = trees$missingGrowthData, method = "proportionnal_mean")
  #     # trees = computeMissingVariableProportionnaly(trees, "hcb2013", "hcb", "sp", variableYisMissing = trees$missingGrowthData, method = "proportionnal_mean")
  #     # 
  #     # # plot dbh vs dbh 2013
  #     # ggplot(trees, aes(x = dbh2013, y = dbh, color = as.factor(sp))) + geom_point(aes(size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), method = "lm") + coord_cartesian(ylim = c(-10, 80), xlim = c(0,90))
  #     # 
  #     # # plot height and hcb vs dbh 2013
  #     # color1 = hsv(0.7, 0.4, 1)
  #     # color2 = hsv(0.1, 0.7, 1)
  #     # ggplot(trees, aes(x = dbh2013, y = height, shape = as.factor(sp))) + geom_point(color = color1, aes(size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), method = "lm", color = color1) + coord_cartesian(ylim = c(0, 40), xlim = c(0,90)) +
  #     #   geom_point(color = color2, aes(x = dbh2013, y = hcb, size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), linetype = "dashed", color = color2, aes(x = dbh2013, y = hcb), method = "lm")
  #     # 
  #     # 
  #     # ggplot(trees, aes(x = height2013, y = height, shape = as.factor(sp))) + geom_point(color = color1, aes(size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), method = "lm", color = color1) +
  #     #   coord_cartesian(ylim = c(0, 35), xlim = c(0,40)) +
  #     #   geom_point(color = color2, aes(x = height2013, y = hcb, size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), linetype = "dashed", color = color2, aes(x = height2013, y = hcb), method = "lm") +
  #     #   geom_abline(slope = 1, alpha = 0.5)
  #     # 
  #     # ggplot(trees, aes(x = hcb2013, y = height, shape = as.factor(sp))) + geom_point(color = color1, aes(size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), method = "lm", color = color1) +
  #     #   coord_cartesian(ylim = c(0, 35), xlim = c(0,40)) +
  #     #   geom_point(color = color2, aes(x = hcb2013, y = hcb, size = missingGrowthData)) +
  #     #   geom_smooth(data = subset(trees, !missingGrowthData), linetype = "dashed", color = color2, aes(x = hcb2013, y = hcb), method = "lm") +
  #     #   geom_abline(slope = 1, alpha = 0.5)
  #     # 
  #     # 
  #     # trees$height - trees$hcb
  #     # 
  #     # # end debug
  #     # 
  #     # # test non negative coefficients
  #     # 
  #     # ggplot(trees, aes(x = dbh2013, y = hcb2013, color = as.factor(sp))) + geom_point() + geom_smooth(method = "lm", data = subset(trees, !missingGrowthData))
  #     # 
  #     # library(colf)
  #     # 
  #     # model1 = lm(data = trees[!trees$missingGrowthData,], hcb2013 ~ dbh2013) ; summary(model1)
  #     # model3 = colf_nlxb(data = trees[!trees$missingGrowthData,], formula = hcb2013 ~ dbh2013, lower = 0) ; summary(model3)
  #     # 
  #     # ggplot(trees, aes(x = dbh2013, y = hcb2013)) + geom_point(size= 2, aes(shape = missingGrowthData)) + coord_cartesian(xlim = c(0,NA), ylim = c(0,NA)) + 
  #     #   geom_abline(intercept = model3$coefficients[1], slope = model3$coefficients[2], alpha = 0.4)
  #     # 
  #     # # end test
  #     
  #     # linear model
  #     # "pour les arbres dont on connait le dbh (calculé à partir des accroissements et dbh 2013), pour un dbh2013 de X, on a un dbh de Y"
  #     trees = computeMissingVariableProportionnaly(trees, "dbh2013", "dbh", "sp", dbhIsMissingList, method = "proportionnal_mean")
  #     trees = computeMissingVariableProportionnaly(trees, "height2013", "height", "sp", dbhIsMissingList, method = "proportionnal_mean")
  #     trees = computeMissingVariableProportionnaly(trees, "hcb2013", "hcb", "sp", dbhIsMissingList, method = "proportionnal_mean")
  #     
  #     # # old way
  #     # # only valid dbh
  #     # validTrees = trees[!dbhIsMissingList, ]
  #     # 
  #     # meanTimeMachineDBHCorrection = mean(validTrees$dbh / validTrees$dbh2013)
  #     # meanTimeMachineBasalAreaCorrection = mean(validTrees$dbh**2 / validTrees$dbh2013**2)
  #     # 
  #     # # should be equivalent to meanTimeMachineDBHCorrection
  #     # meanTimeMachineHeightCorrection = mean(validTrees$height / validTrees$height2013)
  #     # meanTimeMachineHcbCorrection = mean(validTrees$hcb / validTrees$hcb2013)
  #     # 
  #     # trees$dbh[dbhIsNa] = round( trees$dbh2013[dbhIsNa] * meanTimeMachineDBHCorrection, 3 )
  #     # 
  #     # trees$height[dbhIsNa] = round( trees$height2013[dbhIsNa] * meanTimeMachineHeightCorrection, 3 )
  #     # trees$hcb[dbhIsNa] = round( trees$hcb2013[dbhIsNa] * meanTimeMachineHcbCorrection, 3 ) 
  #   }
  # }
  
  
  
  # Commented because these line were replicated below cr-2024.03.21
  # # 2-e CELLS TABLES
  # 
  # ###### Generating cells/ tirage des cellules
  # 
  # ###
  # IdCell=1:nbCell
  # # IdCell=seq(1:nbCell)
  # 
  # clign=rep(0,nbColumns)
  # for (i in 1:(nbLines-1)) {
  #   clign=c(clign,rep(i,nbColumns))
  # }
  # length(clign)
  # 
  # ccolbase = 0:(nbColumns-1)
  # # ccolbase= seq(0:(nbColumns-1))
  # # ccolbase=ccolbase-1
  # ccol=ccolbase
  # for (i in 1:(nbLines-1)) {
  #   ccol=c(ccol,ccolbase)
  # }
  # length(ccol)
  # 
  # soilHeightVec=rnorm(nbCell,mean=soilHeight, sd=VariabSoilHeith*soilHeight)
  # stoneContentVec = rep(stoneContent,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # wfcVec = rep(wfc,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # wiltVec = rep(wilt,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # propMacroVec = rep(propMacro,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # propMacroDeepVec = rep(propMacroDeep, nbCell)
  # bulkVec = rep(bulk,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # SOLCLAYtopVec = rep(SOLCLAYtop,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # SOLCLAYsolVec = rep(SOLCLAYsol,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # SOLFINtopVec = rep(SOLFINtop,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # SOLFINsolVec = rep(SOLFINsol,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # SOLSANDtopVec = rep(SOLSANDtop,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # SOLSANDsolVec = rep(SOLSANDsol,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # deepSoilDepthVec = rep(deepSoilDepth,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # stoneContentDeepVec = rep(stoneContentDeep,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # pracVec = rep(prac,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  # pracDeepVec = rep(pracDeep,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  
  
  
  
  
  
  
  # *****************************************
  # * STEP 3 = Computation / calculs ----
  # *****************************************
  
  library(MCMCpack) 
  
  #### COMPUTATION - SPECIES 1 
  
  # 3 Genetic sp1 ----
  
  ####Initialisation of allelic effects, FCRITBB / Initialisation des effets des alleles FCRITBB
  VAl_FCRITBB = (target_CV_FCRITBB*target_mean_FCRITBB)^2/nbLocus_FCRITBB;
  VAl_FCRITBB = VAl_FCRITBB*herit_FCRITBB;
  effect_FCRITBB = rep(VAl_FCRITBB, nbLocus_FCRITBB)
  
  allele1_FCRITBB = matrix(nrow=nbIndividuals[1], ncol = nbLocus_FCRITBB)
  allele2_FCRITBB = matrix(nrow=nbIndividuals[1], ncol = nbLocus_FCRITBB)
  FCRITBBval=matrix(nrow=nbIndividuals[1], ncol = 1 )
  
  if (herit_FCRITBB>0)
  {
    
    mean_FCRITBBval=0
    sd_FCRITBBval=0
    step=0
    while(( (abs(target_mean_FCRITBB-mean_FCRITBBval)/target_mean_FCRITBB)>0.01) ||  (abs(targetSd_FCRITBB-sd_FCRITBBval)/targetSd_FCRITBB>0.01) )
    {
      step=step+1
      print(step)
      
      freq_FCRITBB=rdirichlet(nbLocus_FCRITBB, c(1,1,1) )[,1]
      effect_FCRITBBdraw=sqrt(effect_FCRITBB/(freq_FCRITBB*(1-freq_FCRITBB)*2))
      
      for(l in 1:nbLocus_FCRITBB) {
        allele1_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_FCRITBB[l])
        allele2_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_FCRITBB[l])
      }
      rm(l)
      
      for (i in 1:nbIndividuals[1])
      {
        genetVal_FCRITBB=0
        for (j in 1:nbLocus_FCRITBB)
        {
          sumAll= allele1_FCRITBB[i,j] + allele2_FCRITBB[i,j]
          if (sumAll==2){
            genetVal_FCRITBB=genetVal_FCRITBB-effect_FCRITBBdraw[j]
          } else if (sumAll==4){
            genetVal_FCRITBB=genetVal_FCRITBB+effect_FCRITBBdraw[j]
          }
        }
        FCRITBBval[i]= genetVal_FCRITBB+target_mean_FCRITBB
      }
      rm(i,j,l)
      #FCRITBBvalOK=FCRITBBval[FCRITBBval>0]
      mean_FCRITBBval=mean(FCRITBBval)
      sd_FCRITBBval=sd(FCRITBBval)
      print(mean_FCRITBBval)
      print(sd_FCRITBBval)
    }
  } else {
    
    
    freq_FCRITBB=rdirichlet(nbLocus_FCRITBB, c(1,1,1) )[,1]
    for(l in 1:nbLocus_FCRITBB) {
      allele1_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_FCRITBB[l])
      allele2_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_FCRITBB[l])
    }
    FCRITBBval=rep(target_mean_FCRITBB, nbIndividuals[1])
    effect_FCRITBBdraw=freq_FCRITBB*0
    mean_FCRITBBval=target_mean_FCRITBB
    sd_FCRITBBval=0
    
  }
  
  # cr-17.01.2022 plot commented
  # par(mfrow=c(1,2))
  # ( FCRITBBval, main="", xlab="FCRITBB", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # abline(v=mean(FCRITBBval), lwd=2, col="red")
  # mean_FCRITBBval=mean(FCRITBBval)
  # mtext(text=paste("mean_FCRITBB=", round(mean_FCRITBBval,2), sep=" "), side= 3, col="red", cex=1.2)
  # 
  # hist( effect_FCRITBBdraw, main="", xlab="Allelic effects", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # 
  # effect_FCRITBB=round(alleleEffectMultiplCoeffFCRITBB*effect_FCRITBBdraw)
  
  ####Initialisation of allelic effects, g1max / Initialisation des effets des alleles g1max
  VAl_g1max = (target_CV_g1max*target_mean_g1max)^2/nbLocus_g1max;
  VAl_g1max=VAl_g1max*herit_g1max;
  effect_g1max=rep(VAl_g1max,nbLocus_FCRITBB)
  
  allele1_g1max=matrix(nrow=nbIndividuals[1], ncol = (nbLocus_g1max))
  allele2_g1max=matrix(nrow=nbIndividuals[1], ncol = (nbLocus_g1max))
  g1maxval=matrix(nrow=nbIndividuals[1], ncol = 1 )
  
  if (herit_g1max>0)
  {
    mean_g1maxVal=0
    sd_g1maxVal=0
    while( (abs(target_mean_g1max-mean_g1maxVal)/target_mean_g1max>0.01) || (abs(targetSd_g1max-sd_g1maxVal)/targetSd_g1max>0.01))
    {
      step=step+1
      print(step)
      
      freq_g1max=rdirichlet(nbLocus_g1max, c(1,1,1) )[,1]
      effect_g1maxdraw=sqrt(effect_g1max/(freq_g1max*(1-freq_g1max)*2))
      
      for(l in 1:nbLocus_g1max) {
        allele1_g1max[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_g1max[l])
        allele2_g1max[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_g1max[l])
      }
      
      for (i in 1:nbIndividuals[1])
      {
        genetVal_g1max=0
        for (j in 1:nbLocus_g1max)
        {
          sumAll= allele1_g1max[i,j] + allele2_g1max[i,j]
          if (sumAll==2){
            genetVal_g1max=genetVal_g1max-effect_g1maxdraw[j]
          } else if (sumAll==4){
            genetVal_g1max=genetVal_g1max+effect_g1maxdraw[j]
          }
        }
        g1maxval[i]= genetVal_g1max+target_mean_g1max
      }
      rm(i,j,l)
      mean_g1maxVal=mean(g1maxval)
      print(mean_g1maxVal)
      sd_g1maxVal=sd(g1maxval)
      print(sd_g1maxVal)
      
    }
  } else {
    freq_g1max=rdirichlet(nbLocus_g1max, c(1,1,1) )[,1]
    for(l in 1:nbLocus_g1max) {
      allele1_g1max[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_g1max[l])
      allele2_g1max[,l]<- 1+ rbinom(nbIndividuals[1], size=1, prob=freq_g1max[l])
    }
    g1maxval=rep(target_mean_g1max, nbIndividuals[1])
    # effect_g1maxdraw=freq_g1maxNbOfGrowingSeedlingsmax*0 # freq_g1maxNbOfGrowingSeedlingsmax introuvable
    effect_g1maxdraw= 0 # freq_g1maxNbOfGrowingSeedlingsmax introuvable
    mean_g1maxVal=target_mean_g1max
    sd_g1maxVal=0
    
  }
  
  # cr-17.01.2022 plot commented
  # par(mfrow=c(1,2))
  # hist( g1maxval, main="", xlab="g1max", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # abline(v=mean(g1maxval), lwd=2, col="red")
  # mean_g1max=mean(g1maxval)
  # sd_g1max=sd(g1maxval)
  # mtext(text=paste("mean_g1max=", round(mean_g1max,2), sep=" "), side= 3, col="red", cex=1.2)
  # hist( effect_g1maxdraw, main="", xlab="Allelic effects", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # 
  # effect_g1max=round(alleleEffectMultiplCoeffg1max*effect_g1maxdraw)
  
  genetVal=as.data.frame(cbind(FCRITBBval,g1maxval))
  names(genetVal)=c("FCRITBB","g1max")
  
  # cr-01.06.2022 write genet value commented
  # write.table(x = genetVal, file="genetValueG0.txt")
  
  # 3 TREES sp1----
  
  ##########Generating individuals/ tirage des individus - species 1
  
  #genotypes
  allele1Neutral=matrix(nrow=nbIndividuals[1], ncol = (nbMsat+nbSNP) )
  allele2Neutral=matrix(nrow=nbIndividuals[1], ncol = (nbMsat+nbSNP) )
  
  if(nbMsat>0){
    for(i in 1:nbMsat) {
      allele1Neutral[,i]<- 11+round(runif(nbIndividuals[1], min=0, max=9))
      allele2Neutral[,i]<- 11+round(runif(nbIndividuals[1], min=0, max=9))
    }
    rm(i)
    #allele1Neutral[,nbLocus_FCRITBB+nbLocusg1max+1]
  }
  
  freqSNP=runif(nbSNP, min=0, max=1)
  for(i in 1:nbSNP) {
    allele1Neutral[,nbMsat+i]<- 1+rbinom(nbIndividuals[1], size=1, prob=freqSNP[i])
    allele2Neutral[,nbMsat+i]<- 1+rbinom(nbIndividuals[1], size=1, prob=freqSNP[i])
  }
  rm(i)
  #allele1Neutral[,nbLocus_FCRITBB+ nbLocusg1max+nbMsat+1]
  
  allele1_1=cbind(allele1_FCRITBB, allele1_g1max,allele1Neutral)
  allele2_1=cbind(allele2_FCRITBB, allele2_g1max,allele2Neutral)
  
  
  
  #### COMPUTATION - SPECIES 2 
  
  # 3 Genetic sp2 ----
  
  ####Initialisation of allelic effects, FCRITBB / Initialisation des effets des alleles FCRITBB
  VAl_FCRITBB = (target_CV_FCRITBB*target_mean_FCRITBB)^2/nbLocus_FCRITBB;
  VAl_FCRITBB=VAl_FCRITBB*herit_FCRITBB;
  effect_FCRITBB=rep(VAl_FCRITBB,nbLocus_FCRITBB)
  
  allele1_FCRITBB=matrix(nrow=nbIndividuals[2], ncol = nbLocus_FCRITBB)
  allele2_FCRITBB=matrix(nrow=nbIndividuals[2], ncol = nbLocus_FCRITBB)
  FCRITBBval=matrix(nrow=nbIndividuals[2], ncol = 1 )
  
  if (herit_FCRITBB>0){
    
    mean_FCRITBBval=0
    sd_FCRITBBval=0
    step=0
    while(( (abs(target_mean_FCRITBB-mean_FCRITBBval)/target_mean_FCRITBB)>0.01) ||  (abs(targetSd_FCRITBB-sd_FCRITBBval)/targetSd_FCRITBB>0.01) )
    {
      step=step+1
      print(step)
      
      freq_FCRITBB=rdirichlet(nbLocus_FCRITBB, c(1,1,1) )[,1]
      effect_FCRITBBdraw=sqrt(effect_FCRITBB/(freq_FCRITBB*(1-freq_FCRITBB)*2))
      
      for(l in 1:nbLocus_FCRITBB) {
        allele1_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_FCRITBB[l])
        allele2_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_FCRITBB[l])
      }
      
      for (i in 1:nbIndividuals[2])
      {
        genetVal_FCRITBB=0
        for (j in 1:nbLocus_FCRITBB)
        {
          sumAll= allele1_FCRITBB[i,j] + allele2_FCRITBB[i,j]
          if (sumAll==2){
            genetVal_FCRITBB=genetVal_FCRITBB-effect_FCRITBBdraw[j]
          } else if (sumAll==4){
            genetVal_FCRITBB=genetVal_FCRITBB+effect_FCRITBBdraw[j]
          }
        }
        FCRITBBval[i]= genetVal_FCRITBB+target_mean_FCRITBB
      }
      #FCRITBBvalOK=FCRITBBval[FCRITBBval>0]
      mean_FCRITBBval=mean(FCRITBBval)
      sd_FCRITBBval=sd(FCRITBBval)
      print(mean_FCRITBBval)
      print(sd_FCRITBBval)
    }
  } else {
    
    
    freq_FCRITBB=rdirichlet(nbLocus_FCRITBB, c(1,1,1) )[,1]
    for(l in 1:nbLocus_FCRITBB) {
      allele1_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_FCRITBB[l])
      allele2_FCRITBB[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_FCRITBB[l])
    }
    FCRITBBval=rep(target_mean_FCRITBB, nbIndividuals[2])
    effect_FCRITBBdraw=freq_FCRITBB*0
    mean_FCRITBBval=target_mean_FCRITBB
    sd_FCRITBBval=0
    
  }
  
  # cr-17.01.2022 plot commented
  # par(mfrow=c(1,2))
  # hist( FCRITBBval, main="", xlab="FCRITBB", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # abline(v=mean(FCRITBBval), lwd=2, col="red")
  # mean_FCRITBBval=mean(FCRITBBval)
  # mtext(text=paste("mean_FCRITBB=", round(mean_FCRITBBval,2), sep=" "), side= 3, col="red", cex=1.2)
  # 
  # hist( effect_FCRITBBdraw, main="", xlab="Allelic effects", cex.lab=1.5, cex.axis=1.2, breaks=20)
  
  effect_FCRITBB=round(alleleEffectMultiplCoeffFCRITBB*effect_FCRITBBdraw)
  
  ####Initialisation of allelic effects, g1max / Initialisation des effets des alleles g1max
  VAl_g1max = (target_CV_g1max*target_mean_g1max)^2/nbLocus_g1max;
  VAl_g1max=VAl_g1max*herit_g1max;
  effect_g1max=rep(VAl_g1max,nbLocus_FCRITBB)
  
  allele1_g1max=matrix(nrow=nbIndividuals[2], ncol = (nbLocus_g1max))
  allele2_g1max=matrix(nrow=nbIndividuals[2], ncol = (nbLocus_g1max))
  g1maxval=matrix(nrow=nbIndividuals[2], ncol = 1 )
  
  if (herit_g1max>0)
  {
    mean_g1maxVal=0
    sd_g1maxVal=0
    while( (abs(target_mean_g1max-mean_g1maxVal)/target_mean_g1max>0.01) || (abs(targetSd_g1max-sd_g1maxVal)/targetSd_g1max>0.01))
    {
      step=step+1
      print(step)
      
      freq_g1max=rdirichlet(nbLocus_g1max, c(1,1,1) )[,1]
      effect_g1maxdraw=sqrt(effect_g1max/(freq_g1max*(1-freq_g1max)*2))
      
      for(l in 1:nbLocus_g1max) {
        allele1_g1max[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_g1max[l])
        allele2_g1max[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_g1max[l])
      }
      
      for (i in 1:nbIndividuals[2])
      {
        genetVal_g1max=0
        for (j in 1:nbLocus_g1max)
        {
          sumAll= allele1_g1max[i,j] + allele2_g1max[i,j]
          if (sumAll==2){
            genetVal_g1max=genetVal_g1max-effect_g1maxdraw[j]
          } else if (sumAll==4){
            genetVal_g1max=genetVal_g1max+effect_g1maxdraw[j]
          }
        }
        g1maxval[i]= genetVal_g1max+target_mean_g1max
      }
      rm(l,i,j)
      mean_g1maxVal=mean(g1maxval)
      print(mean_g1maxVal)
      sd_g1maxVal=sd(g1maxval)
      print(sd_g1maxVal)
      
    }
  } else {
    freq_g1max=rdirichlet(nbLocus_g1max, c(1,1,1) )[,1]
    for(l in 1:nbLocus_g1max) {
      allele1_g1max[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_g1max[l])
      allele2_g1max[,l]<- 1+ rbinom(nbIndividuals[2], size=1, prob=freq_g1max[l])
    }
    g1maxval=rep(target_mean_g1max, nbIndividuals[2])
    # effect_g1maxdraw=freq_g1maxNbOfGrowingSeedlingsmax*0 # freq_g1maxNbOfGrowingSeedlingsmax introuvable
    effect_g1maxdraw=0 # freq_g1maxNbOfGrowingSeedlingsmax introuvable
    mean_g1maxVal=target_mean_g1max
    sd_g1maxVal=0
    
  }
  
  # cr-17.01.2022 plot commented
  # par(mfrow=c(1,2))
  # hist( g1maxval, main="", xlab="g1max", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # abline(v=mean(g1maxval), lwd=2, col="red")
  # mean_g1max=mean(g1maxval)
  # sd_g1max=sd(g1maxval)
  # mtext(text=paste("mean_g1max=", round(mean_g1max,2), sep=" "), side= 3, col="red", cex=1.2)
  # hist( effect_g1maxdraw, main="", xlab="Allelic effects", cex.lab=1.5, cex.axis=1.2, breaks=20)
  # 
  # effect_g1max=round(alleleEffectMultiplCoeffg1max*effect_g1maxdraw)
  
  genetVal=as.data.frame(cbind(FCRITBBval,g1maxval))
  names(genetVal)=c("FCRITBB","g1max")
  # cr-01.06.2022 write genet value commented
  # write.table(x = genetVal, file="genetValueG0.txt")
  
  # 3 TREES sp2 ----
  
  ##########Generating individuals/ tirage des individus - species 2
  
  #genotypes
  allele1Neutral=matrix(nrow=nbIndividuals[2], ncol = (nbMsat+nbSNP) )
  allele2Neutral=matrix(nrow=nbIndividuals[2], ncol = (nbMsat+nbSNP) )
  
  if(nbMsat>0){
    for(i in 1:nbMsat) {
      allele1Neutral[,i]<- 11+round(runif(nbIndividuals[2], min=0, max=9))
      allele2Neutral[,i]<- 11+round(runif(nbIndividuals[2], min=0, max=9))
    }
    rm(i)
    #allele1Neutral[,nbLocus_FCRITBB+nbLocusg1max+1]
  }
  
  freqSNP=runif(nbSNP, min=0, max=1)
  for(i in 1:nbSNP) {
    allele1Neutral[,nbMsat+i]<- 1+rbinom(nbIndividuals[2], size=1, prob=freqSNP[i])
    allele2Neutral[,nbMsat+i]<- 1+rbinom(nbIndividuals[2], size=1, prob=freqSNP[i])
  }
  rm(i)
  #allele1Neutral[,nbLocus_FCRITBB+ nbLocusg1max+nbMsat+1]
  
  allele1_2=cbind(allele1_FCRITBB, allele1_g1max, allele1Neutral)
  allele2_2=cbind(allele2_FCRITBB, allele2_g1max, allele2Neutral)
  
  # end generating individuals 
  
  
  # Genetics container (to improve for multispecies)
  pdgInventoryOptions$genetics_variables = tibble(mean_FCRITBBval = mean_FCRITBBval, sd_FCRITBBval = sd_FCRITBBval,
                                                  mean_g1maxVal = mean_g1maxVal, sd_g1maxVal = sd_g1maxVal)
  
  pdgInventoryOptions$effect_FCRITBBdraw = effect_FCRITBBdraw
  pdgInventoryOptions$effect_g1maxdraw = effect_g1maxdraw
  
  # 3 CELLS ----
  
  ###### Generating cells/ tirage des cellules
  
  ###
  IdCell=1:nbCell
  # IdCell=seq(1:nbCell)
  
  clign=rep(0,nbColumns)
  for (i in 1:(nbLines-1)) {
    clign=c(clign,rep(i,nbColumns))
  }
  rm(i)
  length(clign)
  
  ccolbase = 0:(nbColumns-1)
  # ccolbase= seq(0:(nbColumns-1))
  # ccolbase=ccolbase-1
  ccol=ccolbase
  for (i in 1:(nbLines-1)) {
    ccol=c(ccol,ccolbase)
  }
  rm(i)
  length(ccol)
  
  soilHeightVec=rnorm(nbCell,mean=soilHeight, sd=VariabSoilHeith*soilHeight)
  stoneContentVec = rep(stoneContent,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  wfcVec = rep(wfc,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  wiltVec = rep(wilt,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  propMacroVec = rep(propMacro,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  propMacroDeepVec = rep(propMacroDeep, nbCell)
  bulkVec = rep(bulk,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  SOLCLAYtopVec = rep(SOLCLAYtop,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  SOLCLAYsolVec = rep(SOLCLAYsol,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  SOLFINtopVec = rep(SOLFINtop,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  SOLFINsolVec = rep(SOLFINsol,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  SOLSANDtopVec = rep(SOLSANDtop,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  SOLSANDsolVec = rep(SOLSANDsol,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  deepSoilDepthVec = rep(deepSoilDepth,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  stoneContentDeepVec = rep(stoneContentDeep,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  pracVec = rep(prac,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  pracDeepVec = rep(pracDeep,nbCell)     #percentage of stone in the soil / % de pierre dans le sol
  
  # cr-17.01.2022 plot commented
  # #plotting the distribution of water holding capacity
  # hist(soilHeightVec,main="", xlab="Soil height (mm)", cex.lab=1.5, cex.axis=1.2)
  # WHC=(1-stoneContent)*soilHeightVec*(wfc-wilt)*bulk
  # hist(WHC, main="", xlab="Water holding capacity (mm)", cex.lab=1.5, cex.axis=1.2)
  
  
  
  cells = tibble(IdCell, clign, ccol, soilHeight = soilHeightVec, stoneContent = stoneContentVec, wfc = wfcVec, wilt = wiltVec,
                 propMacro = propMacroVec, propMacroDeep = propMacroDeepVec, bulk = bulkVec, 
                 SOLCLAYtop = SOLCLAYtopVec, SOLCLAYsol = SOLCLAYsolVec, SOLFINtop = SOLFINtopVec, SOLFINsol = SOLFINsolVec,
                 SOLSANDtop = SOLSANDtopVec, SOLSANDsol = SOLSANDsolVec, 
                 deepSoilDepth = deepSoilDepthVec, stoneContentDeep = stoneContentDeepVec, prac = pracVec, pracDeep = pracDeepVec)
  
  
  
  
  # *****************************************
  # * STEP 4 = Print the inventory file / imprimer dans le fichier d'inventaire ----
  # * Creation of regularized plot if asked
  # *****************************************
  
  # Print E2 ----
  
  # write demography-fine plurispecific inventory is all cases
  inventoryName = paste0(code_site,"_PDG_irregdemo_plurisp.inv")

  pdgInventoryDescription = paste0("## Capsis 4.1 - module PDG", 
                                   "\n# PDG Inventory based on GMAP data, code_site: ", code_site, 
         "\n# Date and time: ", paste(Sys.time()), 
         "\n\n# Ignored dead trees: " , dim(subset(dGMAP_hp_site_save_mort, mort == 1))[1], " over ", dim(dGMAP_hp_site_save_mort)[1], " trees",
         "\n# Ignored living trees of other species: " , dim(subset(dGMAP_hp_site_save_vivant_allspecies, ! essence %in% c("hetre", "sapin")))[1],
      " over ",  dim(dGMAP_hp_site_save_vivant_allspecies)[1] , " living trees", 
      "\n# ", dim(dGMAP_hp_site)[1], " selected living trees of species : ", paste(demographic_parameters$speciesFrenchNames, collapse = " "))
  
  if(useTimeMachine){
    pdgInventoryDescription = paste0(pdgInventoryDescription, "\n## Time machine was used")
  }
  
  
  pdgInventoryOptions$pdgPlotParameters = tibble(nlin = nbLines, ncol = nbColumns, cellWidth = cellWidth,
                                                 xorigin = xorigin, yorigin = yorigin,
                                                 LAI = LAI_LIDAR)
  
  
  
  
  writeInventoryPDG(outputInventoryFolder, inventoryName, pdgInventoryDescription,
                    cells, trees, pdgInventoryOptions, altitudeClasses, 
                    allele1_1, allele2_1, allele1_2, allele2_2, thisScriptOptions)
  
  
  
  # fill inventoryTable
  inventoryTableIndex = which(inventoryTable$code_site == code_site & inventoryTable$inventoryType == "irregdemo_plurisp")
  inventoryTable[inventoryTableIndex, ]$created = TRUE
  inventoryTable[inventoryTableIndex, ]$inventoryName = inventoryName
  trees_beech = subset(trees, sp == 1)
  trees_fir = subset(trees, sp == 2)
  inventoryTable[inventoryTableIndex, ]$bee_nprop = dim(trees_beech)[1] / dim(trees)[1]
  inventoryTable[inventoryTableIndex, ]$fir_nprop =  dim(trees_fir)[1] / dim(trees)[1]
  inventoryTable[inventoryTableIndex, ]$bee_BAprop = sum(trees_beech$dbh ** 2) / sum(trees$dbh **2)
  inventoryTable[inventoryTableIndex, ]$fir_BAprop = sum(trees_fir$dbh ** 2) / sum(trees$dbh **2)
  
  standArea = cellWidth**2 * nbLines * nbColumns
  inventoryTable[inventoryTableIndex, ]$GHA = sum(pi*(trees$dbh/2) ** 2) / standArea
  inventoryTable[inventoryTableIndex, ]$bee_dbh = sqrt(mean(trees_beech$dbh**2))
  inventoryTable[inventoryTableIndex, ]$fir_dbh = sqrt(mean(trees_fir$dbh**2))
  
  # fill values_init for all code site
  inventoryTableIndex = which(inventoryTable$code_site == code_site)
  inventoryTable[inventoryTableIndex, ]$bee_BAprop_init = sum(trees_beech$dbh ** 2) / sum(trees$dbh **2)
  inventoryTable[inventoryTableIndex, ]$fir_BAprop_init = sum(trees_fir$dbh ** 2) / sum(trees$dbh **2)
  inventoryTable[inventoryTableIndex, ]$GHA_init = sum(pi*(trees$dbh/2) ** 2) / standArea
  inventoryTable[inventoryTableIndex, ]$bee_dbh_init = sqrt(mean(trees_beech$dbh**2))
  inventoryTable[inventoryTableIndex, ]$fir_dbh_init = sqrt(mean(trees_fir$dbh**2))
  
  
  # Print E0 ----
  if(inventoryMode == "regdemo_monosp" | inventoryMode == "all" | inventoryMode == "allbutE1B"){
    
      
    for(idSpecies in unique(trees$sp)){
      
      # # if species is minority
      # spName = speciesFrenchNames[which(idSp == idSpecies)]
      # if(grepl(code_site, pattern = "_sp_") | grepl(code_site, pattern = "_ph_")){
      #   next
      # }
      
      # get regularized inventories
      res = computeRegDemoMonoSpInventory(treesTable = trees, cellsTable = cells, idSpecies = idSpecies, 
                                          allele1_1, allele2_1, allele1_2, allele2_2, 
                                          pdgInventoryOptions, CASTANEASpeciesFile, thisScriptOptions, LAI_LIDAR)
      
      # is there is only one species
      if (length(unique(trees$sp)) <= 1){
        # write nothing, regdemo_plurisp will be equivalent 
      }else{
        # write regularized PDG inventory
        inventoryName = paste0(code_site, "_", speciesFrenchNames[which(idSp == idSpecies)], "_PDG_regdemo_monosp", ".inv")
        pdgInventoryDescription = paste0("## Capsis 4.1 - module PDG", 
                                         "\n# PDG Inventory based on GMAP data, code_site: ", code_site, 
                                         "\n# Date and time: ", paste(Sys.time()), 
                                         "\n\n# Regularized monospecific inventory with ", demographic_parameters$speciesFullNames[ match(idSpecies, demographic_parameters$idSp)  ])
        if(useTimeMachine){
          pdgInventoryDescription = paste0(pdgInventoryDescription, "\n## Time machine was used")
        }
        
        pdgInventoryOptions$pdgPlotParameters = tibble(nlin = res$nbLines_sp, ncol = res$nbColumns_sp, cellWidth = res$cellWidth_sp,
                                                       xorigin = xorigin, yorigin = yorigin,
                                                       LAI = res$LAI_reg)
        
        writeInventoryPDG(outputInventoryFolder, inventoryName,  pdgInventoryDescription, res$cells_sp, res$trees_sp, pdgInventoryOptions,
                          res$altitudeClasses_sp, 
                          res$allele1_1, res$allele2_1, res$allele1_2, res$allele2_2, thisScriptOptions)
        
        # fill inventoryTable
        inventoryType = paste0("regdemo_monosp_", speciesFrenchNames[which(idSp == idSpecies)] )
        inventoryTableIndex = which(inventoryTable$code_site == code_site & inventoryTable$inventoryType == inventoryType)
        inventoryTable[inventoryTableIndex, ]$created = TRUE
        inventoryTable[inventoryTableIndex, ]$inventoryName = inventoryName
        trees_beech = subset(res$trees_sp, sp == 1)
        trees_fir = subset(res$trees_sp, sp == 2)
        inventoryTable[inventoryTableIndex, ]$bee_nprop = dim(trees_beech)[1] / dim(res$trees_sp)[1]
        inventoryTable[inventoryTableIndex, ]$fir_nprop =  dim(trees_fir)[1] / dim(res$trees_sp)[1]
        inventoryTable[inventoryTableIndex, ]$bee_BAprop = sum(trees_beech$dbh ** 2) / sum(res$trees_sp$dbh **2)
        inventoryTable[inventoryTableIndex, ]$fir_BAprop = sum(trees_fir$dbh ** 2) / sum(res$trees_sp$dbh **2)
        
        standArea = res$cellWidth_sp**2 * res$nbLines_sp * res$nbColumns_sp
        inventoryTable[inventoryTableIndex, ]$GHA = sum(pi*(res$trees_sp$dbh/2) ** 2) / standArea
        
        inventoryTable[inventoryTableIndex, ]$bee_dbh = sqrt(mean(trees_beech$dbh**2))
        inventoryTable[inventoryTableIndex, ]$fir_dbh = sqrt(mean(trees_fir$dbh**2))
      }
      
      # Print CASTANEA ----
      
      inventoryName = paste0(code_site, "_", speciesFrenchNames[which(idSp == idSpecies)], "_CASTANEA", ".inv")
      
      description = paste0("## Capsis 4.1 - module CASTANEAonly",
                           "\n# CASTANEA inventory based on GMAP data, code_site: ", code_site,
                           "\n# Date and time: ", paste(Sys.time()))
      
      if(useTimeMachine){
        description = paste0(description, "\n## Time machine was used")
      }
      
      pdgInventoryOptionsCAST = pdgInventoryOptions
      pdgInventoryOptionsCAST$plotLocationParameters$slope_deg = 0
      pdgInventoryOptionsCAST$plotLocationParameters$aspect_deg = NA
      pdgInventoryOptionsCAST$pdgPlotParameters = NULL
      
      standAreaCAST = dim(res$cells_sp)[1] * res$cellWidth_sp **2
      
      writeInventoryCASTANEAFromRegdemoMonosp(inventoryFolder = outputInventoryFolder, inventoryName = inventoryName, 
                                              description = description,
                                              castaneaSpeciesCode = demographic_parameters$fmIdSp[which(demographic_parameters$idSp == idSpecies)],
                                              pdgInventoryOptions = pdgInventoryOptionsCAST,
                                              cellsTable = res$cells_sp, treesTable = res$trees_sp, 
                                              standArea = standAreaCAST,
                                              initYear = initYear,
                                              outputInConsole = outputInConsole,
                                              LAI = res$LAI_reg)
      
      # fill inventoryTable
      inventoryType = paste0("CASTANEA_", speciesFrenchNames[which(idSp == idSpecies)] )
      inventoryTableIndex = which(inventoryTable$code_site == code_site & inventoryTable$inventoryType == inventoryType)
      inventoryTable[inventoryTableIndex, ]$created = TRUE
      inventoryTable[inventoryTableIndex, ]$inventoryName = inventoryName
      if(idSpecies == idSp[which(speciesFrenchNames == "hetre")]){
        inventoryTable[inventoryTableIndex, ]$bee_nprop = 1
        inventoryTable[inventoryTableIndex, ]$fir_nprop =  0
        inventoryTable[inventoryTableIndex, ]$bee_BAprop = 1
        inventoryTable[inventoryTableIndex, ]$fir_BAprop = 0
        inventoryTable[inventoryTableIndex, ]$bee_dbh = sqrt(mean(res$trees_sp$dbh**2))
        inventoryTable[inventoryTableIndex, ]$fir_dbh = NA
        
      }else if(idSpecies == idSp[which(speciesFrenchNames == "sapin")]){
        inventoryTable[inventoryTableIndex, ]$bee_nprop = 0
        inventoryTable[inventoryTableIndex, ]$fir_nprop =  1
        inventoryTable[inventoryTableIndex, ]$bee_BAprop = 0
        inventoryTable[inventoryTableIndex, ]$fir_BAprop = 1
        inventoryTable[inventoryTableIndex, ]$bee_dbh = NA
        inventoryTable[inventoryTableIndex, ]$fir_dbh = sqrt(mean(res$trees_sp$dbh**2))
      }
      
      inventoryTable[inventoryTableIndex, ]$GHA = sum( pi * (res$trees_sp$dbh/2)**2 ) / standAreaCAST
      
    }
  
  }
  
  # Print E1A ----
  if (inventoryMode == "regdemo_plurisp" | inventoryMode == "all" | inventoryMode == "allbutE1B"){
    
    
    # cr-27.09.2023 : dans tous les cas, je fais regdemo plurisp
    # # check the quantity of trees per species
    # basalArea_sp = rep(x = 0, times = length(unique(trees$sp)))
    # for(itree in 1:dim(trees)[1]){
    #   atree = trees[itree, ]
    #   index_sp = which(unique(trees$sp) == atree$sp)
    #   basalAreaTree = atree$dbh**2  # approx
    #   basalArea_sp[index_sp] = basalArea_sp[index_sp] + basalAreaTree
    # }
    # proportion_sp = basalArea_sp / sum(basalArea_sp)
    # nTreeRound_sp = round(proportion_sp * nTree_regularization)
    # 
    # thereIsSpeciesWithZeroTrees = 0 %in% nTreeRound_sp
    
    res = computeRegDemoPluriSpInventory(trees, cells, allele1_1, allele2_1, allele1_2, allele2_2, 
                                         pdgInventoryOptions, CASTANEASpeciesFile, thisScriptOptions, LAI_LIDAR)
    
    # write regularized PDG inventory
    inventoryName = paste0(code_site, "_PDG_regdemo_plurisp", ".inv")
    pdgInventoryDescription = paste0("## Capsis 4.1 - module PDG", 
                                     "\n# PDG Inventory based on GMAP data, code_site: ", code_site, 
                                     "\n# Date and time: ", paste(Sys.time()), 
                                     "\n\n# Regularized plurispecific inventory")
    
    if(useTimeMachine){
      pdgInventoryDescription = paste0(pdgInventoryDescription, "\n## Time machine was used")
    }
    
    pdgInventoryOptions$pdgPlotParameters = tibble(nlin = res$nbLines_new, ncol = res$nbColumns_new, cellWidth = res$cellWidth_new,
                                                   xorigin = xorigin, yorigin = yorigin,
                                                   LAI = LAI_LIDAR)
    
    writeInventoryPDG(outputInventoryFolder, inventoryName,  pdgInventoryDescription, res$cells_new, res$trees_new, pdgInventoryOptions,
                      res$altitudeClasses_new, 
                      res$allele1_1, res$allele2_1, res$allele1_2, res$allele2_2, thisScriptOptions)
    
    
    # fill inventoryTable
    inventoryTableIndex = which(inventoryTable$code_site == code_site & inventoryTable$inventoryType == "regdemo_plurisp")
    inventoryTable[inventoryTableIndex, ]$created = TRUE
    inventoryTable[inventoryTableIndex, ]$inventoryName = inventoryName
    trees_beech = subset(res$trees_new, sp == 1)
    trees_fir = subset(res$trees_new, sp == 2)
    inventoryTable[inventoryTableIndex, ]$bee_nprop = dim(trees_beech)[1] / dim(res$trees_new)[1]
    inventoryTable[inventoryTableIndex, ]$fir_nprop =  dim(trees_fir)[1] / dim(res$trees_new)[1]
    inventoryTable[inventoryTableIndex, ]$bee_BAprop = sum(trees_beech$dbh ** 2) / sum(res$trees_new$dbh **2)
    inventoryTable[inventoryTableIndex, ]$fir_BAprop = sum(trees_fir$dbh ** 2) / sum(res$trees_new$dbh **2)
    
    standArea = res$cellWidth_new**2 * res$nbLines_new * res$nbColumns_new
    inventoryTable[inventoryTableIndex, ]$GHA = sum(pi*(res$trees_new$dbh/2) ** 2) / standArea
    
    inventoryTable[inventoryTableIndex, ]$bee_dbh = sqrt(mean(trees_beech$dbh**2))
    inventoryTable[inventoryTableIndex, ]$fir_dbh = sqrt(mean(trees_fir$dbh**2))

  }
  
  # Print E1B ----
  
  if (inventoryMode == "irregdemo_monosp" | inventoryMode == "all" ){
    if (length(unique(trees$sp)) <= 1){
      # do nothing, plurisp will be equivalent 
    }else{
      for(idSpecies in unique(trees$sp)){
        
        # if species is minority
        # spName = speciesFrenchNames[which(idSp == idSpecies)]
        # if(grepl(code_site, pattern = "_sp_") | grepl(code_site, pattern = "_ph_")){
        #   next
        # }
        
        
        
        
        # get an irregular monospecific-rendered inventory
        res = computeIrregDemoMonoSpInventory(treesTable = trees, cellsTable = cells, idSpecies = idSpecies, 
                                              allele1_1, allele2_1, allele1_2, allele2_2, pdgInventoryOptions, CASTANEASpeciesFile, LAI_LIDAR)

        pdgInventoryOptions$pdgPlotParameters = tibble(nlin = nbLines, ncol = nbColumns, cellWidth = cellWidth,
                                                       xorigin = xorigin, yorigin = yorigin,
                                                       LAI = res$LAI_sp)
        
        # write irregular monospecific-rendered PDG inventory
        inventoryName = paste0(code_site, "_", speciesFrenchNames[which(idSp == idSpecies)], "_PDG_irregdemo_monosp", ".inv")
        pdgInventoryDescription = paste0("## Capsis 4.1 - module PDG", 
                                         "\n# PDG Inventory based on GMAP data, code_site: ", code_site, 
                                         "\n# Date and time: ", paste(Sys.time()), 
                                         "\n\n# Irregular monospecific inventory with ", demographic_parameters$speciesFullNames[ match(idSpecies, demographic_parameters$idSp)  ])
        if(useTimeMachine){
          pdgInventoryDescription = paste0(pdgInventoryDescription, "\n## Time machine was used")
        }
        
        
        writeInventoryPDG(outputInventoryFolder, inventoryName, pdgInventoryDescription, cells, res$trees_sp, pdgInventoryOptions,
                          altitudeClasses, 
                          res$allele1_1, res$allele2_1, res$allele1_2, res$allele2_2, thisScriptOptions)
        
        
        # fill inventoryTable
        inventoryType = paste0("irregdemo_monosp_", speciesFrenchNames[which(idSp == idSpecies)] )
        inventoryTableIndex = which(inventoryTable$code_site == code_site & inventoryTable$inventoryType == inventoryType)
        inventoryTable[inventoryTableIndex, ]$created = TRUE
        inventoryTable[inventoryTableIndex, ]$inventoryName = inventoryName
        trees_beech = subset(res$trees_sp, sp == 1)
        trees_fir = subset(res$trees_sp, sp == 2)
        inventoryTable[inventoryTableIndex, ]$bee_nprop = dim(trees_beech)[1] / dim(res$trees_sp)[1]
        inventoryTable[inventoryTableIndex, ]$fir_nprop =  dim(trees_fir)[1] / dim(res$trees_sp)[1]
        inventoryTable[inventoryTableIndex, ]$bee_BAprop = sum(trees_beech$dbh ** 2) / sum(res$trees_sp$dbh **2)
        inventoryTable[inventoryTableIndex, ]$fir_BAprop = sum(trees_fir$dbh ** 2) / sum(res$trees_sp$dbh **2)
        
        standArea = cellWidth**2 * nbLines * nbColumns
        inventoryTable[inventoryTableIndex, ]$GHA = sum(pi*(res$trees_sp$dbh/2) ** 2) / standArea
        
        
        trees_true_beech = subset(res$trees_sp, sp == 1 & idTree < 100000)
        trees_true_fir = subset(res$trees_sp, sp == 2 & idTree < 10000)
        inventoryTable[inventoryTableIndex, ]$bee_dbh = sqrt(mean(trees_true_beech$dbh**2))
        inventoryTable[inventoryTableIndex, ]$fir_dbh = sqrt(mean(trees_true_fir$dbh**2))
        
      }
    }
  } # end cases of inventoryMode
  
  
} # end loop on sites
cat("\nscript() finished\n\n")
print(warnings())


# WRITE INVENTORY TABLE ----

inventoryTable[is.nan(inventoryTable$bee_dbh), ]$bee_dbh = NA
inventoryTable[is.nan(inventoryTable$fir_dbh), ]$fir_dbh = NA
inventoryTable[is.nan(inventoryTable$bee_dbh_init), ]$bee_dbh_init = NA
inventoryTable[is.nan(inventoryTable$fir_dbh_init), ]$fir_dbh_init = NA

inventoryTable$comp = ifelse(grepl(inventoryTable$code_site, pattern = "_m"), "m", 
                                    ifelse(grepl(inventoryTable$code_site, pattern = "_ph"), "ph", 
                                           ifelse(grepl(inventoryTable$code_site, pattern = "_sp"), "sp", "NA" 
                                           )))

inventoryTable$diffGHA = abs(inventoryTable$GHA - inventoryTable$GHA_init)
inventoryTable$diffBeeBAprop = abs(inventoryTable$bee_BAprop - inventoryTable$bee_BAprop_init)
inventoryTable$diffFirBAprop = abs(inventoryTable$fir_BAprop - inventoryTable$fir_BAprop_init)
inventoryTable$diffBeeDbh = abs(inventoryTable$bee_dbh - inventoryTable$bee_dbh_init)
inventoryTable$diffFirDbh = abs(inventoryTable$fir_dbh - inventoryTable$fir_dbh_init)


inventoryTable[!inventoryTable$created, ]$diffGHA = NA
inventoryTable[!inventoryTable$created, ]$diffBeeBAprop = NA
inventoryTable[!inventoryTable$created, ]$diffBeeBAprop = NA
inventoryTable[!inventoryTable$created, ]$diffBeeDbh = NA
inventoryTable[!inventoryTable$created, ]$diffFirDbh = NA



# ignore some inventories
# when mixity is too low (see conditions)

# modify the name of some inventories (to ignore)
inventoryTable$fileToTag = FALSE
for(i_inv in 1:dim(inventoryTable)[1]){
  inv_line = inventoryTable[i_inv, ]
  
  if(!inv_line$created){
    next
  }
  
  # if too low mixity, dont do monospecific PDG inventories
  condition1 = ( inv_line$bee_BAprop_init <= 0.03 | inv_line$fir_BAprop_init <= 0.03 ) & grepl(inv_line$inventoryType, pattern = "monosp")

  # if too low alternative species, dont do CASTANEA alternative species inventories
  condition2 = inv_line$bee_BAprop_init <= 0.03 & grepl(inv_line$inventoryType, pattern = "CASTANEA_hetre")
  condition3 = inv_line$fir_BAprop_init <= 0.03 & grepl(inv_line$inventoryType, pattern = "CASTANEA_sapin")
  
  if(condition1 | condition2 | condition3){
    inventoryTable$fileToTag[i_inv] = TRUE
  }
}

if(TRUE %in% inventoryTable$fileToTag){
  # move ignored simulations inventories
  ignoredDirectory = paste0(outputInventoryFolder, "0_ignored/")
  if(!dir.exists(ignoredDirectory)){
    dir.create(path = ignoredDirectory)
  }
  
  file_to_tag = paste0(outputInventoryFolder, subset(inventoryTable, fileToTag)$inventoryName)
  file_to_tag_renamed = paste0(ignoredDirectory, subset(inventoryTable, fileToTag)$inventoryName)
  file.rename(file_to_tag, file_to_tag_renamed)
  
  # move plot of ignored simulation inventories
  ignoredImgDirectory = paste0(outputInventoryFolder, "0_ignored/0_plots_image/")
  if(!dir.exists(ignoredImgDirectory)){
    dir.create(path = ignoredImgDirectory)
  }
  
  imgfile_to_tag = paste0(outputInventoryFolder, "0_plots_image/", subset(inventoryTable, fileToTag)$inventoryName, ".jpg")
  imgfile_to_tag_renamed = paste0(ignoredImgDirectory, subset(inventoryTable, fileToTag)$inventoryName, ".jpg")
  imgfile_to_tag = imgfile_to_tag[!grepl(imgfile_to_tag, pattern = "CASTANEA")]
  imgfile_to_tag_renamed = imgfile_to_tag_renamed[!grepl(imgfile_to_tag_renamed, pattern = "CASTANEA")]
  file.rename(imgfile_to_tag, imgfile_to_tag_renamed)
}

# write summary table
fwrite(inventoryTable, file = paste0(outputInventoryFolder, "0_inventoryTable.csv"), sep = ";")

# end script ----





# while(sink.number() > 0){sink()}