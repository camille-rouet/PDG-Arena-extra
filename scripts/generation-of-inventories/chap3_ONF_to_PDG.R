# This work © 2024 by Camille Rouet is licensed under [CC BY-NC 4.0](http://creativecommons.org/licenses/by-nc/4.0/).

# INIT ----
rm(list = ls()) ; gc()

library(tidyverse)
library(openxlsx)
library(data.table)



source("scripts/define_folders.R")
source("scripts/generation-of-inventories/CR_methods_for_inventory_generation.R")



# PARAMETERS ----

outputInventoryFolder = "simulation-data/chap-onf/2024-07-12_ONF/"

outputInConsole = FALSE
# outputInConsole = TRUE # to comment

targetCellWidth = 1

# Put stand to zero slope
forceNoSlope = FALSE

nbAltitudeClasses = 5

initYear = 2011

signifNumber = 4

# List of RU (in mm) used to make inventories
targetRUvec = c(50,100)



# IMPORT ----

# PDG option file
pdgOptionFileFolder = paste0(myFolder, "01_docs-these/E_ecriture_notes/chapitre_3/arbres_virtuels_ONF/")
pdgOptionFileName = "PDG_options_ONF_2024_03_21.txt"
pdgInventoryOptions = readInventoryOptionFile(filePath = paste0(pdgOptionFileFolder, pdgOptionFileName))


# GMAP soil : ventoux and vercors
pathFolder = paste0(myFolder, "01_docs-these/O_Sites_GMAP/sol/")
fileName = "soil_CASTANEA_GMAP_depth1600.csv" # Soil depth will be recomputed to have the desired RU
dGMAP_soil = read_delim(paste0(pathFolder, fileName), delim = ";", show_col_types = FALSE)


# ONF tree table
ONF_fileFolder = paste0(myFolder, "01_docs-these/E_ecriture_notes/chapitre_3/arbres_virtuels_ONF/")
ONF_filePath = paste0(ONF_fileFolder, "arbres_virtuels.xlsx")
treeTableONF = read.xlsx(xlsxFile = ONF_filePath, sheet = 1)


# ONF stand table
ONFstand_fileFolder = paste0(myFolder, "01_docs-these/E_ecriture_notes/chapitre_3/arbres_virtuels_ONF/")
ONFstand_filePath = paste0(ONFstand_fileFolder, "plac_virtuelles.xlsx")
standTableONF = read.xlsx(xlsxFile = ONFstand_filePath, sheet = 4)


# GMAP sites properties
pathFolder = paste0(myFolder, "01_docs-these/N_LIDAR/__gitlab/donnees/")
fileName = "GMAP_sites.csv"
path = paste0(pathFolder, fileName)
GMAP_sites_properties = read_delim(path, delim = ";", show_col_types = FALSE)


# GMAP sites and trees
pathFolder = paste0(myFolder, "01_docs-these/O_Sites_GMAP/GMAP_donnees/arbres_indiv_horsprod/")
fileName = "Data_arbres_placettes_horsprod.csv"
path = paste0(pathFolder, fileName)
dGMAP_horsprod_all = read_delim(path, delim = "\t", show_col_types = FALSE)




# CLEAN DATASETS

# Clean dGMAP_horsprod_all (dead trees not properly marked)
dGMAP_horsprod_all[grepl(dGMAP_horsprod_all$num_arbre, pattern = "mort") | grepl(dGMAP_horsprod_all$num_tige, pattern = "mort"), ]$mort = 1
dGMAP_horsprod_all = subset(dGMAP_horsprod_all, essence %in% c("hetre", "sapin") & site %in% c("vl", "vtx", "bg"))

# demographic parameters are considered equal between species
pdgInventoryOptions$demographic_parameters[2, c(is.na(pdgInventoryOptions$demographic_parameters[2, ])) ] = pdgInventoryOptions$demographic_parameters[1,  c(is.na(pdgInventoryOptions$demographic_parameters[2, ])) ]


dGMAP_soil = computeAndAddRU(dGMAP_soil)


# Add unique_id column
treeTableONF$unique_id = paste0(treeTableONF$ligne, "_" , treeTableONF$id)

# Add basal area column
treeTableONF$basalArea= pi * (treeTableONF$diam/2)**2


thisScriptOptions = tibble(outputInConsole = outputInConsole,
                           outputInventoryFolder = outputInventoryFolder,
                           forceNoSlope = forceNoSlope,
                           nbAltitudeClasses = nbAltitudeClasses,
                           initYear = initYear)




# Make a stand table ----
countColumn = colnames(treeTableONF) == "id"
sumColumn = colnames(treeTableONF) == "basalArea"
meanNumericalColumn = vapply(treeTableONF[1, ], FUN = is.numeric, FUN.VALUE = T, USE.NAMES =  F) & !countColumn & !sumColumn
nonNumericalColumn = !meanNumericalColumn & !countColumn & !sumColumn

standTableONF_bis = tibble(treeTableONF[1:length(unique(treeTableONF$ligne)), ])
colnames(standTableONF_bis)[colnames(standTableONF_bis) == "id"] = "nbTree"
standTableONF_bis[, countColumn] = aggregate(treeTableONF[, countColumn], FUN = length, by = list(stand = treeTableONF$ligne) )[-1]
standTableONF_bis[, sumColumn] = aggregate(treeTableONF[, sumColumn], FUN = sum, by = list(stand = treeTableONF$ligne) )[-1]
standTableONF_bis[, meanNumericalColumn] = aggregate(treeTableONF[, meanNumericalColumn], FUN = function(x) signif(mean(x), 2), by = list(stand = treeTableONF$ligne) )[-1]
standTableONF_bis[, nonNumericalColumn] = aggregate(treeTableONF[, nonNumericalColumn], FUN = function(x) paste0(unique(x), collapse = " "), by = list(stand = treeTableONF$ligne) )[-1]


colnames(standTableONF_bis)[colnames(standTableONF_bis) == "basalArea"] = "sumBasalArea"
standTableONF_bis$standArea = vapply(standTableONF_bis$stade, FUN = function(x) switch(x, "2-PB" = 900, "3-BM" = 1600, "4-GB" = 2500), FUN.VALUE = 0, USE.NAMES = F)
standTableONF_bis$G = standTableONF_bis$sumBasalArea /standTableONF_bis$standArea 

standTableONF_bis$x = NULL
standTableONF_bis$y = NULL

# Résumé du rapport nbTree, G et diam (=stade)
ggplot(standTableONF_bis, aes(x = nbTree, y = G, size = diam, color = stade)) + geom_point()



# Stand properties (all sites) ----

# Stand properties (altitude, slope...)
GMAP_sites_vl_inter_tmp = subset(GMAP_sites_properties, site == "vl" & sous_site == "inter")
GMAP_sites_vl_inter = GMAP_sites_vl_inter_tmp[1, ]
numCols = vapply(GMAP_sites_vl_inter, FUN = is.numeric, FUN.VALUE = T)
GMAP_sites_vl_inter_tmp$exposition[GMAP_sites_vl_inter_tmp$exposition < 180] = GMAP_sites_vl_inter_tmp$exposition[GMAP_sites_vl_inter_tmp$exposition < 180] + 360
GMAP_sites_vl_inter[, numCols] = signif( t(tibble(apply(GMAP_sites_vl_inter_tmp[, numCols], FUN = mean, MARGIN = 2))) , signifNumber)
GMAP_sites_vl_inter$exposition = GMAP_sites_vl_inter$exposition %%360
GMAP_sites_vl_inter$code_site = NULL
GMAP_sites_vl_inter$code_site_2 = NULL
GMAP_sites_vl_inter$code_triplet = NULL
GMAP_sites_vl_inter$peuplement = NULL
GMAP_sites_vl_inter$num_sous_site = NULL
GMAP_sites_vl_inter$date = NULL
GMAP_sites_vl_inter$hourLIDAR = NULL
GMAP_sites_vl_inter$dateLIDAR = NULL
rm(GMAP_sites_vl_inter_tmp)


# Soil texture and RU
soil_vl_inter_tmp = subset(dGMAP_soil, site == "Vercors" & subsite == "Inter")
soil_vl_inter = soil_vl_inter_tmp[1, ]
numCols = vapply(soil_vl_inter, FUN = is.numeric, FUN.VALUE = T)
soil_vl_inter[, numCols] = signif( t(tibble(apply(soil_vl_inter_tmp[, numCols], FUN = mean, MARGIN = 2))), signifNumber)
soil_vl_inter$code_site = NULL
soil_vl_inter$species = NULL
soil_vl_inter$species_code = NULL
soil_vl_inter = computeAndAddRU(soil_vl_inter)
rm(soil_vl_inter_tmp)






# mean ratio hcb/height (vl inter only)
dGMAP_horsprod_vl_inter = subset(dGMAP_horsprod_all, site == "vl" & sous_site == "inter")
dGMAP_horsprod_vl_inter_noNA = dGMAP_horsprod_vl_inter[ !is.na(dGMAP_horsprod_vl_inter$htot) &!is.na(dGMAP_horsprod_vl_inter$h_base_houppier) ,]
dGMAP_horsprod_vl_inter_noNA_beech = subset(dGMAP_horsprod_vl_inter_noNA, essence == "hetre")
dGMAP_horsprod_vl_inter_noNA_fir = subset(dGMAP_horsprod_vl_inter_noNA, essence == "sapin")
summary(dGMAP_horsprod_vl_inter_noNA$h_base_houppier / dGMAP_horsprod_vl_inter_noNA$htot)
summary(dGMAP_horsprod_vl_inter_noNA_beech$h_base_houppier / dGMAP_horsprod_vl_inter_noNA_beech$htot)
summary(dGMAP_horsprod_vl_inter_noNA_fir$h_base_houppier / dGMAP_horsprod_vl_inter_noNA_fir$htot)

ratio_height_hcb = mean(dGMAP_horsprod_vl_inter_noNA$h_base_houppier / dGMAP_horsprod_vl_inter_noNA$htot)
ratio_height_hcb_beech = mean(dGMAP_horsprod_vl_inter_noNA_beech$h_base_houppier / dGMAP_horsprod_vl_inter_noNA_beech$htot)
ratio_height_hcb_fir = mean(dGMAP_horsprod_vl_inter_noNA_fir$h_base_houppier / dGMAP_horsprod_vl_inter_noNA_fir$htot)

# mean ratio hcb/height (all trees)
dGMAP_horsprod_noNA = dGMAP_horsprod_all[ !is.na(dGMAP_horsprod_all$htot) &!is.na(dGMAP_horsprod_all$h_base_houppier) ,]
summary(dGMAP_horsprod_noNA$h_base_houppier / dGMAP_horsprod_noNA$htot)



# LOOP INVENTORIES ----
i_inv = 0
i_inv_tot = length(targetRUvec) * length(standTableONF$ligne)

# LOOP ON RU
for(i_RU in 1:length(targetRUvec)){
  
  targetRU = targetRUvec[i_RU]
  soil_vl_inter_targetRU = soil_vl_inter
  soil_vl_inter_targetRU$soilHeight = signif(soil_vl_inter_targetRU$soilHeight * targetRU / soil_vl_inter_targetRU$RU, signifNumber)
  soil_vl_inter_targetRU = computeAndAddRU(soil_vl_inter_targetRU, signif = signifNumber)
  
  # LOOP ON CODE SITE
  for(i_code_site in (1:length(standTableONF$ligne))){
    
    i_inv = i_inv + 1
    cat(paste0("\nInventory generation: ",i_inv, " / ", i_inv_tot, "\n" ))
    
    
    code_site = standTableONF$ligne[i_code_site]
    
    treeTableONF_code_site = subset(treeTableONF, ligne == code_site)
    standTable_code_site = subset(standTableONF, ligne == code_site)
    standTable_bis_code_site = subset(standTableONF_bis, ligne == code_site)
    LAI = standTable_code_site$PAI_cut
    standArea = standTable_code_site$Surface # in m2
    standSideSize = sqrt(standTable_code_site$Surface) # in m
    
    # Create the output inventory folder
    if(!dir.exists(outputInventoryFolder)){
      dir.create(outputInventoryFolder, recursive = T)
    }
    
    
    
    # Cells table ----
    
    # Generate cell table
    summary(treeTableONF_code_site$x)
    summary(treeTableONF_code_site$y)
    
    nCellSide = round(standSideSize / targetCellWidth)
    cellWidth = standSideSize / nCellSide
    
    cellTable = computeCellInfos(cell_line = soil_vl_inter_targetRU, nCellSide)
    
    res = find_cell_and_rpositon(x = treeTableONF_code_site$x, y = treeTableONF_code_site$y, cells = cellTable, cellWidth = cellWidth)
    treeTableONF_code_site$cellId = res$cellId
    treeTableONF_code_site$res_x = res$res_x
    treeTableONF_code_site$res_y = res$res_x
    
    # Generate height of crown base
    
    # treeTableONF_code_site$hcb_computed = treeTableONF_code_site$htot * ratio_height_hcb
    treeTableONF_code_site$hcb_computed = treeTableONF_code_site$htot * ifelse(treeTableONF_code_site$esp == "Fagus sylvatica", ratio_height_hcb_beech, ratio_height_hcb_fir)
    
    # Generate tree table ----
    treeTable = tibble(idTree = treeTableONF_code_site$id,
                       pop = 1,
                       cell_IDs = treeTableONF_code_site$cellId,
                       sp = pdgInventoryOptions$demographic_parameters$idSp[match(treeTableONF_code_site$esp, pdgInventoryOptions$demographic_parameters$speciesFullNames)],
                       speciesName = treeTableONF_code_site$esp,
                       x_abs = treeTableONF_code_site$x,
                       y_abs = treeTableONF_code_site$y,
                       xnorm = treeTableONF_code_site$res_x,
                       ynorm = treeTableONF_code_site$res_y,
                       dbh = treeTableONF_code_site$diam,
                       height = treeTableONF_code_site$htot,
                       hcb = treeTableONF_code_site$hcb_computed,
                       age = treeTableONF_code_site$age,
                       LAIinit = 0
    )
    
    
    
    # ALTITUDE CLASSES
    referenceAltitude = GMAP_sites_vl_inter$altitude   # Altitude at which x = y = 0 (or, simply, y = 0 when aspect is 0)
    longitude= GMAP_sites_vl_inter$position_e     # longitude in degree for CASTANEA library / longitude en degré pour la librairie CASTANEA
    latitude= GMAP_sites_vl_inter$position_n     # latitude in degree for CASTANEA library / latitude en degré pour la librairie CASTANEA
    
    # if slope_deg = 45, then 1m along y axis correpond to 1m a=of altitudinal variation
    slope_deg = GMAP_sites_vl_inter$pente # Slope of the plot in degree / pente en degré
    
    aspect_deg = GMAP_sites_vl_inter$exposition # exposition of the slope in degree
    
    if(forceNoSlope){
      slope_deg = 0
      aspect_deg = 0
    }
    
    
    pdgInventoryOptions$plotLocationParameters = tibble(referenceAltitude = referenceAltitude, latitude = latitude, longitude = longitude,
                                                        slope_deg = slope_deg, aspect_deg = aspect_deg)
    
    
    
    altitudeClasses = getAltittudeClases(xMin = 0, xMax = standSideSize, yMin = 0, yMax = standSideSize, slope_deg, aspect_deg, referenceAltitude, nbAltitudeClasses)
    
    
    
    # Allele tables
    nbLocus_FCRITBB = pdgInventoryOptions$genetics_parameters$nbLocus_FCRITBB
    nbLocus_g1max = pdgInventoryOptions$genetics_parameters$nbLocus_g1max
    nbMsat = pdgInventoryOptions$genetics_parameters$nbMsat
    nbSNP = pdgInventoryOptions$genetics_parameters$nbSNP
    
    treeSpeciesTable = table(treeTable$sp)
    
    if(!is.na(treeSpeciesTable["1"])){
      allele1_1 = matrix(1 + rbinom(size = 1, n = 30, prob = 0.5), 
                                   nrow = sapply(treeSpeciesTable["1"], FUN = function(x) ifelse(is.na(x), 0 ,x)), # count the number of Beech
                                   ncol = nbLocus_FCRITBB + nbLocus_g1max + nbMsat+nbSNP)
      allele2_1 = matrix(1 + rbinom(size = 1, n = 30, prob = 0.5), 
                         nrow = sapply(treeSpeciesTable["1"], FUN = function(x) ifelse(is.na(x), 0 ,x)), # count the number of Beech
                         ncol = nbLocus_FCRITBB + nbLocus_g1max + nbMsat+nbSNP)
    }else{
      allele1_1 = NA
      allele2_1 = NA
    }
    
    if(!is.na(treeSpeciesTable["2"])){
      allele1_2 = matrix(1 + rbinom(size = 1, n = 30, prob = 0.5), 
                         nrow = sapply( treeSpeciesTable["2"], FUN = function(x) ifelse(is.na(x), 0 ,x)), # count the number of Fir
                         ncol = nbLocus_FCRITBB + nbLocus_g1max + nbMsat+nbSNP)
      allele2_2 = matrix(1 + rbinom(size = 1, n = 30, prob = 0.5), 
                         nrow = sapply(treeSpeciesTable["2"], FUN = function(x) ifelse(is.na(x), 0 ,x)), # count the number of Fir
                         ncol = nbLocus_FCRITBB + nbLocus_g1max + nbMsat+nbSNP)
    }else{
      allele1_2 = NA
      allele2_2 = NA
    }
    
    
    # PDG plot parameters
    pdgInventoryOptions$pdgPlotParameters = tibble(nlin = nCellSide, ncol = nCellSide, cellWidth = cellWidth,
                                                   xorigin = 0, yorigin = 0,
                                                   LAI = LAI)
    
    # PDG genetics variables
    pdgInventoryOptions$genetics_variables = tibble(mean_FCRITBBval = NA, sd_FCRITBBval = NA,
                                                    mean_g1maxVal = NA, sd_g1maxVal = NA)
    
    
    
    # Comment of introduction in inventory
    description = paste0("## Capsis 4.1 - module PDG",
                         "\n# Date and time: ", paste(Sys.time()), 
                         "\n# Inventory generated from stand ", code_site, " made by ONF",
                         "\n# Generated by Camille Rouet 2024")
    
    
    
    # Write inventory ----
    inventoryName = paste0("RU", targetRU, "_", standTableONF_bis$code[standTableONF_bis$ligne == code_site], ".inv")
    writeInventoryPDG(inventoryFolder = outputInventoryFolder, inventoryName = inventoryName, description = description,
                      cellsTable = cellTable,
                      treesTable = treeTable,
                      pdgInventoryOptions = pdgInventoryOptions,
                      altitudeClasses_ = altitudeClasses,
                      allele1_1 = allele1_1,
                      allele2_1 = allele2_1,
                      allele1_2 = allele1_2,
                      allele2_2 = allele2_2,
                      thisScriptOptions = thisScriptOptions, ratioOnRadiusPlot = 0.9)
  } # end loop on codesite
} # end loop on target RU

# Execute to put output back to console
while(sink.number() > 0){sink()}
