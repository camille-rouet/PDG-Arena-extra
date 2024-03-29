# Author: Camille Rouet 2021-2024
# Script for simulation analysis of PDG-Arena and CASTANEA
# Original paper (submitted to PCI Forest and Wood Sciences) : Rouet, Davi, Druel, Fady & Morin (2024): PDG-Arena: An eco-physiological model for characterizing tree-tree interactions in heterogeneous stands 

# Recommandation for external user :
#   1. Go to QUICK CONFIGURATION: load the .RData file & source "CRMethodsForSimulations.R" 
#   2. Go directly to the graph sections (skip CONFIGURATION, IMPORT and CONVERSION sections)

# Vocabulary :
# PDG-Light : the previous name of PDG-Arena
# E0 = RegdemoMonosp = RN inventories
# E1A = RegdemoPlurisp = RN inventories
# E1B = IrregdemoMonosp = irregular and monospecific inventories (not used in the publication)
# E2 = IrregdemoPlurisp = O inventories


# QUICK CONFIGURATION   --------------------------------------------------------
load("simulation-data/paper-PDG-Arena/2023-12-15_goodHeightGoodAlignment.RData")
source("scripts/analysis-of-simulations/CRMethodsForSimulations.R")


# CONFIGURATION  ---------------------------------------------------------------
rm(list = ls()) ; gc() # clear

capsisPath = ""
varPath = paste0(capsisPath, "var/") 
workfilesPath = ""
PROGRAM_ON_SERVER = FALSE

source("scripts/define_folders.R")

source("scripts/analysis-of-simulations/CRMethodsForSimulations.R")



# IMPORT -----------------------------------------------------------------------


# IMPORT GMAP GROWTH DATA 
pathFolder = paste0(workfilesPath, "data_GMAP/")
fileName = "matt_dendro_densite_croiss_interdat_age_tot.csv"
path = paste0(pathFolder, fileName)
dGMAP_dendro = read_delim(path, delim = ";", show_col_types = FALSE)


# IMPORT GMAP INVENTORY DATA 
pathFolder = paste0(workfilesPath, "data_GMAP/")
fileName = "Data_arbres_placettes_horsprod.csv"
path = paste0(pathFolder, fileName)
dGMAP_horsprod = read_delim(path, delim = "\t", show_col_types = FALSE)
# remove dead trees (num_arbre is like "mort")
num_arbre_conversionNumeric = suppressWarnings(as.numeric(dGMAP_horsprod$num_arbre))
dGMAP_horsprod = subset(dGMAP_horsprod, !is.na(num_arbre_conversionNumeric))


# import all GMAP tree informations
dGMAP_dendro$num_index = getNumIndex(dGMAP_dendro$num_arbre, dGMAP_dendro$num_tige)
dGMAP_dendro$treeGlobalId = paste0(dGMAP_dendro$code_site, '_', dGMAP_dendro$num_index)
dGMAP_horsprod$num_index = getNumIndex(as.numeric(dGMAP_horsprod$num_arbre), dGMAP_horsprod$num_tige)
dGMAP_horsprod$treeGlobalId = paste0(dGMAP_horsprod$code_site, '_', dGMAP_horsprod$num_index)

# filter by site
dGMAP_dendro = subset(dGMAP_dendro, site %in% c("vtx", "bg", "vl"))
dGMAP_horsprod = subset(dGMAP_horsprod, site %in% c("vtx", "bg", "vl"))

# GMAP table with all trees (ie even those who were not simulated)
treeYearTableGMAP_all = makeGMAP_TreeYearTable(dGMAP_dendro, dGMAP_horsprod)



# IMPORT SIMULATIONS
fmOutput = 1 # similar to fmSettings.output
logList = logListList[[fmOutput]]
import_irreg_demo_monosp = FALSE
import_CAST = TRUE
keepFilter = ""

currentSimulation = "2023-12-15_goodHeightGoodAlignment/"
simulationFolderGlobal = paste0("2023_simu_article/", currentSimulation)
folderPlot = paste0("plots/", currentSimulation, "bazardeplot/")


# import CASTANEA simulation
if(import_CAST)
  simuListCAST = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "CAST/"),
                                logList = logList, keepFilter = keepFilter)

# import PDG-Arena irregdemo_plurisp simulation (E2 mode)
simuListIrregdemoPlurisp = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "irregdemo_plurisp/"),
                                          logList = logList, keepFilter = keepFilter)

# import preliminary CASTANEA simulation
simuListIrregdemoPlurisp_preli = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "irregdemo_plurisp/"), 
                                                logList = logList, preliminarySimulation = TRUE, keepFilter = keepFilter)


# import PDG-Arena irregdemo_monosp simulation (E1B mode)
if(import_irreg_demo_monosp)
  simuListIrregdemoMonosp = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "irregdemo_monosp/"),
                                           logList = logList)


# import PDG-Arena regdemo_plurisp simulation (E1A mode)
simuListRegdemoPlurisp = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "regdemo_plurisp/"),
                                        logList = logList, keepFilter = keepFilter)

# import PDG-Arena regdemo_monosp simulation (E0 mode)
simuListRegdemoMonosp = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "regdemo_monosp/"),
                                       logList = logList, keepFilter = keepFilter)


# rename some list names to have only code_site as name
# eg : PDGL2_vtx_haut_sp_1b_sapin_PDG_regdemo_monosp to vtx_haut_sp_1b

for(i in 1:length(names(simuListIrregdemoPlurisp))){
  names(simuListIrregdemoPlurisp)[i] = getCodeSiteFromSimulationName(names(simuListIrregdemoPlurisp)[i])
}

for(i in 1:length(names(simuListRegdemoPlurisp))){
  names(simuListRegdemoPlurisp)[i] = getCodeSiteFromSimulationName(names(simuListRegdemoPlurisp)[i])
}

# add Wood Volume Increment
simuListIrregdemoPlurisp = addWVIfromBAIYearlyResultsAndGMAPheight(simuListIrregdemoPlurisp, dGMAP_horsprod)
if(import_irreg_demo_monosp)
  simuListIrregdemoMonosp = addWVIfromBAIYearlyResultsAndGMAPheight(simuListIrregdemoMonosp, dGMAP_horsprod)
simuListRegdemoPlurisp = addWVIfromBAIYearlyResultsAndGMAPheight(simuListRegdemoPlurisp, dGMAP_horsprod, averagePerSpecies = TRUE)
simuListRegdemoMonosp = addWVIfromBAIYearlyResultsAndGMAPheight(simuListRegdemoMonosp, dGMAP_horsprod, averagePerSpecies = TRUE)

if(import_CAST)
  simuListCAST = addWVIfromBAIYearlyResultsAndGMAPheight_standScale(simuListCAST, dGMAP_horsprod)


# Compute stand scale simulations results
simuListIrregdemoPlurisp_st = getStandScaleSimuList(simuListIrregdemoPlurisp)
simuListRegdemoPlurisp_st = getStandScaleSimuList(simuListRegdemoPlurisp)
simuListRegdemoMonosp_st = getStandScaleSimuList(simuListRegdemoMonosp)




# [MONOSP COMBINAISON]
# Combine the results from monospecific simulations based on the same stand

sp_Beech = 3

# proportion of beech for each code_site, in proportion of basal area
proportionBeechBA = getInitialSpeciesBasalAreaProportionInSimulationList(simuListIrregdemoPlurisp, sp_Beech )
# proportion of beech for each code_site, in proportion of tree
proportionBeechN = getInitialSpeciesTreeProportionInSimulationList(simuListIrregdemoPlurisp, sp_Beech )



# MultiCombine PDG-Arena regdemo simulations

# code_site of simulations that have been carried out as monospecific beech
beechOnlyInventories = grepl(names(simuListRegdemoMonosp_st), pattern = "hetre")
# code_site of simulations that have been carried out as monospecific fir
firOnlyInventories = grepl(names(simuListRegdemoMonosp_st), pattern = "sapin")

# selection of inventories that need combination of monospecific simulations
prefixInventories_monosp = vapply(names(simuListRegdemoMonosp_st)[beechOnlyInventories], FUN = getCodeSiteFromSimulationName, FUN.VALUE = "" , USE.NAMES = F)
prefixInventories_plurisp = vapply(names(simuListIrregdemoPlurisp), FUN = function(x) strsplit(x, split = "_PDG")[[1]][1], FUN.VALUE = "" , USE.NAMES = F)
inventoriesDoneInMonosp = prefixInventories_plurisp %in% prefixInventories_monosp

# Combine simulations
combSimuListRegdemoMonosp_st = combineStandScaleSimulationsList(data1 = simuListRegdemoMonosp_st[beechOnlyInventories],
                                                                data2 = simuListRegdemoMonosp_st[firOnlyInventories],
                                                                proportionBA1 = proportionBeechBA[inventoriesDoneInMonosp],
                                                                proportionN1 = proportionBeechN[inventoriesDoneInMonosp], 
                                                                outputDepth = fmOutput)


# MultiCombine CASTANEA simulations
if(import_CAST){
  # code_site of simulations that have been carried out as monospecific beech
  beechOnlyInventories = names(simuListCAST)[grepl(names(simuListCAST), pattern = "hetre")]
  beechOnlyInventories_code_site = sapply(beechOnlyInventories, FUN = getCodeSiteFromSimulationName, USE.NAMES = F)
  
  # code_site of simulations that have been carried out as monospecific fir
  firOnlyInventories = names(simuListCAST)[grepl(names(simuListCAST), pattern = "sapin")]
  firOnlyInventories_code_site = sapply(firOnlyInventories, FUN = getCodeSiteFromSimulationName, USE.NAMES = F)
  
  # code_site of simulations that have been carried out as monospecific beech and monospecific fir
  # Simulations based on these site therefore need to be combined
  beechANDfir_code_site = beechOnlyInventories_code_site[beechOnlyInventories_code_site %in% firOnlyInventories_code_site]
  
  # for each stand that need combination, find in simuListCAST the inventory index corresponding to its monospecific beech simulation
  index_beech_only_andtocombine = sapply(beechANDfir_code_site, 
                                         FUN = function(xcodesite) which( sapply(names(simuListCAST), 
                                                                                 FUN = function(xname) grepl(x = xname, pattern = xcodesite) & grepl(x = xname, pattern = "hetre"), USE.NAMES = F) ))
  
  # for each stand that need combination, find in simuListCAST the inventory index corresponding to its monospecific fir simulation
  index_fir_only_andtocombine = sapply(beechANDfir_code_site, 
                                       FUN = function(xcodesite) which( sapply(names(simuListCAST), 
                                                                               FUN = function(xname) grepl(x = xname, pattern = xcodesite) & grepl(x = xname, pattern = "sapin"), USE.NAMES = F) ))
  
  # combine the simulations
  combSimuListCAST = combineStandScaleSimulationsList(data1 = simuListCAST[index_beech_only_andtocombine],
                                                      data2 = simuListCAST[index_fir_only_andtocombine],
                                                      proportionBA1 = proportionBeechBA[beechANDfir_code_site],
                                                      proportionN1 = proportionBeechN[beechANDfir_code_site], 
                                                      outputDepth = fmOutput)
}


# MultiCombine preliminary CASTANEA simulations

# Recreate a list of preliminary simulation 
simuListIrregdemoPlurisp_preli_bis = list()
for(i in 1:length(simuListIrregdemoPlurisp_preli)){
  aSimuName = names(simuListIrregdemoPlurisp_preli)[i]
  aSimu = simuListIrregdemoPlurisp_preli[[aSimuName]]
  
  # extract the preliminary simulation corresponding to beech and fir respectively
  aSimuHetre = sfiltering(aSimu, fmSpeciesIDs = 3)
  aSimuSapin = sfiltering(aSimu, fmSpeciesIDs = 1)
  
  if(dim(aSimuHetre$yearlyResults)[1] > 0){
    simuListIrregdemoPlurisp_preli_bis[[paste0(aSimuName, "_hetre")]] = aSimuHetre
  }
  
  if(dim(aSimuSapin$yearlyResults)[1] > 0){
    simuListIrregdemoPlurisp_preli_bis[[paste0(aSimuName, "_sapin")]] = aSimuSapin
  }
}

# code_site of simulations that have been carried out as monospecific beech
beechOnlyInventories = names(simuListIrregdemoPlurisp_preli_bis)[grepl(names(simuListIrregdemoPlurisp_preli_bis), pattern = "hetre")]
beechOnlyInventories_code_site = sapply(beechOnlyInventories, FUN = getCodeSiteFromSimulationName, USE.NAMES = F)

# code_site of simulations that have been carried out as monospecific fir
firOnlyInventories = names(simuListIrregdemoPlurisp_preli_bis)[grepl(names(simuListIrregdemoPlurisp_preli_bis), pattern = "sapin")]
firOnlyInventories_code_site = sapply(firOnlyInventories, FUN = getCodeSiteFromSimulationName, USE.NAMES = F)


# code_site of simulations that have been carried out as monospecific beech and monospecific fir
# Simulations based on these site therefore need to be combined
beechANDfir_code_site = beechOnlyInventories_code_site[beechOnlyInventories_code_site %in% firOnlyInventories_code_site]

# for each stand that need combination, find in simuListIrregdemoPlurisp_preli_bis the inventory index corresponding to its monospecific beech simulation
index_beech_only_andtocombine = sapply(beechANDfir_code_site, 
                                       FUN = function(xcodesite) which( sapply(names(simuListIrregdemoPlurisp_preli_bis), 
                                                                               FUN = function(xname) grepl(x = xname, pattern = xcodesite) & grepl(x = xname, pattern = "hetre"), USE.NAMES = F) ))

# for each stand that need combination, find in simuListIrregdemoPlurisp_preli_bis the inventory index corresponding to its monospecific fir simulation
index_fir_only_andtocombine = sapply(beechANDfir_code_site, 
                                     FUN = function(xcodesite) which( sapply(names(simuListIrregdemoPlurisp_preli_bis), 
                                                                             FUN = function(xname) grepl(x = xname, pattern = xcodesite) & grepl(x = xname, pattern = "sapin"), USE.NAMES = F) ))

# combine the simulations
combSimuListPreli = combineStandScaleSimulationsList(data1 = simuListIrregdemoPlurisp_preli_bis[index_beech_only_andtocombine],
                                                     data2 = simuListIrregdemoPlurisp_preli_bis[index_fir_only_andtocombine],
                                                     proportionBA1 = proportionBeechBA[beechANDfir_code_site],
                                                     proportionN1 = proportionBeechN[beechANDfir_code_site])
# end multiCombine Simu preliminaire 



# MultiCombine irregdemo monosp
# Combines individual beeches and firs that were simulated separatly in irregdemo monosp simulation
if(import_irreg_demo_monosp){
  beechOnlyInventories = grepl(names(simuListIrregdemoMonosp), pattern = "hetre")
  firOnlyInventories = grepl(names(simuListIrregdemoMonosp), pattern = "sapin")
  combSimuListIrregdemoMonosp = combineIndividualScaleSimulationList(data1 = simuListIrregdemoMonosp[beechOnlyInventories],
                                                                     data2 = simuListIrregdemoMonosp[firOnlyInventories])
}

# [\COMBINAISON]





# [ASSEMBLE SIMULATION THAT WERE NOT SIMULATED AS MONOSP]
# Some inventories were not simulated in the (regular and irregular) monosp inventory sets because they were 100% monospecific and hence already done as monospecific in the "plurisp" inventory set
# Here, they are reassembled in combSimuList

# Irregdemo monosp PDG-Arena
if(import_irreg_demo_monosp){
  notSimulatedInIrregdemoMonosp = names(simuListIrregdemoPlurisp)[!names(simuListIrregdemoPlurisp) %in% names(combSimuListIrregdemoMonosp)]
  for(inventory in notSimulatedInIrregdemoMonosp){
    combSimuListIrregdemoMonosp[[inventory]] = simuListIrregdemoPlurisp[[inventory]]
  }
  combSimuListIrregdemoMonosp = combSimuListIrregdemoMonosp[sort(names(combSimuListIrregdemoMonosp))]
  combSimuListIrregdemoMonosp_st = getStandScaleSimuList(combSimuListIrregdemoMonosp)
}

# Regdemo monosp PDG-Arena
notSimulatedInRegdemoMonosp = names(simuListRegdemoPlurisp_st)[!names(simuListRegdemoPlurisp_st) %in% names(combSimuListRegdemoMonosp_st)]
if(length(notSimulatedInRegdemoMonosp) > 0){
  for(inventory in notSimulatedInRegdemoMonosp){
    combSimuListRegdemoMonosp_st[[inventory]] = simuListRegdemoPlurisp_st[[inventory]]
  }
  combSimuListRegdemoMonosp_st = combSimuListRegdemoMonosp_st[sort(names(combSimuListRegdemoMonosp_st))]
}


# CASTANEA
if(import_CAST){
  nameSimuListCAST_code_site = sapply(names(simuListCAST), FUN = getCodeSiteFromSimulationName, USE.NAMES = F)
  notSimulatedInCASTANEAMonosp_bool = !nameSimuListCAST_code_site %in% names(combSimuListCAST)
  notSimulatedInCASTANEAMonosp = nameSimuListCAST_code_site[notSimulatedInCASTANEAMonosp_bool]
  
  if(length(notSimulatedInCASTANEAMonosp) > 0 ){
    for(i in 1:length(notSimulatedInCASTANEAMonosp)){
      inventory_codesite = notSimulatedInCASTANEAMonosp[i]
      inventory_simuListCAST = names(simuListCAST)[notSimulatedInCASTANEAMonosp_bool][i]
      combSimuListCAST[[inventory_codesite]] = simuListCAST[[inventory_simuListCAST]]
    }
    combSimuListCAST = combSimuListCAST[sort(names(combSimuListCAST))]
  }
}

# preliminary CASTANEA simulation
nameSimuListPreli_code_site = sapply(names(simuListIrregdemoPlurisp_preli_bis), FUN = getCodeSiteFromSimulationName)
notAssembledInPreliMonosp_bool = !nameSimuListPreli_code_site %in% names(combSimuListPreli)
notAssembledInPreliMonosp = nameSimuListPreli_code_site[notAssembledInPreliMonosp_bool]

if(length(notAssembledInPreliMonosp) > 0 ){
  for(i in 1:length(notAssembledInPreliMonosp)){
    inventory_codesite = notAssembledInPreliMonosp[i]
    inventory_simuListPreli = names(simuListIrregdemoPlurisp_preli_bis)[notAssembledInPreliMonosp_bool][i]
    combSimuListPreli[[inventory_codesite]] = simuListIrregdemoPlurisp_preli_bis[[inventory_simuListPreli]]
  }
  combSimuListPreli = combSimuListPreli[sort(names(combSimuListPreli))]
}



# [\ASSEMBLE MONOSP SIMU LIST]




# [SYNTHESIS OF SIMULATIONS IMPORTATION]

# We should have one simulist containing each 39 simulations for :
# E2 irregdemo plurisp (individual and stand-scale simulations)
length(simuListIrregdemoPlurisp)
length(simuListIrregdemoPlurisp_st)

# E1A regdemo plurisp (individual and stand-scale)
length(simuListRegdemoPlurisp)
length(simuListRegdemoPlurisp_st)

# E0 regdemo monosp recombined (stand-scale)
length(combSimuListRegdemoMonosp_st)

# E1B irregdemo monosp recombined (individual and stand-scale)
if(import_irreg_demo_monosp){
  length(combSimuListIrregdemoMonosp)
  length(combSimuListIrregdemoMonosp_st)
}

# CASTANEA (stand-scale)
if(import_CAST){
  length(combSimuListCAST) 
}

# preliminary CASTANEA simulation (stand-scale)
length(combSimuListPreli) 


# [\SYNTHESIS]




# find the simulated sites and filter GMAP DATA
siteList = NULL
for(sitePotential in c("vl_", "vtx_", "bg_")){
  if(sum(grepl(names(simuListIrregdemoPlurisp), pattern = sitePotential)) > 0){
    
    siteCode = substr(sitePotential, 1, nchar(sitePotential) - 1)
    siteList = c(siteList, siteCode)
  }
}

dGMAP_dendro = subset(dGMAP_dendro, site %in% siteList)
dGMAP_horsprod = subset(dGMAP_horsprod, site %in% siteList)


# [\IMPORT]









# CONVERSION INTO TABLES ------------------------------
# Convert simulation lists into TREE-YEAR tables (for GMAP, E1B and E2) and STAND-YEAR tables (for GMAP, CASTANEA, E0, E1A, E1B and E2)
# TREE-PERIOD and STAND-PERIOD tables are also computed (period are range of years)

# E2 mode
# i. make the primary tree-year table
treeYearTable_E2 = makeTreeYearTable(simuListIrregdemoPlurisp)


treeid_simulated = unique(treeYearTable_E2$treeGlobalId)
treeid_measured = unique(dGMAP_dendro$treeGlobalId)
treeid_simulated_notmeasured = treeid_simulated[!treeid_simulated %in% treeid_measured]
treeid_measured_notsimulated = treeid_measured[!treeid_measured %in% treeid_simulated]
warning(paste0("These trees were measured but not simulated:\n", paste0(treeid_measured_notsimulated, collapse = "\n")))

treeYearTableGMAP = makeGMAP_TreeYearTable(dGMAP_dendro, subset(dGMAP_horsprod, treeGlobalId %in% treeid_simulated))
treeYearTableGMAP = extrapoleMissingMeasuredBAI(treeYearTableGMAP)
treeYearTable_E2 = left_join(treeYearTable_E2, treeYearTableGMAP)


# Check : ALL should be FALSE
table(treeYearTable_E2$BAI_mes == -999)
table(is.na(treeYearTable_E2$BAI_mes))

# ii. make stand-year table
standYearTable_E2 = makeStandYearTable(treeYearTable_E2)

# iii. make tree-period and stand-period tables

periodList = list(c(1996, 2002), c(2003, 2006), c(2007, 2013), c(1996, 2013))
treePeriodTable_E2 = makeTreePeriodTable(treeYearTable = treeYearTable_E2, periodList)
standPeriodTable_E2 = makeStandPeriodTable(standYearTable_E2, periodList)
# end E2 mode

# standYearTable_GMAP (based on standYearTable_E2)
standYearTable_GMAP = standYearTable_E2[, colnames(standYearTable_E2)[! colnames(standYearTable_E2) %in% c("GPPabs_sim", "BAI_sim", "WVI_sim")]]
# treeYearTableGMAP



# E1B mode (irregdemo monosp recombiné)
# i. make the primary tree-year table
if(import_irreg_demo_monosp){
  treeYearTable_E1B = makeTreeYearTable(combSimuListIrregdemoMonosp)
  
  
  treeid_simulated = unique(treeYearTable_E1B$treeGlobalId)
  treeid_measured = unique(dGMAP_dendro$treeGlobalId)
  treeid_simulated_notmeasured = treeid_simulated[!treeid_simulated %in% treeid_measured]
  treeid_measured_notsimulated = treeid_measured[!treeid_measured %in% treeid_simulated]
  
  warning(paste0("These trees were measured but not simulated:\n", paste0(treeid_measured_notsimulated, collapse = "\n")))
  
  # treeYearTableGMAP = makeGMAP_TreeYearTable(dGMAP_dendro, subset(dGMAP_horsprod, treeGlobalId %in% treeid_simulated))
  # treeYearTableGMAP = extrapoleMissingMeasuredBAI(treeYearTableGMAP)
  treeYearTable_E1B = left_join(treeYearTable_E1B, treeYearTableGMAP)
  
  
  # Check : ALL should be FALSE
  table(treeYearTable_E1B$BAI_mes == -999)
  table(is.na(treeYearTable_E1B$BAI_mes))
  
  # ii. make stand-year table
  standYearTable_E1B = makeStandYearTable(treeYearTable_E1B)
  
  # iii. make tree-period and stand-period tables
  
  periodList = list(c(1996, 2002), c(2003, 2006), c(2007, 2013), c(1996, 2013))
  treePeriodTable_E1B = makeTreePeriodTable(treeYearTable = treeYearTable_E1B, periodList)
  standPeriodTable_E1B = makeStandPeriodTable(standYearTable_E1B, periodList)
}
# end E1B mode


# add stand area to standYearTable_GMAP
standYearTable_GMAP$standArea_m2 = -1
for(a_code_site in unique(standYearTable_GMAP$code_site)){
  standYearTable_GMAP[standYearTable_GMAP$code_site == a_code_site, ]$standArea_m2 = simuListIrregdemoPlurisp[[a_code_site]]$inventory$standArea_m2
}


# STAND-YEAR tables from stand scale simulation
standYearTable_E0 = makeStandYearTableFromStandScale(combSimuListRegdemoMonosp_st, standYearTable_GMAP)
if(import_CAST)
  standYearTable_CAST = makeStandYearTableFromStandScale(combSimuListCAST, standYearTable_GMAP)
standYearTable_E1A = makeStandYearTableFromStandScale(simuListRegdemoPlurisp_st, standYearTable_GMAP)
standYearTable_E0PRELI = makeStandYearTableFromStandScale(combSimuListPreli, standYearTable_GMAP)

# # CASTANEA "bis" case
# standYearTable_E0PRELIbis = makeCASTANEAstandYearTable(simuListIrregdemoPlurisp_init, simuListIrregdemoPlurisp_preli)
# standYearTable_E0PRELIbis$code_site = getCodeSiteFromSimulationName(standYearTable_E0PRELIbis$code_site)
# standYearTable_E0PRELIbis$BAI_sim = standYearTable_E0PRELIbis$BAI_simCAST
# standYearTable_E0PRELIbis$GPPabs_sim = standYearTable_E0PRELIbis$GPPabs_simCAST
# standYearTable_E0PRELIbis = inner_join(standYearTable_E0PRELIbis, standYearTable_GMAP)


# STAND-PERIOD tables from stand scale simulation
standPeriodTable_E0 = makeStandPeriodTable(standYearTable_E0, periodList)
if(import_CAST)
  standPeriodTable_CAST = makeStandPeriodTable(standYearTable_CAST, periodList)
standPeriodTable_E1A = makeStandPeriodTable(standYearTable_E1A, periodList)
standPeriodTable_E0PRELI = makeStandPeriodTable(standYearTable_E0PRELI, periodList)
# standPeriodTable_E0PRELIbis = makeStandPeriodTable(standYearTable_E0PRELIbis, periodList)

standPeriodTable_GMAP = standPeriodTable_E2
standPeriodTable_GMAP$BAIy_sim = standPeriodTable_GMAP$BAIy_mes
standPeriodTable_GMAP$GPPy_abs_sim = 0


# E2 and E1B modes, computed directly from stand scale simulation lists
standYearTable_E2_fromStand = makeStandYearTableFromStandScale(simuListIrregdemoPlurisp_st, standYearTable_GMAP)
standPeriodTable_E2_fromStand = makeStandPeriodTable(standYearTable_E2_fromStand, periodList)
standPeriodTable_E2_from_indiv = standPeriodTable_E2
standPeriodTable_E2 = standPeriodTable_E2_fromStand

if(import_irreg_demo_monosp){
  standYearTable_E1B_fromStand = makeStandYearTableFromStandScale(combSimuListIrregdemoMonosp_st, standYearTable_GMAP)
  standPeriodTable_E1B_fromStand = makeStandPeriodTable(standYearTable_E1B_fromStand, periodList)
}



# remove all units of a table
removeUnits = function(a_standPeriodTable){
  for(col in names(a_standPeriodTable)[sapply(a_standPeriodTable[1, ], is.numeric)]){
    a_standPeriodTable[[col]] = as.vector(a_standPeriodTable[[col]])
  }
  return(a_standPeriodTable)
}

# add some columns for a table
addNewColumns = function(a_standPeriodTable, simuList_st){
  a_standPeriodTable$code_triplet_cut = sapply(a_standPeriodTable$triplet, FUN = function(x){splitted = strsplit(x, split = "_")[[1]] ; return(paste0(splitted[-1], collapse =  "."))})
  a_standPeriodTable$code_site_cut = sapply(a_standPeriodTable$code_site, FUN = function(x){splitted = strsplit(x, split = "_")[[1]] ; return(paste0(splitted[-1], collapse =  "."))})
  a_standPeriodTable$LAI = 0
  a_standPeriodTable$GhaInit = 0
  a_standPeriodTable$incident_yearlyMJm2 = 0
  a_standPeriodTable$standArea_m2 = 0
  
  for(a_code_site in unique(a_standPeriodTable$code_site)){
    simu_st = simuList_st[[a_code_site]]
    simu_st_yR = simu_st$yearlyResults
    simu_st_yR_firstYear = simu_st_yR[simu_st_yR$year == min(simu_st_yR$year), ]
    a_standPeriodTable[a_standPeriodTable$code_site == a_code_site, ]$LAI = unique(simu_st_yR$LAImaxThisYear)
    a_standPeriodTable[a_standPeriodTable$code_site == a_code_site, ]$GhaInit = simu_st_yR_firstYear$Gha
    a_standPeriodTable[a_standPeriodTable$code_site == a_code_site, ]$incident_yearlyMJm2 = mean(simu_st_yR_firstYear$incident_yearlyMJm2)
    a_standPeriodTable[a_standPeriodTable$code_site == a_code_site, ]$standArea_m2 = simu_st$inventory$standArea_m2
    
  }
  
  # for(col in names(a_standPeriodTable)[sapply(a_standPeriodTable[1, ], is.numeric)]){
  #   a_standPeriodTable[[col]] = as.vector(a_standPeriodTable[[col]])
  # }
  
  # set units
  for(col in c("BAIy_sim", "BAIy_mes")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], cm2/yr)
  }
  
  # set units
  for(col in c("WVIy_sim", "WVIy_mes")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], m3/yr)
  }
  
  # set units
  for(col in c("GPPy_abs_sim", "NPPy_abs_sim", "Rautoy_abs_sim")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], g/yr)
  }
  
  # set units
  for(col in c("standArea_m2")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], m2)
  }
  
  a_standPeriodTable$WVIy_m2_mes = a_standPeriodTable$WVIy_mes / a_standPeriodTable$standArea_m2
  a_standPeriodTable$WVIy_m2_sim = a_standPeriodTable$WVIy_sim / a_standPeriodTable$standArea_m2
  a_standPeriodTable$GPPy_m2_sim = a_standPeriodTable$GPPy_abs_sim / a_standPeriodTable$standArea_m2
  a_standPeriodTable$NPPy_m2_sim = a_standPeriodTable$NPPy_abs_sim / a_standPeriodTable$standArea_m2
  a_standPeriodTable$BAIy_m2_mes = a_standPeriodTable$BAIy_mes / a_standPeriodTable$standArea_m2
  a_standPeriodTable$BAIy_m2_sim = a_standPeriodTable$BAIy_sim / a_standPeriodTable$standArea_m2
  
  # set units
  for(col in c("GPPy_m2_sim", "NPPy_m2_sim")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], g/m2/yr) # conversion
  }
  
  
  # set new units
  for(col in c("WVIy_m2_mes", "WVIy_m2_sim")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], m3/m2/yr)
  }
  
  # set new units
  for(col in c("BAIy_m2_mes", "BAIy_m2_sim")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], cm2/m2/yr)
  }
  
  return(a_standPeriodTable)
}

# convert the unit for some colymn of a table
changeUnit_TableList = function(a_standPeriodTable, simuList_st){
  
  for(col in c("GPPy_abs_sim", "NPPy_abs_sim", "Rautoy_abs_sim")){
    a_standPeriodTable[[col]] = set_units(a_standPeriodTable[[col]], kg/yr)
  }
  
  return(a_standPeriodTable)
}


# Add columns for stand-period table and add units
standPeriodTable_E2 = addNewColumns(standPeriodTable_E2, simuListIrregdemoPlurisp_st)
standPeriodTable_E2_fromStand = addNewColumns(standPeriodTable_E2_fromStand, simuListIrregdemoPlurisp_st)
if(import_irreg_demo_monosp){
  standPeriodTable_E1B = addNewColumns(standPeriodTable_E1B, simuListIrregdemoPlurisp_st)
  standYearTable_E1B_fromStand = addNewColumns(standYearTable_E1B_fromStand, simuListIrregdemoPlurisp_st)
}
standPeriodTable_E1A = addNewColumns(standPeriodTable_E1A, simuListIrregdemoPlurisp_st)
standPeriodTable_E0 = addNewColumns(standPeriodTable_E0, simuListIrregdemoPlurisp_st)
if(import_CAST)
  standPeriodTable_CAST = addNewColumns(standPeriodTable_CAST, simuListIrregdemoPlurisp_st)
standPeriodTable_E0PRELI = addNewColumns(standPeriodTable_E0PRELI, simuListIrregdemoPlurisp_st)
standPeriodTable_GMAP = addNewColumns(standPeriodTable_GMAP, simuListIrregdemoPlurisp_st)



# CHANGE UNIT /!\ DO NOT RE-EXECUTE THESE LINES, this would change the values
standPeriodTable_E2 = changeUnit_TableList(standPeriodTable_E2, simuListIrregdemoPlurisp_st)
standPeriodTable_E2_fromStand = changeUnit_TableList(standPeriodTable_E2_fromStand, simuListIrregdemoPlurisp_st)
if(import_irreg_demo_monosp){
  standPeriodTable_E1B = changeUnit_TableList(standPeriodTable_E1B, simuListIrregdemoPlurisp_st)
  standYearTable_E1B_fromStand = changeUnit_TableList(standYearTable_E1B_fromStand, simuListIrregdemoPlurisp_st)
}
standPeriodTable_E1A = changeUnit_TableList(standPeriodTable_E1A, simuListIrregdemoPlurisp_st)
standPeriodTable_E0 = changeUnit_TableList(standPeriodTable_E0, simuListIrregdemoPlurisp_st)
if(import_CAST)
  standPeriodTable_CAST = changeUnit_TableList(standPeriodTable_CAST, simuListIrregdemoPlurisp_st)
standPeriodTable_E0PRELI = changeUnit_TableList(standPeriodTable_E0PRELI, simuListIrregdemoPlurisp_st)
standPeriodTable_GMAP = changeUnit_TableList(standPeriodTable_GMAP, simuListIrregdemoPlurisp_st)







suffix_list = c("E2", "E1A", "E0"
                # , "E0PRELI"
                # , "E0PRELIbis"
)

title_list = list(ggtitle("st_per_E2"), ggtitle("st_per_E1A"), 
                  ggtitle("st_per_E0")
                  #, ggtitle("st_per_E0PRELI")
                  # , ggtitle("st_per_E0PRELIbis")
)
standPeriodTable_list = list(standPeriodTable_E2, standPeriodTable_E1A, 
                             standPeriodTable_E0
                             #, standPeriodTable_E0PRELI
                             # , standPeriodTable_E0PRELIbis
)

if(import_irreg_demo_monosp){
  title_list$st_per_E1B = ggtitle("st_per_E1B")
  standPeriodTable_list$st_per_E1B = standPeriodTable_E1B
  suffix_list = c(suffix_list, "E1B")
}

if(import_CAST){
  title_list$st_per_CAST = ggtitle("st_per_CAST")
  standPeriodTable_list$st_per_CAST = standPeriodTable_CAST
  suffix_list = c(suffix_list, "CAST")
}




# [\CONVERSION]

# SYNTHESE
# We have table with simulated and measured Wood Volume Increment for :
# tree-year : E2
# tree-year : E1B (irregdemo monosp)
# stand-year : E2
# stand-year : E1B
# stand-year : E1A
# stand-year : E0
# stand-year : CASTANEA
# + the same tables using periods







# GRAPH SHOWING ALL STANDS --------------------------------------------------------------------
# (Appendix) Graph with all stands 

size = 4
xValMin = 0 ; xValMax = 4.5e-3 ; yValMin = 0 ; yValMax = 4.5e-3 ; coord_limits = coord_cartesian(xlim = c(xValMin,xValMax), ylim = c(yValMin,yValMax))
folderPlot = paste0("plots/", currentSimulation, "bazardeplot/")
for(i in 1:length(standPeriodTable_list)){
  standPeriodTable = standPeriodTable_list[[i]]
  atitle = title_list[[i]]
  ggplot(subset(standPeriodTable, period == "1996_2013"), 
         aes(x = WVIy_m2_mes, y = WVIy_m2_sim, color = composition, label = code_site_cut, shape = factor(site))
  ) + 
    # atitle +
    # facet_grid( . ~ composition) +
    geom_text_repel(size = 3.5, alpha = 0.7) +
    geom_point(size = size) + geom_abline(slope = 1, alpha = 0.25) + 
    coord_limits+
    ylab("WVI simulated")+
    xlab("WVI measured")+
    scale_color_discrete(name = "Composition",
                         labels = c("Mixed", "Beech", "Fir"))+
    scale_shape_manual(name = "Site",
                       labels = c("Bauges", "Vercors", "Ventoux"),
                       values = c(15,19,17)) +
    annotate("text", x=1.3125, y=1.375, label= "1:1", size = 6)
  
  saveLastGgPlot(folderPlot, plot_width = 1280, ratio = 1.1, fileName = paste0("WVIm_m2", "_", atitle$title, ".pdf"))
}



# other graphs
ggplot(subset(standPeriodTable_E2, period == "1996_2013"), aes(x = BAIy_mes, y = BAIy_sim, shape = factor(triplet), color = site)) + scale_shape_manual(values = 0:100) + geom_point(size = size) + geom_abline(slope = 1, alpha = 0.25) + title_4 + coord_limits

ggplot(subset(standPeriodTable_E2, period == "1996_2013"), aes(x = BAIy_mes, y = BAIy_sim, shape = factor(triplet), color = site)) + scale_shape_manual(values = 0:100) + geom_point(size = size) + geom_abline(slope = 1, alpha = 0.25) + title_4 + facet_grid(. ~composition) + coord_limits
ggplot(subset(standPeriodTable_E2, period == "1996_2013"), aes(x = BAIy_mes, y = BAIy_sim, shape = factor(triplet), color = composition)) + scale_shape_manual(values = 0:100) + geom_point(size = size) + geom_abline(slope = 1, alpha = 0.25) + title_4 + facet_grid(. ~site) + coord_limits

ggplot(standPeriodTable_E2, aes(x = BAIy_mes, y = BAIy_sim)) + geom_point() + geom_abline(slope = 1, alpha = 0.25) + facet_grid(composition ~ period) + title_4 + coord_limits
ggplot(standPeriodTable_E2, aes(x = BAIy_mes_scaled, y = GPPy_abs_sim_scaled)) + geom_point()+ geom_abline(slope = 1, alpha = 0.25) + facet_grid(. ~ composition) + title_4 




# GRAPHE GPP
ggplot(subset(standPeriodTable_E2, period == "1996_2013"), 
       aes(x = BAIy_mes, y = BAIy_sim, color = composition, label = code_site_cut, shape = factor(site))
) + 
  atitle + 
  facet_grid(composition ~ site) +
  geom_text_repel(size = 3.5, alpha = 0.7) +
  geom_point(size = size) + geom_abline(slope = 1, alpha = 0.25) + 
  coord_limits + scale_shape_manual(values = c(15,19,17)) 

folderPlot = paste0("plots/", currentSimulation, "resultats/")
saveLastGgPlot(folderPlot, plot_width = 1440, ratio = 4/3, fileName = paste0("bysitecompo_", atitle$title, ".png"))














# 3.1 CORRELATION MATRIX -----------------------------------------------------------
# (Part 3.1) Correlations between modelling situations


listOfTables = list(subset(standPeriodTable_CAST, period == "1996_2013"),
                    subset(standPeriodTable_E0, period == "1996_2013"),
                    subset(standPeriodTable_E1A, period == "1996_2013"),
                    # subset(standPeriodTable_E1B, period == "1996_2013"),
                    subset(standPeriodTable_E2, period == "1996_2013"))

# Initialize and fill a correlation matrix
correlationSimulationMatrix = as.data.frame(tibble( CASTANEA = 0, E0 = 0, E1A = 0, 
                                                    # E1B = 0, 
                                                    E2 = 0, .rows = length(listOfTables)))
colnames(correlationSimulationMatrix) = c("CASTANEA (RN)", "PDG-Arena (RN)", "PDG-Arena (RS)", "PDG-Arena (O)")
rownames(correlationSimulationMatrix) = names(correlationSimulationMatrix)
correlationMode = "correlation" # RMSE 1-r2 correlation
variableOfInterest = "GPPy_m2_sim"

folderPlot = paste0("plots/", currentSimulation, "bazardeplot/")

nI = dim(correlationSimulationMatrix)[1]
nJ = dim(correlationSimulationMatrix)[2]
# for every simulation mode couple
for(i in 1:nI){
  for(j in 1:nJ){
    
    # extract simulation results
    simulationTable1 = listOfTables[[i]]
    simulationTable2 = listOfTables[[j]]
    
    # isolate variable of interest
    variableOfInterestRename1 = paste0(variableOfInterest, "1")
    variableOfInterestRename2 = paste0(variableOfInterest, "2")
    simulationTable1[[variableOfInterestRename1]] = simulationTable1[[variableOfInterest]]
    simulationTable2[[variableOfInterestRename2]] = simulationTable2[[variableOfInterest]]
    
    # created a coupled table
    coupledTable = inner_join(simulationTable1, simulationTable2, by = c("code_site", "site", "triplet", "composition", "period", 
                                                                         "code_triplet_cut", "code_site_cut") )
    
    # get correlation for the variable of interest
    correlationOfCouple = getComparisonCoefficient(coupledTable, nameVar1 = variableOfInterestRename1, nameVar2 = variableOfInterestRename2)
    
    correlationSimulationMatrix[i, j] = correlationOfCouple[correlationOfCouple$test == correlationMode,]$PDGvsMES
    
    
    # # graphes 2 à 2
    # a_title = paste0(names(correlationSimulationMatrix)[i], "_v_", names(correlationSimulationMatrix)[j])
    # ggplot(coupledTable,
    #        aes(x = BAIy_sim1, y = BAIy_sim2, color = composition, label = code_site_cut, shape = factor(site))
    # ) +
    #   ggtitle(a_title) +
    #   # facet_grid(composition ~ site) +
    #   geom_text_repel(size = 3.5, alpha = 0.7) +
    #   geom_point(size = size) + geom_abline(slope = 1, alpha = 0.25) +
    #   coord_limits + scale_shape_manual(values = c(15,19,17))
    # saveLastGgPlot(folderPlot, plot_width = 1440, ratio = 4/3, fileName = paste0(a_title, ".png"))
    
  }
}


# Plot the correlation matrix
textSize = 15

myMat = correlationSimulationMatrix ; matrixValues = myMat[lower.tri(myMat)]
minMatrixValue = floor(min(matrixValues)*100)/100 ; maxMatrixValue = ceiling(max(matrixValues)*1000)/1000 ; sdMatrix = sd(matrixValues) ; meanMatrix = mean(matrixValues)
paste0(correlationMode, " for ", variableOfInterest)
ggcorrplot(myMat, legend.title = correlationMode, 
           lab = T, lab_size = 6,
           show.diag = T,
           digits = 3,
           show.legend = F,
           method = "square", type = "upper",
           outline.color = "black") +
  # ggtitle(paste0(correlationMode, " for ", variableOfInterest)) + 
  scale_fill_gradient(name = "r^2", 
                      breaks = round(c(min(matrixValues), mean(range(matrixValues)), max(matrixValues)),2),
                      limit = c(minMatrixValue - 0.5 * sdMatrix, maxMatrixValue + 0.5*sdMatrix), 
                      low = hsv(0, 0, 1), high = hsv(0.65, 0.8, 0.7) )+
  theme(axis.text.x = element_text(size = textSize, 
                                   angle=0, hjust = 0.5, vjust = 0.75), # horizontal and centered bottom text
        axis.text.y = element_text(size = textSize), # horizontal and centered bottom text
        # legend.key.size = unit(1, 'cm'), # legend size
        text=element_text(size=textSize) # text size
  ) 
folderPlot = paste0("plots/", currentSimulation, "bazardeplot/")
saveLastGgPlot(folderPlot, plot_height = 720, ratio = 4/3, fileSuffix = ".pdf")










# 3.1 COMPARISON BETWEEN SIMULATIONS  ----------------------------------------
# (Part 3.1) Graph of a variable from one simulation vs another


# Table simulation A vs simulation B
suffix1 = "_E0"
suffix2 = "_CAST"
simulationComparaisonTable = get(paste0("standPeriodTable", suffix1))
variableOfInterest = "GPPy_m2_sim"
variableList = c("BAIy_sim", "BAIy_m2_sim", "BAIy_m2_sim", "WVIy_sim", "WVIy_m2_mes", "WVIy_m2_sim", "GPPy_abs_sim", "GPPy_m2_sim", "NPPy_abs_sim", "NPPy_m2_sim", "Rautoy_abs_sim", 
                 "vegAbso_MJm2Y", "vegAbsorbance", "vegAbsoPAR_MJm2Y", "vegAbsorbancePAR", 
                 "REWmin", "RU_level_min", "RU_shortage_max")

for(var in variableList){
  colnames(simulationComparaisonTable)[colnames(simulationComparaisonTable) == var] = paste0(var, suffix1)
}
simulationComparaisonTable = inner_join(simulationComparaisonTable, get(paste0("standPeriodTable", suffix2))) # rename all new columns if they are added

for(var in variableList){
  colnames(simulationComparaisonTable)[colnames(simulationComparaisonTable) == var] = paste0(var, suffix2)
}


# Plot 

pointSize = 3 ; textSize = 14
xValMin = 0 ; xValMax = 1780 ; yValMin = 0 ; yValMax = 1780 ; coord_limits = coord_cartesian(xlim = c(xValMin,xValMax), ylim = c(yValMin,yValMax))

simulationComparaisonTable_sub = subset(simulationComparaisonTable, period == "1996_2013")
simulationComparaisonTable_sub = removeUnits(simulationComparaisonTable_sub)
ggplot(simulationComparaisonTable_sub, 
       aes_string(x = paste0(variableOfInterest, suffix2), y = paste0(variableOfInterest, suffix1), color = "composition", shape = "factor(site)"
                  # label = code_site_cut, shape = factor(site)
       )
) + 
  facet_grid( . ~ .) +
  # geom_text_repel(size = 3.5, alpha = 0.7) +
  geom_point(size = pointSize, alpha = 0.75) + geom_abline(slope = 1, alpha = 0.25) + 
  coord_limits + 
  scale_shape_manual(name = "Site",
                     labels = c("Bauges", "Vercors", "Ventoux"),
                     values = c(15,19,17))+
  # ggtitle("PDG-Arena (R+N) vs CASTANEA (R+N) ") + 
  scale_color_discrete(name = "Composition",
                       labels = c("Mixed", "Beech", "Fir"))+ 
  xlab("GPP simulated with CASTANEA [g/m2/yr]")+
  ylab("GPP simulated with PDG-Arena [g/m2/yr]")+
  theme(
    # text=element_text(size=textSize), #change font size of all text
    # axis.text=element_text(size=textSize), #change font size of axis text
    axis.title=element_text(size=textSize), #change font size of axis titles
    # plot.title=element_text(size=textSize), #change font size of plot title
    legend.text=element_text(size=textSize), #change font size of legend text
    legend.title=element_text(size=textSize) #change font size of legend title
  )+
  annotate("text", x=1700, y= 1600, label= "1:1", size = 4.5) +
  annotate("text", x=920, y= 1250, label= paste0("r = ", round(cor(simulationComparaisonTable_sub[[paste0(variableOfInterest, suffix2)]], simulationComparaisonTable_sub[[paste0(variableOfInterest, suffix1)]]), 3)
  ) , size = 4.5) 
folderPlot = paste0("plots/", currentSimulation, "bazardeplot/")
saveLastGgPlot(folderPlot, plot_width = 720, ratio = 1.20, fileSuffix = ".pdf")













# 3.2 STATISTICS -------------------------------------------------------------------
# (Part 3.2) Correlation and other coefficient on simulated versus measured variables

# Defines the variables of interest and coefficient to comptute
var2 = "WVIy_m2_sim"
var1 = "WVIy_m2_mes"
coefficients_list = c("r2", "MAPE")
# coefficients_list = c("correlation", "r2", "1-r2", "RMSE", "MAPE")

# Compute the coefficients
stat_list = list()
index = 1
for(standPeriodTable in standPeriodTable_list){
  stat_list[[suffix_list[index]]] = getComparisonCoefficientPerSiteAndComposition(subset(standPeriodTable, period == "1996_2013"), var1, var2, coefficients_list)
  index = index + 1
}

# Show the coefficients
printList = list("STAT GLOBALES", 
                 "Effet structure : E1A vs E2")
for(stat_item_name in names(stat_list)){
  printList[[stat_item_name]] = stat_list[[stat_item_name]][["global"]]
}
print(printList)

printList = list("STAT MELANGE",
                 "Effet diversité : E0 vs E1A", 
                 "Effet structure : E1A vs E2")
for(stat_item_name in names(stat_list)){
  printList[[stat_item_name]] = stat_list[[stat_item_name]][["comp_m"]]
}
print(printList)

printList = list("STAT HETRES", 
                 "Effet structure : E1A vs E2")
for(stat_item_name in names(stat_list)){
  printList[[stat_item_name]] = stat_list[[stat_item_name]][["comp_ph"]]
}
print(printList)

printList = list("STAT SAPINS", 
                 "Effet structure : E1A vs E2")
for(stat_item_name in names(stat_list)){
  printList[[stat_item_name]] = stat_list[[stat_item_name]][["comp_sp"]]
}
print(printList)













# 3.3 BOXPLOT OF VARIABLES ---------------------------------------------------------
# (Part 3.3) Comparison of simulated variables between modelling situations

standPeriodTable_E2 = standPeriodTable_E2_fromStand

# Filter on a composition set
# all: "", mixed: "m", pure beech: "ph", pure fir: "sp"
filter_composition = "m"

# Table simulation A vs simulation B
variableList = colnames(standPeriodTable_E0)[! colnames(standPeriodTable_E0) %in% c("code_site", "site", "triplet", "composition", "period", "code_triplet_cut", "code_site_cut", "LAI", "GhaInit", "standArea_m2")]


# create a table with one line per code_site and one column per variable per mode
lollypopTable = standPeriodTable_E0
for(var in variableList){
  colnames(lollypopTable)[colnames(lollypopTable) == var] = paste0(var, "_E0")
}

lollypopTable = inner_join(lollypopTable, standPeriodTable_E1A) # rename all new columns if they are added
for(var in variableList){
  colnames(lollypopTable)[colnames(lollypopTable) == var] = paste0(var, "_E1A")
}

# ADD E2 and CASTANEA
lollypopTable = inner_join(lollypopTable, standPeriodTable_E2) # rename all new columns if they are added
for(var in variableList){
  colnames(lollypopTable)[colnames(lollypopTable) == var] = paste0(var, "_E2")
}

lollypopTable = inner_join(lollypopTable, standPeriodTable_CAST) # rename all new columns if they are added
for(var in variableList){
  colnames(lollypopTable)[colnames(lollypopTable) == var] = paste0(var, "_CAST")
}


# filter on period and composition
lollypopTable = subset(lollypopTable, period == "1996_2013")

if(filter_composition %in% c("m", "ph", "sp")){
  lollypopTable = subset(lollypopTable, composition == filter_composition)
}


# code_site renaming : from bg_bas_m_4, bg_bas_m_5.. to Bauges 1, Bauges 2 etc.
new_code_site = sapply(unique(lollypopTable$code_site), FUN = getSiteFromCodeSite, USE.NAMES = F)
new_code_site = ifelse(new_code_site == "bg", "Bauges", ifelse(new_code_site == "vtx", "Ventoux", ifelse(new_code_site == "vl", "Vercors", "")))
new_code_site = paste0(new_code_site, " ", ave(new_code_site == new_code_site, new_code_site, FUN = cumsum))
lollypopTable$code_site_old = lollypopTable$code_site
lollypopTable$code_site = new_code_site
# lollypopTable_long$code_site = new_code_site[match(lollypopTable_long$code_site, unique(lollypopTable_long$code_site))]


# choose variable for sorting
variableOfInterest = "GPPy_abs_sim" ; variableOfInterest_1 = paste0(variableOfInterest, "_E0") ; variableOfInterest_2 = paste0(variableOfInterest, "_E1A")

# sort
lollypopTable = lollypopTable[sort(drop_units(lollypopTable)[[variableOfInterest_1]], index.return = T)$ix, ]

# transform into a long table
lollypopTable_long = pivot_longer(data = drop_units(lollypopTable), 
                                  grep(colnames(lollypopTable), pattern = paste0(variableList, collapse = "|")), 
                                  names_to = "variable", 
                                  values_to = "value")
# lollypopTable_long = pivot_longer(drop_units(lollypopTable), grep(colnames(lollypopTable), pattern = "GPP|BAI|NPP|Rauto|WVI|Abso"), "variable", "value")
lollypopTable_long$code_site = as.factor(lollypopTable_long$code_site)

# add variable_short column (remove suffixes from variable)
lollypopTable_long$variable_short = vapply(lollypopTable_long$variable, FUN = function(x) strsplit(x, split = paste0("_E0|_E1A|_E2|_CAST"))[[1]][1], FUN.VALUE = "", USE.NAMES = F)
lollypopTable_long$simu = str_match(lollypopTable_long$variable, pattern = "E0|E1A|E2|CAST")[,1]




# PLOT UNE LIGNE PAR STAND (lollypop)
# color1 = hsv(0.65, 0.5, 1, 1) ; color2 = hsv(0.05, 0.8, 0.9)

# pointSize = 3 ; textSize = 11
# sub_lollypopTable_long = subset(lollypopTable_long, variable %in% c(variableOfInterest_1, variableOfInterest_2))

# ggplot(sub_lollypopTable_long) +
#   geom_segment( aes(x = factor(code_site, levels = unique(code_site), ordered = T), 
#                     xend = code_site, y = 0, yend = value),
#                 color = hsv(0, 0, 0, 0.2) ) +
#   geom_point(aes(x = code_site, y = value, color = variable, shape = variable),  size = pointSize) +
#   theme_light() +
#   coord_flip() +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks.y = element_blank(),
#     aspect.ratio = 0.8 * length(unique(sub_lollypopTable_long$code_site)) / 39,
#     # text=element_text(size=textSize), #change font size of all text
#     # axis.text=element_text(size=textSize), #change font size of axis text
#     axis.title=element_text(size=textSize), #change font size of axis titles
#     # plot.title=element_text(size=textSize), #change font size of plot title
#     legend.text=element_text(size=textSize), #change font size of legend text
#     legend.title=element_text(size=textSize) #change font size of legend title
#   ) + 
#   scale_color_manual(labels=c(label_1, label_2),
#                          values=c(color1, color2)) + 
#   scale_shape_manual(labels=c(label_1, label_2),
#                          values=c(16, 18))+ 
#   labs(color=legendTitle, shape = legendTitle) + 
#   ylab(variableOfInterest) + 
#   # ylab("Canopy absorbance") +
#   xlab("Placette")
# 
# folderPlot = paste0("plots/", currentSimulation, "bazardeplot/")
# saveLastGgPlot(folderPlot, plot_height = 720, ratio = 4/3, fileSuffix = ".pdf")



# BOXPLOT ON PHYSIOLOGICAL VARIABLES

# rename simu set
lollypopTable_long$simu[lollypopTable_long$simu == "E0"] = "RN"
lollypopTable_long$simu[lollypopTable_long$simu == "E1A"] = "RS"
lollypopTable_long$simu[lollypopTable_long$simu == "E2"] = "O"

# reorder simu
lollypopTable_long$simu = factor(lollypopTable_long$simu, levels = c('RN', 'RS', 'O'),ordered = TRUE)




# which pair to show p-value?
stat_comparison_pairs = list( c("E0", "E1A"), c("E0", "E2"), c("E1A", "E2") )
stat_comparison_pairs = list( c("RN", "RS"), c("RN", "O"), c("RS", "O") )


# List of plotable variables: 
# "vegAbso_MJm2Y"       "vegAbsorbance"       "vegAbsoPAR_MJm2Y"    "vegAbsorbancePAR"    "GPPy_abs_sim"        "NPPy_abs_sim"        "Rautoy_abs_sim"     
# "BAIy_sim"            "BAIy_mes"            "WVIy_sim"            "WVIy_mes"            "REWmin"              "RU_level_min"        "RU_shortage_max"    
# "incident_yearlyMJm2" "WVIy_m2_mes"         "WVIy_m2_sim"         "GPPy_m2_sim"         "NPPy_m2_sim"         "BAIy_m2_mes"         "BAIy_m2_sim"  
# "transpiration"

# Choice of variables to plot:
variablesPlot = c(
  "transpiration",
  # "REWmin", 
  # "RU_level_min",
  # "GPPy_m2_sim", 
  # "vegAbsorbance", 
  "RU_shortage_max")

# Plot boxplot
ggplot(subset(lollypopTable_long, variable_short %in% variablesPlot &
                grepl(simu, pattern = "RN|RS|O")), 
       aes(x = simu, y = value, fill = simu)) + 
  geom_boxplot() + 
  facet_wrap(. ~ variable_short, scales="free", 
             labeller = labeller(variable_short = c(
               "transpiration" = "Transpiration (mm)",
               # "GPPy_m2_sim" = "GPP (gC/m2/yr)",
               # "vegAbsorbance" = "Absorbance (no unit)",
               "REWmin" = "Min yearly water level (%)",
               "RU_level_min" = "Min yearly water level (mm)",
               "RU_shortage_max" = "Maximum water shortage (mm)"))) +
  ylab("")+
  # coord_cartesian(ylim = c(0,NA))+
  scale_fill_manual(name = "Inventories", values=c("#ddeeff", "#a2cffd", "#1a74d2"))+
  theme(axis.title.x = element_blank(),
        # axis.text.x=element_blank(), axis.ticks.x=element_blank()
  )+
  theme(strip.text.x = element_text(size = 10))+
  # theme(legend.position="none")+
  stat_compare_means(comparisons =  stat_comparison_pairs, 
                     method = "wilcox.test", paired = T, 
                     label = "p.signif")



saveLastGgPlot(folderPlot, plot_height = 720, ratio = 4/3, scale = 0.7, fileSuffix = ".pdf")








# SOME LAST STATISTICS
lollypopTable = removeUnits(lollypopTable)

t.test(lollypopTable$GPPy_abs_sim_E0, lollypopTable$GPPy_abs_sim_E1A, var.equal = T)
t.test(lollypopTable$GPPy_abs_sim_E0, lollypopTable$GPPy_abs_sim_E1A, paired = T)
ks.test(lollypopTable$GPPy_abs_sim_E0, lollypopTable$GPPy_abs_sim_E1A)

(mean(lollypopTable$GPPy_abs_sim_E1A) -  mean(lollypopTable$GPPy_abs_sim_E0)) / mean(lollypopTable$GPPy_abs_sim_E0)
wilcox.test(lollypopTable$GPPy_m2_sim_E0, lollypopTable$GPPy_m2_sim_E1A, paired = T)

(mean(lollypopTable$vegAbsorbance_E1A) -  mean(lollypopTable$vegAbsorbance_E0)) / mean(lollypopTable$vegAbsorbance_E0)
# IF ERROR " impossible de calculer une p-value exacte avec des zéros", it is because values the variables were taken from the exact same simulation (for example in the case of monospecific stands, simulations were the same between regdemo_monosp and regdemo_plurisp)
wilcox.test(lollypopTable$vegAbsorbance_E0, lollypopTable$vegAbsorbance_E1A, paired = T)




