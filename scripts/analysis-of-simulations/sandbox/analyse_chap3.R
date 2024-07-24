# Camille Rouet 2024 CC-BY-ND

# CONFIGURATION ----
rm(list = ls()) ; gc() # clear

capsisPath = ""
varPath = paste0(capsisPath, "var/") 
workfilesPath = ""
PROGRAM_ON_SERVER = TRUE

source("scripts/define_folders.R")
source("scripts/analysis-of-simulations/CRMethodsForSimulations.R")

library(ggpattern)


# IMPORT YEARLY ----

# IMPORT SIMULATIONS YEARLY
fmOutput = 1 # similar to fmSettings.output
logList = logListList[[fmOutput]]
keepFilter = ""

currentSimulation = "2024-07-12-ONF/"
simulationFolderGlobal = paste0("2024_simu_article/", currentSimulation)
folderPlot = paste0("local_plots/", currentSimulation, "bazardeplot/")

simuListRU50 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU50/"),
                                          logList = logList, keepFilter = keepFilter)
simuListRU50 = cleanSimuListNames(simuListRU50)


simuListRU100 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU100/"),
                              logList = logList, keepFilter = keepFilter)
simuListRU100 = cleanSimuListNames(simuListRU100)


# Assemble RU50 & RU100
simuList = list()
for(a_code_site in names(simuListRU50)){
  simuList[[a_code_site]] = simuListRU50[[a_code_site]]
}
for(a_code_site in names(simuListRU100)){
  simuList[[a_code_site]] = simuListRU100[[a_code_site]]
}

rm(simuListRU50, simuListRU100)






# # # IMPORT DAILY ----
# 
# fmOutput = 2 # similar to fmSettings.output
# logList = logListList[[fmOutput]]
# keepFilter = "2-PB__1-HETpur"
# 
# currentSimulation = "2024-04-10_ONF/"
# simulationFolderGlobal = paste0("2024_simu_article/", currentSimulation)
# folderPlot = paste0("local_plots/", currentSimulation, "bazardeplot/")
# 
# dsimuListRU50 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU50/"),
#                               logList = logList, keepFilter = keepFilter)
# dsimuListRU50 = cleanSimuListNames(dsimuListRU50)
# 
# 
# dsimuListRU100 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU100/"),
#                                logList = logList, keepFilter = keepFilter)
# dsimuListRU100 = cleanSimuListNames(dsimuListRU100)
# 
# 
# 
# # Stand scale daily
# dsimuListRU50_st = getStandScaleSimuList(dsimuListRU50)
# dsimuListRU100_st = getStandScaleSimuList(dsimuListRU100)




# CONVERT YEARLY ----

# Species filtering
simuList_ph = sfiltering(simuList, fmSpeciesIDs = 3)
simuList_sp = sfiltering(simuList, fmSpeciesIDs = 1)

# Remove empty simulations
for(name in names(simuList_ph)){
  if("isEmpty" %in% names(simuList_ph[[name]])){
    simuList_ph[[name]] = NULL
  }
}
for(name in names(simuList_sp)){
  if("isEmpty" %in% names(simuList_sp[[name]])){
    simuList_sp[[name]] = NULL
  }
}


# Stand scale yearly
simuList_st = getStandScaleSimuList(simuList)
simuList_ph_st = getStandScaleSimuList(simuList_ph)
simuList_sp_st = getStandScaleSimuList(simuList_sp)




# Stand-year table


yvarList = c("GPP", "TR", "REWmin", "LAImaxThisYear", "dbh", "RU_shortage_max", "veg_yearlyMJm2", "incident_yearlyMJm2")

# make a stand year table from simuList and add meta info for each simulation in simuList
makeStandYearTableWithMetaInfo = function(simuList_st, simuList, yvarList){
  
  standYearTable = makeStandYearTableFromStandScale_universal(simuList_st, yvarList)
  
  standYearTable$RU = getRUfromCode_site(standYearTable$code_site)
  standYearTable$RU_shortage_max_relative = standYearTable$RU_shortage_max / standYearTable$RU
  standYearTable$classSize = getClassSizefromCode_site(standYearTable$code_site)
  standYearTable$composition = getCompositionfromCodesite_onf(standYearTable$code_site)
  standYearTable$compositionIndex = getCompositionIndexfromCode_site(standYearTable$code_site)
  standYearTable$isMixed = standYearTable$compositionIndex %in% c(2,3)
  
  standYearTable$beechCompositionRate = ifelse(standYearTable$composition == "HETpur", 1, 
                                               ifelse(standYearTable$composition == "HETsap", 0.67, 
                                                      ifelse(standYearTable$composition == "SAPhet", 0.33, 
                                                             ifelse(standYearTable$composition == "SAPpur", 0, NA))))
  standYearTable$firCompositionRate = 1 - standYearTable$beechCompositionRate
  
  
  standYearTable$nTree = getNTreefromCode_site(standYearTable$code_site)
  
  # meta infos that need to go in simuList
  standYearTable$nha = 0
  standYearTable$standArea_m2 = 0
  standYearTable$nTreeSp = 0
  for(a_code_site in names(simuList_st)){
    standYearTable[standYearTable$code_site == a_code_site, ]$nha = mean(simuList_st[[a_code_site]]$yearlyResults$Nha) # nha of existing tree only (if beech of fir only)
    standYearTable[standYearTable$code_site == a_code_site, ]$standArea_m2 = simuList_st[[a_code_site]]$inventory$standArea_m2
    standYearTable[standYearTable$code_site == a_code_site, ]$nTreeSp = length(unique(simuList[[a_code_site]]$yearlyResults$idFmCell))
  }
  
  colnames(standYearTable)[colnames(standYearTable) == "LAImaxThisYear"] = "LAI"
  standYearTable$gha = standYearTable$nha * pi * (standYearTable$dbh/100 / 2)**2
  
  standYearTable$LAI = round(standYearTable$LAI, 2)
  
  # in m2_basalarea/m2_soil
  gha_m2m2 = standYearTable$gha / 10000
  
  # from gC/m2_soil to gC/m2_leaf
  standYearTable$GPPperLAI = standYearTable$GPP / standYearTable$LAI
  
  # from gC/m2_soil to gC/m2_basalarea
  standYearTable$GPPperGha = standYearTable$GPP / gha_m2m2
  
  # from mm/m2_soil to mm/m2_leaf
  standYearTable$TRperLAI = standYearTable$TR / standYearTable$LAI
  
  # from mm/m2_soil to mm/m2_basalarea
  standYearTable$TRperGha = standYearTable$TR / gha_m2m2
  
  
  standYearTable$vegAbsorbed_MJm2 = standYearTable$veg_yearlyMJm2
  standYearTable$incident_MJm2 = standYearTable$incident_yearlyMJm2
  standYearTable$ABS = standYearTable$veg_yearlyMJm2 / standYearTable$incident_yearlyMJm2 # absorbance, no unit (or, per unit)
  
  standYearTable$ABSperGha = standYearTable$ABS / gha_m2m2
  standYearTable$ABSperLAI = standYearTable$ABS / standYearTable$LAI
  
  
  return(standYearTable)
}

standYearTable = makeStandYearTableWithMetaInfo(simuList_st, simuList, yvarList)
standYearTable_ph = makeStandYearTableWithMetaInfo(simuList_ph_st, simuList_ph, yvarList)
standYearTable_sp = makeStandYearTableWithMetaInfo(simuList_sp_st, simuList_sp, yvarList)




# Make a stand-year table with average value for all years
computeMeanStandYearTable = function(standYearTable){
  
  standYearTable_meanYear = aggregate(standYearTable, FUN = mean, by = list(code_site_bis = standYearTable$code_site))
  
  standYearTable_meanYear$code_site = standYearTable_meanYear$code_site_bis
  standYearTable_meanYear$code_site_bis = NULL
  standYearTable_meanYear$RU = getRUfromCode_site(standYearTable_meanYear$code_site)
  standYearTable_meanYear$RU_shortage_max_relative = standYearTable_meanYear$RU_shortage_max / standYearTable_meanYear$RU
  standYearTable_meanYear$classSize = getClassSizefromCode_site(standYearTable_meanYear$code_site)
  standYearTable_meanYear$composition = getCompositionfromCodesite_onf(standYearTable_meanYear$code_site)
  standYearTable_meanYear$compositionIndex = getCompositionIndexfromCode_site(standYearTable_meanYear$code_site)
  standYearTable_meanYear$isMixed = standYearTable_meanYear$compositionIndex %in% c(2,3)
  # standYearTable_meanYear$nTree = getNTreefromCode_site(standYearTable_meanYear$code_site)
  # standYearTable_meanYear$nha = getNhafromCode_site(standYearTable_meanYear$code_site)
  
  
  # no need to recompute the need per year
  # standYearTable_meanYear$GPPperLAI = standYearTable_meanYear$GPP / standYearTable_meanYear$LAI #
  # standYearTable_meanYear$GPPperGha = standYearTable_meanYear$GPP / standYearTable_meanYear$gha
  # 
  # standYearTable_meanYear$TRperLAI = standYearTable_meanYear$TR / standYearTable_meanYear$LAI
  # standYearTable_meanYear$TRperGha = standYearTable_meanYear$TR / standYearTable_meanYear$gha
  
  return(standYearTable_meanYear)
}

standYearTable_meanYear = computeMeanStandYearTable(standYearTable)
standYearTable_ph_meanYear = computeMeanStandYearTable(standYearTable_ph)
standYearTable_sp_meanYear = computeMeanStandYearTable(standYearTable_sp)



# Assemble GPP, TR and ABS variables for each species in standYearTable_meanYear
standYearTable_meanYear$GPP_beech = 0
standYearTable_meanYear$GPPperLAI_beech = 0
standYearTable_meanYear$GPPperGha_beech = 0
standYearTable_meanYear$GPP_fir = 0
standYearTable_meanYear$GPPperLAI_fir = 0
standYearTable_meanYear$GPPperGha_fir = 0
standYearTable_meanYear$TR_fir = 0
standYearTable_meanYear$TRperLAI_fir = 0
standYearTable_meanYear$TRperGha_fir = 0
standYearTable_meanYear$TR_beech = 0
standYearTable_meanYear$TRperLAI_beech = 0
standYearTable_meanYear$TRperGha_beech = 0
standYearTable_meanYear$ABS_beech = 0
standYearTable_meanYear$ABSperLAI_beech = 0
standYearTable_meanYear$ABSperGha_beech = 0
standYearTable_meanYear$ABS_fir = 0
standYearTable_meanYear$ABSperLAI_fir = 0
standYearTable_meanYear$ABSperGha_fir = 0

standYearTable_meanYear$RU_shortage_max_beech = 0
standYearTable_meanYear$RU_shortage_max_fir = 0

# For each site
for(a_code_site in standYearTable_ph_meanYear$code_site){
  sel1 = standYearTable_meanYear$code_site == a_code_site
  sel2 = standYearTable_ph_meanYear$code_site == a_code_site
  
  standYearTable_meanYear[sel1, ]$GPP_beech = standYearTable_ph_meanYear[sel2, ]$GPP
  standYearTable_meanYear[sel1, ]$GPPperLAI_beech = standYearTable_ph_meanYear[sel2, ]$GPPperLAI
  standYearTable_meanYear[sel1, ]$GPPperGha_beech = standYearTable_ph_meanYear[sel2, ]$GPPperGha
  
  standYearTable_meanYear[sel1, ]$TR_beech = standYearTable_ph_meanYear[sel2, ]$TR
  standYearTable_meanYear[sel1, ]$TRperLAI_beech = standYearTable_ph_meanYear[sel2, ]$TRperLAI
  standYearTable_meanYear[sel1, ]$TRperGha_beech = standYearTable_ph_meanYear[sel2, ]$TRperGha
  
  standYearTable_meanYear[sel1, ]$ABS_beech = standYearTable_ph_meanYear[sel2, ]$ABS
  standYearTable_meanYear[sel1, ]$ABSperLAI_beech = standYearTable_ph_meanYear[sel2, ]$ABSperLAI
  standYearTable_meanYear[sel1, ]$ABSperGha_beech = standYearTable_ph_meanYear[sel2, ]$ABSperGha
  
  standYearTable_meanYear[sel1, ]$RU_shortage_max_beech = standYearTable_ph_meanYear[sel2, ]$RU_shortage_max
}

for(a_code_site in standYearTable_sp_meanYear$code_site){
  sel1 = standYearTable_meanYear$code_site == a_code_site
  sel2 = standYearTable_sp_meanYear$code_site == a_code_site
  
  standYearTable_meanYear[sel1, ]$GPP_fir = standYearTable_sp_meanYear[sel2, ]$GPP
  standYearTable_meanYear[sel1, ]$GPPperLAI_fir = standYearTable_sp_meanYear[sel2, ]$GPPperLAI
  standYearTable_meanYear[sel1, ]$GPPperGha_fir = standYearTable_sp_meanYear[sel2, ]$GPPperGha
  
  standYearTable_meanYear[sel1, ]$TR_fir = standYearTable_sp_meanYear[sel2, ]$TR
  standYearTable_meanYear[sel1, ]$TRperLAI_fir = standYearTable_sp_meanYear[sel2, ]$TRperLAI
  standYearTable_meanYear[sel1, ]$TRperGha_fir = standYearTable_sp_meanYear[sel2, ]$TRperGha
  
  standYearTable_meanYear[sel1, ]$ABS_fir = standYearTable_sp_meanYear[sel2, ]$ABS
  standYearTable_meanYear[sel1, ]$ABSperLAI_fir = standYearTable_sp_meanYear[sel2, ]$ABSperLAI
  standYearTable_meanYear[sel1, ]$ABSperGha_fir = standYearTable_sp_meanYear[sel2, ]$ABSperGha
  
  
  standYearTable_meanYear[sel1, ]$RU_shortage_max_fir = standYearTable_sp_meanYear[sel2, ]$RU_shortage_max
}


# correct code_site
standYearTable_meanYear$code_site = vapply(standYearTable_meanYear$code_site, FUN = function(x) paste0(strsplit(x, split = "_")[[1]][-1][-1], collapse = "_"), FUN.VALUE = "", USE.NAMES = F)

# add label for water stressed
standYearTable_meanYear$waterStressed = ifelse(standYearTable_meanYear$RU == 50, "Water stressed", "Normal")




# Compute NBE (Net Biodiversity Effects) ----

NBE_variables = c("RU_shortage_max_relative", "REWmin", 
                  "RU_shortage_max", "RU_shortage_max_beech", "RU_shortage_max_fir",
                  "GPP", "GPPperLAI", "GPPperGha",
                  "GPP_fir", "GPPperLAI_fir", "GPPperGha_fir", 
                  "GPP_beech", "GPPperLAI_beech", "GPPperGha_beech", 
                  "TR", "TRperLAI", "TRperGha",
                  "TR_beech", "TRperLAI_beech", "TRperGha_beech", 
                  "TR_fir", "TRperLAI_fir", "TRperGha_fir",
                  "ABS", "ABSperLAI", "ABSperGha",
                  "ABS_fir", "ABSperLAI_fir", "ABSperGha_fir", 
                  "ABS_beech", "ABSperLAI_beech", "ABSperGha_beech")
for(var in NBE_variables){
  standYearTable_meanYear[[paste0("NBE_", var)]] = 0
}

for(i in 1:dim(standYearTable_meanYear)[1]){
  aSimu = standYearTable_meanYear[i, ]
  
  if(!grepl(aSimu$composition, pattern = "pur")){ # exclude pure stand
    nTreeSimu = aSimu$nTree
    classSizeSimu = aSimu$classSize
    RUSimu = aSimu$RU
    
    # Find pure fir and pure beech simulations that have the same nTree, classSize, RU and composition
    beechMono = subset(standYearTable_meanYear, nTree == nTreeSimu & classSize == classSizeSimu & RU == RUSimu & composition == "HETpur")
    firMono = subset(standYearTable_meanYear, nTree == nTreeSimu & classSize == classSizeSimu & RU == RUSimu & composition == "SAPpur")
    
    beechProportion = aSimu$beechCompositionRate
    
    numericalColumns = colnames(aSimu)[vapply(aSimu[1,], FUN = is.numeric, FUN.VALUE = FALSE, USE.NAMES = F)]
    mixedAdditiveASimu = beechMono[, numericalColumns] * beechProportion + firMono[, numericalColumns] * (1-beechProportion)
    
    perGhaPerLAIColumns = colnames(aSimu)[grepl(colnames(aSimu), pattern = "perGha|perLAI")]
    
    
    # prediction vers -> observation
    # NBE = (V_observed - V_predicted) / V_predicted
    for(var in NBE_variables){
      
      # For variables related to a species basal area of leaf area, NBE is computed based on this species monospecific simulation 
      if(grepl(var, pattern = "perGha|perLAI")){
        if(grepl(var, pattern = "beech")){
          standYearTable_meanYear[i, paste0("NBE_", var)] = (standYearTable_meanYear[i, var] - beechMono[[var]]) / beechMono[[var]]
          next
        }else if(grepl(var, pattern = "fir")){
          standYearTable_meanYear[i, paste0("NBE_", var)] = (standYearTable_meanYear[i, var] - firMono[[var]]) / firMono[[var]]
          next
        }
      }
      # For normal variables : NBE is computed based on predictions on monospecific simulations
      standYearTable_meanYear[i, paste0("NBE_", var)] = (standYearTable_meanYear[i, var] - mixedAdditiveASimu[[var]]) / mixedAdditiveASimu[[var]]
      
    }
    
  }
}


# Filter final table
standYearTable_meanYear_mixed = subset(standYearTable_meanYear, isMixed)
standYearTable_meanYear_mixed_bm = subset(standYearTable_meanYear_mixed, classSize == "3BM")





# Plots Simu yearly ----
folderPlot = paste0("local_plots/", currentSimulation, "allSimuCheck_yearly/")
yvars = c("Gha", "GPP", "NPP", "LAImaxThisYear", "RU_shortage_max", "BiomassOfReservesBeforeRefill")

for(i in 1:length(simuList_st)){
  name_simu = names(simuList_st)[i]
  aplot = simuPlots(simuList_st[[i]], tableName = "y", xvar = "y", yvar = yvars)

  saveGgPlot(aplot, folderPlot,
             plot_height = 960, plot_width = NULL,
             ratio = 4/3,
             scale = 1, fileName = name_simu, fileSuffix = ".pdf")
  rm(aplot)
  rm(name_simu)
}



# Plot organization of data ----
folderPlot = paste0("local_plots/", currentSimulation, "simulation_plan_check/")

# VARIABLES = nTree (LAI) x classSize x RU x composition
# nTree (3) x classSize (3) x RU (2) x composition (4) = 72

# LAI depends solely on nTree (not on classSize, nor composition).
# gha increases with classSize and nTree (density as equal classSize)
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = LAI, x = nTree, size = gha)) + geom_point() + facet_grid(classSize ~ composition)

saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "LAI_nTree", fileSuffix = ".pdf")

# nha (stem/m2) decreases with on classSize (which decreases standArea and maintain nTree) and increases with nTree (density for the same classSize) 
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = nha, x = nTree, size = classSize, color = standArea_m2)) + geom_point() + facet_grid( ~ composition)
# ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = nha, x = nTree, size = standArea_m2)) + geom_point() + facet_grid(classSize ~ composition)
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "nha_nTree", fileSuffix = ".pdf")

# class Size increases gha and lower nha, and has no effect on LAI and nTree
# (LAI is linked with gha and nha but ultimately relies on nTree, aka density at equal classSize)
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = gha, x = nha, size = LAI, color = nTree)) + geom_point() + facet_grid(composition ~ classSize)
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "gha_nha", fileSuffix = ".pdf")

# gha is approximately the same for each composition, gha depends on density (nTree) and class size
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = gha, x = nTree, size = classSize)) + geom_point() + facet_grid( ~ composition)

saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "gha_nTree", fileSuffix = ".pdf")


# Physiological check ----

folderPlot = paste0("local_plots/", currentSimulation, "physiology_basic/")


standYearTable_meanYear_RU_wide = standYearTable_meanYear[ , c("code_site", "RU", "classSize", "nTree","composition", "firCompositionRate", "isMixed", "ABS", "GPP", "TR") ]

standYearTable_meanYear_RU_wide = pivot_wider(data = standYearTable_meanYear_RU_wide, 
                                                    names_from = RU, 
                                                    values_from = c(ABS, GPP, TR))
standYearTable_meanYear_RU_wide$LAI = standYearTable_meanYear$LAI[match(standYearTable_meanYear_RU_wide$code_site, standYearTable_meanYear$code_site)]


# LAI / nTree and fir proportion drives GPP and TR.
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = firCompositionRate)) + geom_point() + 
  expand_limits(y = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_LAI", fileSuffix = ".pdf")

# AND RU
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = firCompositionRate)) + geom_point() + 
  expand_limits(y = 0)+
  geom_segment(data = standYearTable_meanYear_RU_wide, 
               aes(x = LAI, y = GPP_100, xend = LAI, yend = GPP_50, color = firCompositionRate),
               arrow = arrow(length = unit(0.2, "cm")))
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_LAI_RU100to50", fileSuffix = ".pdf")

# AND ABSORBANCE
ggplot(standYearTable_meanYear, aes(x = ABS, y = GPP, color = firCompositionRate)) + geom_point() +
  expand_limits(y = 0) +
  geom_segment(data = standYearTable_meanYear_RU_wide, 
               aes(x = ABS_100, y = GPP_100, xend = ABS_50, yend = GPP_50, color = firCompositionRate),
               arrow = arrow(length = unit(0.2, "cm")))
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_ABS_RU100to50", fileSuffix = ".pdf")



# LAI / nTree and fir proportion drives TR
ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = firCompositionRate)) + geom_point() + expand_limits(y = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "TR_LAI", fileSuffix = ".pdf")


# AND RU
ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = firCompositionRate)) + geom_point() + 
  expand_limits(y = 0)+
  geom_segment(data = standYearTable_meanYear_RU_wide, 
               aes(x = LAI, y = TR_100, xend = LAI, yend = TR_50, color = firCompositionRate),
               arrow = arrow(length = unit(0.2, "cm")))
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "TR_LAI_RU100to50", fileSuffix = ".pdf")

# AND ABSORBANCE
ggplot(standYearTable_meanYear, aes(x = ABS, y = TR, color = firCompositionRate)) + geom_point() +
  expand_limits(y = 0) +
  geom_segment(data = standYearTable_meanYear_RU_wide, 
               aes(x = ABS_100, y = TR_100, xend = ABS_50, yend = TR_50, color = firCompositionRate),
               arrow = arrow(length = unit(0.2, "cm")))
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "TR_ABS_RU100to50", fileSuffix = ".pdf")


# ClassSize has a negligeable effect
# on GPP (fir grows more)
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# on TR
ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "TR_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# on RU_shortage_max
ggplot(standYearTable_meanYear, aes(x = LAI, y = RU_shortage_max, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "RU_shMax_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# on REWmin
ggplot(standYearTable_meanYear, aes(x = LAI, y = REWmin, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "REWmin_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# Same plot, but with code_site shown and saved as file
# myplot = ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = classSize, size = 4 -as.numeric(as.factor(classSize)), label = code_site)) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) +
#   geom_text_repel(size=2, alpha=0.6, color = "black") + guides(size = FALSE)
# saveGgPlot(myplot, plot_height = 1080, plot_width = 1920, scale = 1, plotfolderPlot = folderPlot, fileName = "TR")

# More RU increases maximum water shortage (more possibility to transpire)
ggplot(standYearTable_meanYear, aes(x = as.factor(RU), y = RU_shortage_max, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "RU_shMax_RU", fileSuffix = ".pdf")

# Less RU (hydric stress) increases relative maximum water shortage (more close to the limit)
ggplot(standYearTable_meanYear, aes(x = as.factor(RU), y = RU_shortage_max_relative, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "RU_shMaxRel_RU", fileSuffix = ".pdf")


# Absorbance check
ggplot(subset(standYearTable_meanYear, RU == 100), aes(x = nTree, y = ABS, color = composition, linetype = classSize)) + geom_line() + facet_grid(. ~ composition) + expand_limits(y = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "ABS__nTree", fileSuffix = ".pdf")

ggplot(subset(standYearTable_meanYear, RU == 100), aes(x = nTree, y = GPPperLAI_beech, color = composition, linetype = classSize)) + 
  geom_line() + facet_grid(. ~ composition) + expand_limits(y = 0, x = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPPperLAI_beech__nTree", fileSuffix = ".pdf")


ggplot(subset(standYearTable_meanYear, RU == 100), aes(x = nTree, y = GPPperLAI, color = composition, linetype = classSize)) + 
  geom_line() + facet_grid(. ~ composition) + expand_limits(y = 0, x = 0)
saveGgPlot(folder = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_fir__nTree", fileSuffix = ".pdf")








# PLOTS NBE ----
folderPlot = paste0("local_plots/", currentSimulation, "NBE/")

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = 9, # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
i_pb = 0


variableTypeSuffix = "" ; compoSuffix = ""

for(variableTypeSuffix in c("", "perGha", "perLAI")){
  for(compoSuffix in c("", "_fir", "_beech")){
    setTxtProgressBar(pb, i_pb) # change progress ba
    
    NBE_GPP_var = paste0("NBE_GPP", variableTypeSuffix, compoSuffix)
    NBE_TR_var = paste0("NBE_TR", variableTypeSuffix, compoSuffix)
    NBE_ABS_var = paste0("NBE_ABS", variableTypeSuffix, compoSuffix)
    
    
    NBE_RU_shortage_max_var = paste0("NBE_RU_shortage_max", "", compoSuffix)
    
    standYearTable_meanYear$NBE_GPP_tmp = standYearTable_meanYear[[NBE_GPP_var]]
    standYearTable_meanYear$NBE_TR_tmp = standYearTable_meanYear[[NBE_TR_var]]
    standYearTable_meanYear$NBE_ABS_tmp = standYearTable_meanYear[[NBE_ABS_var]]
    standYearTable_meanYear$NBE_RU_shortage_max_tmp = standYearTable_meanYear[[NBE_RU_shortage_max_var]]
    
    standYearTable_meanYear_mixed = subset(standYearTable_meanYear, isMixed)
    
    
    # widen data set
    standYearTable_meanYear_mixed_RU_wide = standYearTable_meanYear_mixed[ , c("code_site", "RU", "classSize", "nTree","composition", "isMixed", 
                                                                               "NBE_ABS_tmp", "NBE_GPP_tmp", "NBE_TR_tmp", "NBE_RU_shortage_max_tmp") ]
      
    standYearTable_meanYear_mixed_RU_wide = pivot_wider(data = standYearTable_meanYear_mixed_RU_wide, 
                                                        names_from = RU, 
                                                        values_from = c(NBE_ABS_tmp, NBE_GPP_tmp, NBE_TR_tmp, NBE_RU_shortage_max_tmp))
    
    

    
    
    #H0 : Evaluation of NBE on absorbance
    
    
    # NBE_GPP vs NBE_ABS
    # RU 100
    ggplot(subset(standYearTable_meanYear_mixed, RU == 100)) +
      xlab(NBE_ABS_var) + ylab(NBE_GPP_var) +
      geom_point(aes(x = NBE_ABS_tmp, y = NBE_GPP_tmp, color = composition, size = nTree), shape =  10) + expand_limits(x = 0) + expand_limits(y = 0) +
      scale_size_continuous(range = c(2,5)) +
      theme(legend.position = "right")
    
    saveGgPlot(folderPlot = folderPlot, 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_", NBE_ABS_var, "__RU100"), fileSuffix = ".pdf")
    
    # RU 100 to 50
    ggplot(standYearTable_meanYear_mixed) +
      xlab(NBE_ABS_var) + ylab(NBE_GPP_var) +
      geom_point(aes(x = NBE_ABS_tmp, y = NBE_GPP_tmp, color = composition, size = nTree), shape =  10) +
      expand_limits(x = 0) + expand_limits(y = 0) +
      scale_size_continuous(range = c(2,5)) +
      theme(legend.position = "right") + 
      geom_segment(data = standYearTable_meanYear_mixed_RU_wide, 
                   aes(x = NBE_ABS_tmp_100, y = NBE_GPP_tmp_100, xend = NBE_ABS_tmp_50, yend = NBE_GPP_tmp_50),
                   arrow = arrow(length = unit(0.2, "cm")), 
                   color = "black") 
    
    saveGgPlot(folderPlot = folderPlot, 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_", NBE_ABS_var, "__RU100to50"), fileSuffix = ".pdf")
    
    # NBE_TR vs NBE_ABS
    ggplot(standYearTable_meanYear_mixed) +
      xlab(NBE_ABS_var) + ylab(NBE_TR_var) +
      geom_point(aes(x = NBE_ABS_tmp, y = NBE_TR_tmp, color = composition, size = nTree), shape =  10) + expand_limits(x = 0) + expand_limits(y = 0) +
      scale_size_continuous(range = c(2,5)) +
      theme(legend.position = "right") + 
      geom_segment(data = standYearTable_meanYear_mixed_RU_wide, 
                   aes(x = NBE_ABS_tmp_100, y = NBE_TR_tmp_100, xend = NBE_ABS_tmp_50, yend = NBE_TR_tmp_50),
                   arrow = arrow(length = unit(0.2, "cm")), 
                   color = "black") 
    
    saveGgPlot(folderPlot = folderPlot, 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_TR_var, "__f_", NBE_ABS_var, "__RU100to50"), fileSuffix = ".pdf")
    
    
    # absorbance f(density)
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(nTree), y = NBE_ABS_tmp, color = composition)) + geom_boxplot() + 
      expand_limits(y = 0) +
      ylab(NBE_ABS_var)
    
    saveGgPlot(folderPlot = paste0(folderPlot, "res_nTree/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_ABS_var, "__f_nTree"), fileSuffix = ".pdf")
    
    # absorbance f(age)
    ggplot(standYearTable_meanYear_mixed, aes(x = classSize, y = NBE_ABS_tmp, color = composition)) + geom_boxplot() + 
      expand_limits(y = 0) +
      ylab(NBE_ABS_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_ABS_var, "__f_classSize"), fileSuffix = ".pdf")
    
    
    # absorbance f(RU)
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(RU), y = NBE_ABS_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ .) + expand_limits(y = 0) +
      ylab(NBE_ABS_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_RU/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_ABS_var, "__f_RU"), fileSuffix = ".pdf")
    
    
    
    
    # H1.1 et H1.2 : Evaluation of NBE on growth and transpiration
    # NBE is globally positive on growth (3.5%) and transpiration (10%). NBE transpiration is more pronounced, especially on sapin-dominant mixture
    # RU 100
    ggplot(subset(standYearTable_meanYear_mixed, RU == 100), aes(x = NBE_TR_tmp, y = NBE_GPP_tmp, color = composition, size = nTree, shape = as.factor(RU))) +
      xlab(NBE_TR_var) + ylab(NBE_GPP_var) +
      geom_point(shape = 10) + expand_limits(x = 0) + expand_limits(y = 0)  +
      theme(legend.position = "right") + 
      # scale_shape_manual(values = c(2,1 )) +
      scale_size_continuous(range = c(2,5)) # +facet_wrap(. ~ RU)+
    
    saveGgPlot(folderPlot = folderPlot, 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_", NBE_TR_var, "__RU100"), fileSuffix = ".pdf")
    
    
    
    
    
    # RU 50 and 100
    # Hydric stress reduces NBE on growth, not on transpiration
    ggplot(standYearTable_meanYear_mixed) +
      xlab(NBE_TR_var) + ylab(NBE_GPP_var) +
      geom_point(aes(x = NBE_TR_tmp, y = NBE_GPP_tmp, color = composition, size = nTree), shape =  10) + expand_limits(x = 0) + expand_limits(y = 0) +
      scale_size_continuous(range = c(2,5)) +
      theme(legend.position = "right") + 
      geom_segment(data = standYearTable_meanYear_mixed_RU_wide, 
                   aes(x = NBE_TR_tmp_100, y = NBE_GPP_tmp_100, xend = NBE_TR_tmp_50, yend = NBE_GPP_tmp_50),
                   arrow = arrow(length = unit(0.2, "cm")), 
                   color = "black") 
    
    saveGgPlot(folderPlot = folderPlot, 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_", NBE_TR_var, "__RU100to50"), fileSuffix = ".pdf")
    
    # # bonus
    # ggplot(standYearTable_meanYear_mixed, aes(x = NBE_TR_tmp, y = NBE_GPP_tmp, color = composition, size = as.factor(RU))) + geom_point() + 
    #   expand_limits(x = 0) + expand_limits(y = 0) +
    #   xlab(NBE_TR_var) + ylab(NBE_GPP_var)
    # 
    # 
    # # NBE is not exactly identical on transpiration and RU shortage max
    # ggplot(standYearTable_meanYear_mixed, aes(x = NBE_ru_shortage_max, y = NBE_TR_tmp, color = composition)) + geom_point() + 
    #   expand_limits(x = 0) + expand_limits(y = 0) +
    #   ylab(NBE_TR_var)
    
    
    
    # H2.1.a La densité augmente le NBE-growth sur les hetre-dominant, et réduit le NBE-growth sur les sapin-dominant
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(nTree), y = NBE_GPP_tmp, color = composition)) + geom_boxplot() + 
      expand_limits(y = 0) +
      ylab(NBE_GPP_var)
    
    saveGgPlot(folderPlot = paste0(folderPlot, "res_nTree/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_nTree"), fileSuffix = ".pdf")
    
    # bonus
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(nTree), y = NBE_GPP_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ classSize) + expand_limits(y = 0) +
      ylab(NBE_GPP_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_nTree_perClassSize"), fileSuffix = ".pdf")
    
    
    
    # H2.1.b La densité augmente le NBE-transpiration sur tous les types de peuplements
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(nTree), y = NBE_TR_tmp, color = composition)) + geom_boxplot() + 
      expand_limits(y = 0) +
      ylab(NBE_TR_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_nTree/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_TR_var, "__f_nTree"), fileSuffix = ".pdf")
    
    # bonus
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(nTree), y = NBE_TR_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ classSize) + expand_limits(y = 0) +
      ylab(NBE_TR_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_TR_var, "__f_nTree_perClassSize"), fileSuffix = ".pdf")
    
    
    # H2.2.a L'aĝe a un effet positif sur la NBE-gpp sur les hetre-dominant et négatif sur la NBE-gpp sur les sapin-dominant
    ggplot(standYearTable_meanYear_mixed, aes(x = classSize, y = NBE_GPP_tmp, color = composition)) + geom_boxplot() + 
      expand_limits(y = 0) +
      ylab(NBE_GPP_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_classSize"), fileSuffix = ".pdf")
    # bonus
    ggplot(standYearTable_meanYear_mixed, aes(x = classSize, y = NBE_GPP_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ nTree) + expand_limits(y = 0) +
      ylab(NBE_GPP_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_classSize_perNTree"), fileSuffix = ".pdf")
    
    
    
    # H2.2.b L'aĝe a un effet nul sur la NBE-water, sauf pour les petites densité (effet positif) et grosse densité (effet négatif chez les hetre-dominant)
    ggplot(standYearTable_meanYear_mixed, aes(x = classSize, y = NBE_TR_tmp, color = composition)) + geom_boxplot() + 
      expand_limits(y = 0) +
      ylab(NBE_TR_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_TR_var, "__f_classSize"), fileSuffix = ".pdf")
    # bonus
    ggplot(standYearTable_meanYear_mixed, aes(x = classSize, y = NBE_TR_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ nTree) + expand_limits(y = 0) +
      ylab(NBE_TR_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_classSize/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_TR_var, "__f_classSize_perNTree"), fileSuffix = ".pdf")
    
    
    
    # H3.1 Le stress hydrique réduit la NBE growth, surtout pour les SAP-dominant
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(RU), y = NBE_GPP_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ .) + expand_limits(y = 0) +
      ylab(NBE_GPP_var)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_RU/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_GPP_var, "__f_RU"), fileSuffix = ".pdf")
    
    # H3.2 le stress joue peu sur la NBE tr
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(RU), y = NBE_TR_tmp, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ .) + expand_limits(y = 0) +
      ylab(NBE_TR_var)
    # --> comment rendre compte de l'effet du stress hydrique sur les dégats du peuplement ?
    saveGgPlot(folderPlot = paste0(folderPlot, "res_RU/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_TR_var, "__f_RU"), fileSuffix = ".pdf")
    
    # NBE_RU_shortage_max
    ggplot(subset(standYearTable_meanYear_mixed, RU == 100)) +
      xlab(NBE_TR_var) + ylab(NBE_RU_shortage_max_var) +
      geom_point(aes(x = NBE_TR_tmp, y = NBE_RU_shortage_max_tmp, color = composition, size = nTree), shape =  10) + 
      expand_limits(x = 0) + expand_limits(y = 0) +
      scale_size_continuous(range = c(2,5)) +
      theme(legend.position = "right")
    saveGgPlot(folderPlot = folderPlot, 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = paste0(NBE_RU_shortage_max_var, "__f_", NBE_TR_var, "__RU100"), fileSuffix = ".pdf")
    
    # NBE_REW ~ RU
    # ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(RU), y = NBE_REWmin, color = composition)) + geom_boxplot() + 
    #   facet_grid( . ~ .) + expand_limits(y = 0)
    
    
    ggplot(standYearTable_meanYear_mixed, aes(x = as.factor(RU), y = NBE_RU_shortage_max, color = composition)) + geom_boxplot() + 
      facet_grid( . ~ .) + expand_limits(y = 0)
    saveGgPlot(folderPlot = paste0(folderPlot, "res_RU/"), 
               plot_height = 480, plot_width = NULL,
               ratio = 3/2,
               scale = 1, fileName = "NBE_ru_shortage_max_RU", fileSuffix = ".pdf")
    
    i_pb = i_pb +1
  }
}

close(pb)












# Effet propres des variables ----

# a. Plot effet propres sur NBE
variables = c("NBE_GPP", "NBE_TR", "NBE_ABS", "GPP", "TR", "ABS", "REWmin",  
              "NBE_GPP_fir", "NBE_TR_fir", "NBE_ABS_fir", "GPP_fir", "TR_fir", "ABS_fir", 
              "NBE_GPP_beech", "NBE_TR_beech", "NBE_ABS_beech", "GPP_beech", "TR_beech", "ABS_beech",
              "NBE_GPPperLAI", "NBE_TRperLAI", "NBE_ABSperLAI", "GPPperLAI", "TRperLAI", "ABSperLAI", "REWmin",  
              "NBE_GPPperLAI_fir", "NBE_TRperLAI_fir", "NBE_ABSperLAI_fir", "GPPperLAI_fir", "TRperLAI_fir", "ABSperLAI_fir", 
              "NBE_GPPperLAI_beech", "NBE_TRperLAI_beech", "NBE_ABSperLAI_beech", "GPPperLAI_beech", "TRperLAI_beech", "ABSperLAI_beech",
              "NBE_GPPperGha", "NBE_TRperGha", "NBE_ABSperGha", "GPPperGha", "TRperGha", "ABSperGha", "REWmin",  
              "NBE_GPPperGha_fir", "NBE_TRperGha_fir", "NBE_ABSperGha_fir", "GPPperGha_fir", "TRperGha_fir", "ABSperGha_fir", 
              "NBE_GPPperGha_beech", "NBE_TRperGha_beech", "NBE_ABSperGha_beech", "GPPperGha_beech", "TRperGha_beech", "ABSperGha_beech")
for(var in variables){
  mult = 1
  if(grepl(var, pattern = "NBE")){
    mult = 100
  }
  
  folderPlotA = paste0("local_plots/", currentSimulation, "variables_effets_propres/")
  folderPlotB = paste0(folderPlotA, "parComposition/")
  if(grepl(var, pattern = "perLAI")){
    folderPlotA = paste0(folderPlotA, "perLAI/")
    folderPlotB = paste0(folderPlotB, "perLAI/")
  }else if(grepl(var, pattern = "perGha")){
    folderPlotA = paste0(folderPlotA, "perGha/")
    folderPlotB = paste0(folderPlotB, "perGha/")
  }
  
  standYearTable_meanYear_mixed_bm$var_tmp = standYearTable_meanYear_mixed_bm[[var]]
  ggplot(standYearTable_meanYear_mixed_bm, aes(y = mult * var_tmp, x = as.factor(nTree), color = composition)) + geom_point() + 
    ylab(var) +
    expand_limits(x = 0) + expand_limits(y = 0) +
    facet_grid(. ~ waterStressed)
  saveGgPlot(folderPlot = folderPlotA, 
             plot_height = 480, plot_width = NULL,
             ratio = 3/2,
             scale = 1, fileName = paste0(var), fileSuffix = ".pdf")
  
  
  if(grepl(var, pattern = "NBE")){
    next  
  }
  standYearTable_meanYear$var_tmp = standYearTable_meanYear[[var]]
  ggplot(subset(standYearTable_meanYear, classSize == "3BM"), aes(x = composition, y = var_tmp, color = composition, size = nTree)) +
    ylab(var) +
    geom_jitter(width = 0.15) + 
    expand_limits(y = 0) + facet_grid(. ~ waterStressed)
  saveGgPlot(folderPlot = folderPlotB, 
             plot_height = 480, plot_width = NULL,
             ratio = 3/2,
             scale = 1, fileName = paste0(var, "_perComposition"), fileSuffix = ".pdf")
}




# b. Modèle linéaire
# cela permet d'évaluer tous les effets et leurs interactions
lm1 = lm(data = standYearTable_meanYear_mixed_bm, 1000 * NBE_GPP ~ composition + waterStressed * as.factor(nTree))
summary(lm1)
mean(subset(standYearTable_meanYear_mixed_bm, composition == "SAPhet")$NBE_TR) / mean(subset(standYearTable_meanYear_mixed_bm, composition == "HETsap")$NBE_TR)
                                                         


# c. Calcul de variation de NBE à partir d'une simulation de référence : RU100, Densité moyenne
# (une simu de référence pour hetre dominant et pour sapin dominant)
refSimuBeechDom = subset(standYearTable_meanYear_mixed_bm, RU == 100 & nTree == 54 & composition == "HETsap")
refSimuFirDom = subset(standYearTable_meanYear_mixed_bm, RU == 100 & nTree == 54 & composition == "SAPhet")

standYearTable_meanYear_relativeto_refSimuBeechDom = standYearTable_meanYear_mixed_bm
standYearTable_meanYear_relativeto_refSimuFirDom = standYearTable_meanYear_mixed_bm
numericalCol = vapply(standYearTable_meanYear_mixed_bm[1, ], FUN = function(x) is.numeric(x), FUN.VALUE = T, USE.NAMES = F)
numericalCol = ! colnames(standYearTable_meanYear_mixed_bm) %in% c("code_site", "classSize", "composition", "isMixed", "year", "LAI", "dbh", "RU",
                                                                               "compositionIndex", "beechCompositionRate", "firCompositionRate",
                                                                               "nTree", "nha", "standArea_m2", "nTreeSp", "gha", "waterStressed")
# Pour chaque simulation et variable, on divise le résultat par celui de la valeur de la simulation de référence
for(i in 1:dim(standYearTable_meanYear_mixed_bm)[1]){
  standYearTable_meanYear_relativeto_refSimuBeechDom[i, numericalCol] = round((standYearTable_meanYear_mixed_bm[i, numericalCol] / refSimuBeechDom[, numericalCol] - 1) * 100, 1)
  standYearTable_meanYear_relativeto_refSimuFirDom[i, numericalCol] = round((standYearTable_meanYear_mixed_bm[i, numericalCol] / refSimuFirDom[, numericalCol] - 1) * 100, 1)
}

# selection de variables à afficher
variables = c("nTree", "RU", "NBE_GPP", "NBE_TR", "NBE_ABS", "GPP", "TR", "ABS", "REWmin")
variables_beechfir = c("nTree", "NBE_GPP_fir", "NBE_TR_fir", "NBE_ABS_fir", "GPP_fir", "TR_fir", "ABS_fir", "REWmin", "NBE_GPP_beech", "NBE_TR_beech", "NBE_ABS_beech", "GPP_beech", "TR_beech", "ABS_beech")


# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas HETsap) (en %)
selection1 = standYearTable_meanYear_relativeto_refSimuBeechDom$composition == "HETsap" &
  (standYearTable_meanYear_relativeto_refSimuBeechDom$RU == 100 | standYearTable_meanYear_relativeto_refSimuBeechDom$nTree == 54)
standYearTable_meanYear_relativeto_refSimuBeechDom[selection1,
                                                            variables]


# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas SAPhet) (en %)
selection2 = standYearTable_meanYear_relativeto_refSimuFirDom$composition == "SAPhet" &
  (standYearTable_meanYear_relativeto_refSimuFirDom$RU == 100 | standYearTable_meanYear_relativeto_refSimuFirDom$nTree == 54)
standYearTable_meanYear_relativeto_refSimuFirDom[selection2,
                                                          variables]



# Check effect on species
# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas HETsap) (en %)
selection1 = standYearTable_meanYear_relativeto_refSimuBeechDom$composition == "HETsap" &
  (standYearTable_meanYear_relativeto_refSimuBeechDom$RU == 100 | standYearTable_meanYear_relativeto_refSimuBeechDom$nTree == 54)
standYearTable_meanYear_relativeto_refSimuBeechDom[selection1,
                                                            variables_beechfir]


# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas SAPhet) (en %)
selection2 = standYearTable_meanYear_relativeto_refSimuFirDom$composition == "SAPhet" &
  (standYearTable_meanYear_relativeto_refSimuFirDom$RU == 100 | standYearTable_meanYear_relativeto_refSimuFirDom$nTree == 54)
standYearTable_meanYear_relativeto_refSimuFirDom[selection2,
                                                          variables_beechfir]


# Résumé : on a les effets isolés (par rapport à une simulation de référence) de la densité et l'effet de la RU 
# sur les valeurs de GPP, TR et ABS ainsi que sur les effets mélanges sur ces variables
# On a aussi les valeurs pour chaque espèce
# On a un  effet pour les mélanges à hêtre dominant, et un autre pour les mélanges à sapins dominants






# Check NBE ~ RU par espece ----
summary(standYearTable_meanYear$GPP - (standYearTable_meanYear$GPP_beech + standYearTable_meanYear$GPP_fir) )



# SimuRef is
# aComposition is HETsap or SAPhet
getValues = function(fourSimuRef, aComposition){
  
  if(aComposition == "HETsap"){
    compositionRates = c(0.67, 0.33)
  }else{
    compositionRates = c(0.33, 0.67)
  }
  
  # Global
  GPP = subset(fourSimuRef, composition == aComposition)$GPP
  GPP_pred =  compositionRates[1] * subset(fourSimuRef, composition == "HETpur")$GPP + compositionRates[2] *  subset(fourSimuRef, composition == "SAPpur")$GPP
  NBE_GPP = ( GPP - GPP_pred ) / GPP_pred
  
  NBE_GPP
  subset(fourSimuRef, composition == aComposition)$NBE_GPP
  
  # Beech
  GPP_beech = subset(fourSimuRef, composition == aComposition)$GPP_beech
  GPP_pred_beech =  compositionRates[1] * subset(fourSimuRef, composition == "HETpur")$GPP_beech + compositionRates[2] *  subset(fourSimuRef, composition == "SAPpur")$GPP_beech
  NBE_GPP_beech = ( GPP_beech - GPP_pred_beech ) / GPP_pred_beech
  
  NBE_GPP_beech
  subset(fourSimuRef, composition == aComposition)$NBE_GPP_beech
  
  
  # Fir
  GPP_fir = subset(fourSimuRef, composition == aComposition)$GPP_fir
  GPP_pred_fir =  compositionRates[1] * subset(fourSimuRef, composition == "HETpur")$GPP_fir + compositionRates[2] *  subset(fourSimuRef, composition == "SAPpur")$GPP_fir
  NBE_GPP_fir = ( GPP_fir - GPP_pred_fir ) / GPP_pred_fir
  
  NBE_GPP_fir
  subset(fourSimuRef, composition == aComposition)$NBE_GPP_fir
  
  
  tab = tibble(composition= aComposition,
               RU = unique(fourSimuRef$RU),
               nTree = unique(fourSimuRef$nTree),
               GPP = c(GPP), 
               GPP_pred = c(GPP_pred),
               NBE_GPP = c(NBE_GPP), 
               GPP_beech = c(GPP_beech), 
               GPP_pred_beech = c(GPP_pred_beech),
               NBE_GPP_beech = c(NBE_GPP_beech), 
               GPP_fir = c(GPP_fir), 
               GPP_pred_fir = c(GPP_pred_fir),
               NBE_GPP_fir = c(NBE_GPP_fir))
}



# Donner une RU et une densité pour comparer les résultats
tab_RU100 = getValues(subset(standYearTable_meanYear, RU == 100 & nTree == 30 & classSize == "3BM"), "HETsap")
tab_RU50 = getValues(subset(standYearTable_meanYear, RU == 100 & nTree == 30 & classSize == "3BM"), "SAPhet")

rbind(tab_RU100, tab_RU50)





# Bar plot ----

plot_variables = c("NBE_GPP", "NBE_TR", "NBE_ABS")
pivot_variables = unique(c(plot_variables, 
                           "GPP", "TR", "ABS", "NBE_GPP", "NBE_TR", "NBE_ABS",
                          "GPP_beech", "TR_beech", "ABS_beech", "NBE_GPP_beech", "NBE_TR_beech", "NBE_ABS_beech", 
                          "GPP_fir", "TR_fir", "ABS_fir", "NBE_GPP_fir", "NBE_TR_fir", "NBE_ABS_fir"))
standYearTable_meanYear_mixed_bm_longer = pivot_longer(standYearTable_meanYear_mixed_bm[, c("nTree", "RU", "composition", pivot_variables)], cols = pivot_variables, names_to = "variable", values_to = "value")

# Add a species column
standYearTable_meanYear_mixed_bm_longer$species = sub(".*_(beech|fir).*", "\\1", standYearTable_meanYear_mixed_bm_longer$variable)

# Set species at "all" if not specifically beech or fir
standYearTable_meanYear_mixed_bm_longer$species[!standYearTable_meanYear_mixed_bm_longer$species %in% c("beech", "fir")]  = "all"

# Rename the variable column to get reed of species
standYearTable_meanYear_mixed_bm_longer = mutate(.data = standYearTable_meanYear_mixed_bm_longer, 
                                                 variable = sub("_(beech|fir)", "", variable))


# Ajouter une colonne pour la taille des triangles basée sur l'indice
standYearTable_meanYear_mixed_bm_longer$triangle_size <- as.numeric(factor(standYearTable_meanYear_mixed_bm_longer$nTree, levels = c(30, 54, 75)))







# Create a beautiful barplot
# VerticalXX input should be use only if plot_variables are POSITIVES
barPlot = function(title, standYearTable_meanYear_mixed_bm_longer, plot_variables, alphaVar, alphaVarName, alphaVarLabels, 
                   verticalVariable = NULL, verticalVariableUpDownLabels = NULL, y_ticks = NULL, 
                   hashInsteadOfAlpha = F, addPoint = F, hashType = "circle", inverseHashDensity = F,
                   addOutline = F){
  
  standYearTable_meanYear_mixed_bm_longer$alphaVariable = standYearTable_meanYear_mixed_bm_longer[[alphaVar]]
  y_ticks_labels = y_ticks
  
  # Assign a negative value if data of verticalVariable should be display below
  if(!is.null(verticalVariable)){
    standYearTable_meanYear_mixed_bm_longer$value = standYearTable_meanYear_mixed_bm_longer$value * ifelse(standYearTable_meanYear_mixed_bm_longer[[verticalVariable]] == names(verticalVariableUpDownLabels)[1], 1, -1)
    y_ticks_labels = paste0("+ ", abs(y_ticks))
  }
  
  subtable = subset(standYearTable_meanYear_mixed_bm_longer, variable %in% plot_variables)
  
  # (optional) Use pattern instead of alpha to highlight a variable
  if(hashInsteadOfAlpha){
    theAes = aes(y = 100 * value, x = as.factor(nTree), fill = variable, pattern_density = as.factor(alphaVariable))
    theGeomCol = geom_col_pattern(position = "dodge", width = 0.7, alpha = 0.4,
                                  pattern = hashType, pattern_size = 0, pattern_fill = "#414141", pattern_angle = -60, 
                                  pattern_spacing = 0.03, pattern_alpha = 0.4, color = if(addOutline){"#414141"}else{rgb(0,0,0,0)}, size = if(addOutline){0.2}else{0})
    densities = c(0.5, 0)
    pattern_opt = scale_pattern_density_discrete(range = ifelse(rep(inverseHashDensity, 2), rev(densities), densities), 
                                                 name = alphaVarName, labels = alphaVarLabels)
    alpha_opt = NULL
  }else{
    # Normal
    theAes = aes(y = 100 * value, x = as.factor(nTree), fill = variable, alpha = as.factor(alphaVariable), shape = as.factor(alphaVariable))
    theGeomCol = geom_col(position = "dodge", width = 0.7, color = if(addOutline){"#414141"}else{rgb(0,0,0,0)}, size = if(addOutline){0.2}else{0})
    pattern_opt = NULL
    alpha_opt = scale_alpha_manual(values = c(0.4, 1), 
                       name = alphaVarName, labels = alphaVarLabels)
  }
  
  # (optional) Add points at the edge of the bars
  if(addPoint){
    colorPoint = rgb(0.4,0.4,0.4)
    shapesPoint = c(1, 13)
    theGeomPoint = geom_point(show.legend = F, theAes, position=position_dodge(width= 0.7), color = colorPoint)
    shape_opt = scale_shape_manual(values = shapesPoint)
  }else{
    theGeomPoint = NULL
    shape_opt = NULL
  }
  
  # Base ggplot
  theplot = ggplot(subtable, theAes) + 
    theGeomCol + pattern_opt +
    theGeomPoint + shape_opt +
    geom_point(aes(y = 0, size = triangle_size), shape = 17, color = rgb(0.2,0.2,0.2)) +  # Ajouter des triangles avec tailles différentes
    geom_abline(slope = 0, intercept = 0, alpha = 0.25, size = 0.25) +
    scale_size_continuous(range = c(2,6), breaks = c(1, 2, 3), 
                          name = "Densité", labels = c("Faible", "Intermédiaire", "Forte")) +
    scale_fill_manual(values = c('NBE_GPP' = 'grey', 'NBE_TR' = 'deepskyblue', 'NBE_ABS' = "orange")) +
    alpha_opt +
    labs(x = title,
         y = "") +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.x = element_blank(), 
          # strip.text = element_blank() # retirer les noms de facettes
    ) +
    facet_grid(.~variable, labeller = labeller(variable = c("NBE_ABS" = "NBE(Absorbance)",  "NBE_GPP" = "NBE(GPP)", "NBE_TR" = "NBE(Transpiration)"))) +
    guides( fill = guide_none(),
            alpha = guide_legend(order = 1),
            pattern_density = guide_legend(order = 1, override.aes = list(fill = "orange", color = "#414141", shape = NA)))
  
  # Add customized y-axis ticks
  if(!is.null(y_ticks)){
    theplot = theplot + scale_y_continuous(breaks = y_ticks, labels = y_ticks_labels) 
  }

  
  # Add annotation to axes
  if(!is.null(verticalVariable)){
    
    # Ajouter l'annotation uniquement pour le premier panneau 
    theplot = theplot + 
      geom_text(data = data.frame(
        nTree= -Inf,  # Positionner à l'infini pour se situer au bord du graphique
        value = Inf,  # Positionner à l'infini pour se situer au bord du graphique
        variable = "NBE_ABS",  # Indiquer le panneau spécifique
        alphaVariable = subtable$alphaVariable[1],
        label = verticalVariableUpDownLabels[1]
      ), aes(x = nTree, y = value, label = label), hjust = 0, vjust = 1.1, show.legend = F) + 
      geom_text(data = data.frame(
        nTree= -Inf,  # Positionner à l'infini pour se situer au bord du graphique
        value = -Inf,  # Positionner à l'infini pour se situer au bord du graphique
        variable = "NBE_ABS",  # Indiquer le panneau spécifique
        alphaVariable = subtable$alphaVariable[1],
        label = verticalVariableUpDownLabels[2]
      ), aes(x = nTree, y = value, label = label), hjust = 0, vjust = -1.1, show.legend = F)
  }
  
  if(addPoint){
    theplot = theplot + guides( alpha = guide_legend(order = 1, override.aes = list(fill = "orange", shape = shapesPoint, color = colorPoint, linetype = 0)))
  }else{
    theplot = theplot + guides( alpha = guide_legend(order = 1, override.aes = list(fill = "orange", shape = NA)))
  }
  
  return(theplot)
}


folderPlot = paste0("local_plots/", currentSimulation, "variables_effets_propres/NBE_barplot/")



# All NBE on HETsap, SAPhet, per density
barPlot("Effet du mélange net (%)", 
        subset(standYearTable_meanYear_mixed_bm_longer, RU == 100 & species == "all"), plot_variables, 
        alphaVar = "composition", alphaVarName = "Mélange\ndominé par le", alphaVarLabels =c("hêtre","sapin"))

saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_barplot_base", fileSuffix = ".pdf")

# RU effect on NBE(GPP) 
barPlot("Effet du mélange net (%)", 
        subset(standYearTable_meanYear_mixed_bm_longer, species == "all"), plot_variables, 
        alphaVar = "RU", alphaVarName = "RU", alphaVarLabels = c("50" = "50", "100" = "100"),
        verticalVariable = "composition", verticalVariableUpDownLabels = c("HETsap" = "Hêtre dom", "SAPhet" = "Sapin dom"),
        hashInsteadOfAlpha = T, hashType = "circle", y_ticks = c(-15, -10, -5, 0, 5, 10, 15, 20), addOutline = T) 
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_barplot_RU_effect", fileSuffix = ".pdf")

# RU effect on NBE(GPP) (HETsap)
barPlot("Effet du mélange net (%)\n(hêtre majoritaire)", 
        subset(standYearTable_meanYear_mixed_bm_longer, species == "all" & composition == "HETsap"), plot_variables, 
        alphaVar = "RU", alphaVarName = "RU", alphaVarLabels = c("50" = "50", "100" = "100"),
        #verticalVariable = "composition", verticalVariableUpDownLabels = c("HETsap" = "Hêtre dom", "SAPhet" = "Sapin dom"),
        hashInsteadOfAlpha = T, hashType = "circle", addOutline = T) 
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_barplot_RU_effect_HETsap", fileSuffix = ".pdf")

# RU effect on NBE(GPP) (SAPhet)
barPlot("Effet du mélange net (%)\n(sapin majoritaire)", 
        subset(standYearTable_meanYear_mixed_bm_longer, species == "all" & composition == "SAPhet"), plot_variables, 
        alphaVar = "RU", alphaVarName = "RU", alphaVarLabels = c("50" = "50", "100" = "100"),
        #verticalVariable = "composition", verticalVariableUpDownLabels = c("HETsap" = "Hêtre dom", "SAPhet" = "Sapin dom"),
        hashInsteadOfAlpha = T, hashType = "circle", addOutline = T) 
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_barplot_RU_effect_SAPhet", fileSuffix = ".pdf")


# NBE(GPP) per species
barPlot("Effet du mélange net (%)\n(hêtre majoritaire)", 
        subset(standYearTable_meanYear_mixed_bm_longer, RU == 100 & composition == "HETsap" & species %in% c("beech", "fir")), plot_variables, 
        alphaVar = "species", alphaVarName = "Espèce", alphaVarLabels =c("beech" = "Hêtre", "fir" = "Sapin"),
        addPoint = F, hashInsteadOfAlpha = T, hashType = "stripe", inverseHashDensity = T)
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_barplot_bySpecies_HETsap", fileSuffix = ".pdf")


barPlot("Effet du mélange net (%)\n(sapin majoritaire)",
        subset(standYearTable_meanYear_mixed_bm_longer, RU == 100 & composition == "SAPhet" & species %in% c("beech", "fir")), plot_variables, 
        alphaVar = "species", alphaVarName = "Espèce", alphaVarLabels =c("beech" = "Hêtre", "fir" = "Sapin"),
        addPoint = F, hashInsteadOfAlpha = T, hashType = "stripe", inverseHashDensity = T)
saveGgPlot(folderPlot = folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_barplot_bySpecies_SAPhet", fileSuffix = ".pdf")




# old version
# ggplot(subset(standYearTable_meanYear_mixed_bm_longer, variable %in% plot_variables & RU ==100), 
#        aes(x = 100 * value * ifelse(composition == "HETsap", -1 ,1), y = as.factor(nTree), fill = variable, alpha = composition)) +
#   geom_bar(stat = 'identity', position = 'identity', width = 0.7) +
#   geom_point(aes(x = triangle_xpos, size = triangle_size), shape = 5, color = 'black') +  # Ajouter des triangles avec tailles différentes
#   scale_size_continuous(range = c(2,5), breaks = c(1, 2, 3), 
#                         name = "Densité", labels = c("Faible", "Intermédiaire", "Forte")) +
#   scale_fill_manual(values = c('NBE_GPP' = 'grey', 'NBE_TR' = 'deepskyblue', 'NBE_ABS' = "orange"),
#                     name = "NBE", labels = c("NBE_ABS" = "Absorbance",  "NBE_GPP" = "GPP", "NBE_TR" = "Transpiration")) +
#   scale_alpha_manual(values = c('HETsap' = 0.5, 'SAPhet' = 1), 
#                     name = "Espèce\ndominante", labels = c("Hêtre","Sapin")) +
#   labs(x = "Effet du mélange net (%)",
#        y = "") +
#   theme_minimal() +
#   theme(legend.position = "top",
#         axis.text.y = element_blank(), 
#         # strip.text = element_blank() # retirer les noms de facettes
#         ) +
#   facet_grid(variable ~ ., labeller = labeller(variable = c("NBE_ABS" = "NBE(Absorbance)",  "NBE_GPP" = "NBE(GPP)", "NBE_TR" = "NBE(Transpiration)"))) +
#   guides( alpha = guide_legend(override.aes = list(fill = "orange", shape = NA)),  # Retirer les points de la légende
#           fill = guide_none() ) +
#   scale_x_continuous(breaks = c(-10, -5, 0, 5,10), labels = c(10, 5, 0, 5, 10))

# Next : 
#   Faire un graphe avec les NBE beech et les NBE sapin ? (utile pour montrer la répartition du NBE, mais compliqué en raison de la catégorie "espèce dominante hetre/sapin" qui va rendre confus) 
#   Faire un graphe avec les NBE RU50  et les NBE RU 100 ? (pas besoin car la RU n'affecte que la NBE GPP)


