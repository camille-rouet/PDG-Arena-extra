# Camille Rouet 2024 CC-BY-ND

# CONFIGURATION ----
rm(list = ls()) ; gc() # clear

capsisPath = ""
varPath = paste0(capsisPath, "var/") 
workfilesPath = ""
PROGRAM_ON_SERVER = FALSE

source("scripts/define_folders.R")
source("scripts/analysis-of-simulations/CRMethodsForSimulations.R")



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




# Stand-year table, average for all years

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



# Assemble GPP and TR variables for each species in standYearTable_meanYear
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
      # xlab(NBE_ABS_var) + ylab(NBE_GPP_var) +
      geom_point(aes(x = NBE_RU_shortage_max_tmp, y = NBE_GPP_tmp, color = composition, size = nTree), shape =  10) + expand_limits(x = 0) + expand_limits(y = 0) +
      scale_size_continuous(range = c(2,5)) +
      theme(legend.position = "right") + 
      geom_segment(data = standYearTable_meanYear_mixed_RU_wide, 
                   aes(x = NBE_RU_shortage_max_tmp_100, y = NBE_GPP_tmp_100, xend = NBE_RU_shortage_max_tmp_50, yend = NBE_GPP_tmp_50),
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











# 10.04.2024 check daily simulation ----

simuPlots(dsimuListRU50_st$RU50_01_2PB__1HETpur__Narbres_75__Nha_833, tableName = "d", xvar = "d",
          yvar = waterDailyVariables)







# test ----

standYearTable_meanYear_mixed_bm = subset(standYearTable_meanYear_mixed, classSize == "3BM")

standYearTable_meanYear_mixed_bm$waterStressed = standYearTable_meanYear_mixed_bm$RU == 50
ggplot(standYearTable_meanYear_mixed_bm, aes(y = 1000 * NBE_GPP, x = as.factor(nTree), color = composition)) + geom_boxplot() + 
  expand_limits(x = 0) + expand_limits(y = 0) +
  facet_grid(. ~ waterStressed)

# linear model
# cela permet d'évaluer tous les effets et leurs interactions
lm1 = lm(data = standYearTable_meanYear_mixed_bm, 1000 * NBE_GPP ~ composition + waterStressed * as.factor(nTree))
summary(lm1)
mean(subset(standYearTable_meanYear_mixed_bm, composition == "SAPhet")$NBE_TR) / mean(subset(standYearTable_meanYear_mixed_bm, composition == "HETsap")$NBE_TR)
                                                         


# Calcul de variation de NBE à partir d'une simulation de référence : RU100, Densité moyenne
# (une simu de référence pour hetre dominant et pour sapin dominant)
refSimuBeechDom = subset(standYearTable_meanYear_mixed_bm, RU == 100 & nTree == 54 & composition == "HETsap")
refSimuFirDom = subset(standYearTable_meanYear_mixed_bm, RU == 100 & nTree == 54 & composition == "SAPhet")

standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom = standYearTable_meanYear_mixed_bm
standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom = standYearTable_meanYear_mixed_bm
numericalCol = vapply(standYearTable_meanYear_mixed_bm[1, ], FUN = function(x) is.numeric(x), FUN.VALUE = T, USE.NAMES = F)
numericalCol = ! colnames(standYearTable_meanYear_mixed_bm) %in% c("code_site", "classSize", "composition", "isMixed", "year", "LAI", "dbh", "RU",
                                                                               "compositionIndex", "beechCompositionRate", "firCompositionRate",
                                                                               "nTree", "nha", "standArea_m2", "nTreeSp", "gha")

for(i in 1:dim(standYearTable_meanYear_mixed_bm)[1]){
  standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom[i, numericalCol] = round((standYearTable_meanYear_mixed_bm[i, numericalCol] / refSimuBeechDom[, numericalCol] - 1) * 100, 1)
  standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom[i, numericalCol] = round((standYearTable_meanYear_mixed_bm[i, numericalCol] / refSimuFirDom[, numericalCol] - 1) * 100, 1)
}

variables = c("nTree", "RU", "NBE_GPP", "NBE_TR", "NBE_ABS", "GPP", "TR", "ABS", "REWmin")
variables_beechfir = c("nTree", "NBE_GPP_fir", "NBE_TR_fir", "NBE_ABS_fir", "GPP_fir", "TR_fir", "ABS_fir", "REWmin", "NBE_GPP_beech", "NBE_TR_beech", "NBE_ABS_beech", "GPP_beech", "TR_beech", "ABS_beech")
variables_fir = c("nTree", "RU", "NBE_GPP", "NBE_TR", "NBE_ABS", "GPP", "TR", "ABS", "REWmin", )


# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas HETsap) (en %)
selection1 = standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom$composition == "HETsap" &
  (standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom$RU == 100 | standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom$nTree == 54)
standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom[selection1,
                                                            variables]


# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas SAPhet) (en %)
selection2 = standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom$composition == "SAPhet" &
  (standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom$RU == 100 | standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom$nTree == 54)
standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom[selection2,
                                                          variables]



# Check effect on species
# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas HETsap) (en %)
selection1 = standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom$composition == "HETsap" &
  (standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom$RU == 100 | standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom$nTree == 54)
standYearTable_meanYear_mixed_bm_relativeto_refSimuBeechDom[selection1,
                                                            variables_beechfir]


# Effet unique de la densité forte, de la densité faible, de la RU50 sur la NBE (cas SAPhet) (en %)
selection2 = standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom$composition == "SAPhet" &
  (standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom$RU == 100 | standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom$nTree == 54)
standYearTable_meanYear_mixed_bm_relativeto_refSimuFirDom[selection2,
                                                          variables_beechfir]
