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

currentSimulation = "2024-06-22-ONF/"
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






# # IMPORT DAILY ----

fmOutput = 2 # similar to fmSettings.output
logList = logListList[[fmOutput]]
keepFilter = "2-PB__1-HETpur"

currentSimulation = "2024-04-10_ONF/"
simulationFolderGlobal = paste0("2024_simu_article/", currentSimulation)
folderPlot = paste0("local_plots/", currentSimulation, "bazardeplot/")

dsimuListRU50 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU50/"),
                              logList = logList, keepFilter = keepFilter)
dsimuListRU50 = cleanSimuListNames(dsimuListRU50)


dsimuListRU100 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU100/"),
                               logList = logList, keepFilter = keepFilter)
dsimuListRU100 = cleanSimuListNames(dsimuListRU100)



# Stand scale daily
dsimuListRU50_st = getStandScaleSimuList(dsimuListRU50)
dsimuListRU100_st = getStandScaleSimuList(dsimuListRU100)




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


yvarList = c("GPP", "TR", "REWmin", "LAImaxThisYear", "dbh", "RU_shortage_max")

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
                                               ifelse(standYearTable$composition == "HETsap", 0.70, 
                                                      ifelse(standYearTable$composition == "SAPhet", 0.30, 
                                                             ifelse(standYearTable$composition == "SAPpur", 0,NA))))
  standYearTable$firCompositionRate = 1 - standYearTable$beechCompositionRate
  
  
  standYearTable$nTree = getNTreefromCode_site(standYearTable$code_site)
  
  # meta infos that need to go in simuList
  standYearTable$nha = 0
  standYearTable$standArea_m2 = 0
  standYearTable$nTreeSp = 0
  for(a_code_site in names(simuList_st)){
    standYearTable[standYearTable$code_site == a_code_site, ]$nha = mean(simuList_st[[a_code_site]]$yearlyResults$Nha) # nha of existing tree only (is beech of fir only)
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
  
  standYearTable_meanYear$GPPperLAI = standYearTable_meanYear$GPP / standYearTable_meanYear$LAI #
  standYearTable_meanYear$GPPperGha = standYearTable_meanYear$GPP / standYearTable_meanYear$gha
  
  standYearTable_meanYear$TRperLAI = standYearTable_meanYear$TR / standYearTable_meanYear$LAI
  standYearTable_meanYear$TRperGha = standYearTable_meanYear$TR / standYearTable_meanYear$gha
  
  return(standYearTable_meanYear)
}

standYearTable_meanYear = computeMeanStandYearTable(standYearTable)
standYearTable_ph_meanYear = computeMeanStandYearTable(standYearTable_ph)
standYearTable_sp_meanYear = computeMeanStandYearTable(standYearTable_sp)



# Assemble GPP and TR variables for each species in standYearTable_meanYear
standYearTable_meanYear$GPP_beech = 0
standYearTable_meanYear$GPPperLAI_beech = 0
standYearTable_meanYear$GPPperGha_beech = 0
standYearTable_meanYear$TR_beech = 0
standYearTable_meanYear$TRperLAI_beech = 0
standYearTable_meanYear$TRperGha_beech = 0
standYearTable_meanYear$GPP_fir = 0
standYearTable_meanYear$GPPperLAI_fir = 0
standYearTable_meanYear$GPPperGha_fir = 0
standYearTable_meanYear$TR_fir = 0
standYearTable_meanYear$TRperLAI_fir = 0
standYearTable_meanYear$TRperGha_fir = 0

for(a_code_site in standYearTable_ph_meanYear$code_site){
  sel1 = standYearTable_meanYear$code_site == a_code_site
  sel2 = standYearTable_ph_meanYear$code_site == a_code_site
  
  standYearTable_meanYear[sel1, ]$GPP_beech = standYearTable_ph_meanYear[sel2, ]$GPP
  standYearTable_meanYear[sel1, ]$GPPperLAI_beech = standYearTable_ph_meanYear[sel2, ]$GPPperLAI
  standYearTable_meanYear[sel1, ]$GPPperGha_beech = standYearTable_ph_meanYear[sel2, ]$GPPperGha
  
  standYearTable_meanYear[sel1, ]$TR_beech = standYearTable_ph_meanYear[sel2, ]$TR
  standYearTable_meanYear[sel1, ]$TRperLAI_beech = standYearTable_ph_meanYear[sel2, ]$TRperLAI
  standYearTable_meanYear[sel1, ]$TRperGha_beech = standYearTable_ph_meanYear[sel2, ]$TRperGha
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
}




# Next : compute NBE for GPPperGha and GPPperLAI, TRperGha, TRperLAI for all and each species


# Compute Net Biodiversity Effects (NBE)



NBE_variables = c("GPP", "TR", "RU_shortage_max", "RU_shortage_max_relative", "REWmin", 
                  "GPP_beech", "GPPperLAI_beech", "GPPperGha_beech", "TR_beech", "TRperLAI_beech", "TRperGha_beech", "GPP_fir", "GPPperLAI_fir", "GPPperGha_fir", "TR_fir", "TRperLAI_fir", "TRperGha_fir")
for(var in NBE_variables){
  standYearTable_meanYear[[paste0("NBE_", var)]] = 0
}

for(i in 1:dim(standYearTable_meanYear)[1]){
  aSimu = standYearTable_meanYear[i, ]
  
  if(!grepl(aSimu$composition, pattern = "pur")){
    nTreeSimu = aSimu$nTree
    classSizeSimu = aSimu$classSize
    RUSimu = aSimu$RU
    
    # Find pure fir and pure beech simulations that have the same nTree, classSize, RU and composition
    beechMono = subset(standYearTable_meanYear, nTree == nTreeSimu & classSize == classSizeSimu & RU == RUSimu & composition == "HETpur")
    firMono = subset(standYearTable_meanYear, nTree == nTreeSimu & classSize == classSizeSimu & RU == RUSimu & composition == "SAPpur")
    
    beechProportion = aSimu$beechCompositionRate
    
    numericalColumns = colnames(aSimu)[vapply(aSimu[1,], FUN = is.numeric, FUN.VALUE = FALSE, USE.NAMES = F)]
    mixedAdditiveASimu = beechMono[, numericalColumns] * beechProportion + firMono[, numericalColumns] * (1-beechProportion) 
    
    # prediction vers -> observation
    # NBE = (V_observed - V_predicted) / V_predicted
    for(var in NBE_variables){
      standYearTable_meanYear[i, paste0("NBE_", var)] = (standYearTable_meanYear[i, var] - mixedAdditiveASimu[[var]]) / mixedAdditiveASimu[[var]]
    }
    
  }
}




# Plots Simu yearly ----
folderPlot = paste0("local_plots/", currentSimulation, "allSimuCheck_yearly/")
yvars = c("Gha", "GPP", "NPP", "LAImaxThisYear", "RU_shortage_max", "BiomassOfReservesBeforeRefill")

for(i in 1:length(simuListRU50_st)){
  name_simu = names(simuListRU50_st)[i]
  aplot = simuPlots(simuListRU50_st[[i]], tableName = "y", xvar = "y", yvar = yvars)

  saveGgPlot(aplot, folderPlot,
             plot_height = 960, plot_width = NULL,
             ratio = 4/3,
             scale = 1, fileName = name_simu, fileSuffix = ".pdf")
  rm(aplot)
  rm(name_simu)
}

for(i in 1:length(simuListRU100_st)){
  name_simu = names(simuListRU100_st)[i]
  aplot = simuPlots(simuListRU100_st[[i]], tableName = "y", xvar = "y", yvar = yvars)

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

saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "LAI_nTree", fileSuffix = ".pdf")

# nha (stem/m2) decreases with on classSize (which decreases standArea and maintain nTree) and increases with nTree (density for the same classSize) 
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = nha, x = nTree, size = classSize, color = standArea_m2)) + geom_point() + facet_grid( ~ composition)
# ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = nha, x = nTree, size = standArea_m2)) + geom_point() + facet_grid(classSize ~ composition)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "nha_nTree", fileSuffix = ".pdf")

# class Size increases gha and lower nha, and has no effect on LAI and nTree
# (LAI is linked with gha and nha but ultimately relies on nTree, aka density at equal classSize)
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = gha, x = nha, size = LAI, color = nTree)) + geom_point() + facet_grid(composition ~ classSize)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "gha_nha", fileSuffix = ".pdf")

# gha is approximately the same for each composition, gha depends on density (nTree) and class size
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = gha, x = nTree, size = classSize)) + geom_point() + facet_grid( ~ composition)

saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 4/3,
           scale = 1, fileName = "gha_nTree", fileSuffix = ".pdf")


# Physiological check ----

folderPlot = paste0("local_plots/", currentSimulation, "physiology_basic/")


# LAI / nTree and fir proportion drives GPP and TR.
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = firCompositionRate)) + geom_point() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_LAI", fileSuffix = ".pdf")

# LAI / nTree and fir proportion drives TR
ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = firCompositionRate)) + geom_point() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "TR_LAI", fileSuffix = ".pdf")


# ClassSize has a negligeable effect
# on GPP (fir grows more)
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "GPP_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# on TR
ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "TR_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# on RU_shortage_max
ggplot(standYearTable_meanYear, aes(x = LAI, y = RU_shortage_max, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "RU_shMax_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# on REWmin
ggplot(standYearTable_meanYear, aes(x = LAI, y = REWmin, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "REWmin_LAI_compo_RU_classSize", fileSuffix = ".pdf")

# Same plot, but with code_site shown and saved as file
# myplot = ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = classSize, size = 4 -as.numeric(as.factor(classSize)), label = code_site)) + geom_point() + expand_limits(y = 0) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) +
#   geom_text_repel(size=2, alpha=0.6, color = "black") + guides(size = FALSE)
# saveGgPlot(myplot, plot_height = 1080, plot_width = 1920, scale = 1, plotfolderPlot = folderPlot, fileName = "TR")

# More RU increases maximum water shortage (more possibility to transpire)
ggplot(standYearTable_meanYear, aes(x = as.factor(RU), y = RU_shortage_max, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "RU_shMax_RU", fileSuffix = ".pdf")

# Less RU (hydric stress) increases relative maximum water shortage (more close to the limit)
ggplot(standYearTable_meanYear, aes(x = as.factor(RU), y = RU_shortage_max_relative, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "RU_shMaxRel_RU", fileSuffix = ".pdf")







# PLOTS NBE ----
folderPlot = paste0("local_plots/", currentSimulation, "results/")



# H1.1 et H1.2 : Evaluation of NBE on growth and transpiration
# NBE is globally positive on growth (3.5%) and transpiration (10%). NBE transpiration is more pronounced, especially on sapin-dominant mixture
# RU 100
ggplot(subset(standYearTable_meanYear, isMixed & RU == 100), aes(x = NBE_TR, y = NBE_GPP, color = composition, size = nTree, shape = as.factor(RU))) +
  geom_point(shape = 1) + expand_limits(x = 0) + expand_limits(y = 0)  +
  theme(legend.position = "right") + 
  # scale_shape_manual(values = c(2,1 )) +
  scale_size_continuous(range = c(3,7)) # +facet_wrap(. ~ RU)

saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_NBEtr_RU100", fileSuffix = ".pdf")



# RU 50 and 100
# Hydric stress reduces NBE on growth, not on transpiration
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_TR, y = NBE_GPP, color = composition, size = nTree, shape = as.factor(RU))) +
  geom_point() + expand_limits(x = 0) + expand_limits(y = 0)  +
  theme(legend.position = "right") + 
  scale_shape_manual(values = c(2,1 )) +
  scale_size_continuous(range = c(3,7)) # +facet_wrap(. ~ RU)

saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_NBEtr_RU50RU100", fileSuffix = ".pdf")

# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_TR, y = NBE_GPP, color = composition, size = as.factor(RU))) + geom_point() + expand_limits(x = 0) + expand_limits(y = 0)



# NBE is not exactly identical on transpiration and RU shortage max
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_ru_shortage_max, y = NBE_TR, color = composition)) + geom_point() + expand_limits(x = 0) + expand_limits(y = 0)



# H2.1.a La densité augmente le NBE-growth sur les hetre-dominant, et réduit le NBE-growth sur les sapin-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_GPP, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_nTree", fileSuffix = ".pdf")

# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_GPP, color = composition)) + geom_boxplot() + facet_grid( . ~ classSize) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_nTree_perClassSize", fileSuffix = ".pdf")



# H2.1.b La densité augmente le NBE-transpiration sur tous les types de peuplements
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_TR, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_nTree", fileSuffix = ".pdf")

# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_TR, color = composition)) + geom_boxplot() + facet_grid( . ~ classSize) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_nTree_perClassSize", fileSuffix = ".pdf")


# H2.2.a L'aĝe a un effet positif sur la NBE-gpp sur les hetre-dominant et négatif sur la NBE-gpp sur les sapin-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_GPP, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize", fileSuffix = ".pdf")
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_GPP, color = composition)) + geom_boxplot() + facet_grid( . ~ nTree) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize_perNTree", fileSuffix = ".pdf")



# H2.2.b L'aĝe a un effet nul sur la NBE-water, sauf pour les petites densité (effet positif) et grosse densité (effet négatif chez les hetre-dominant)
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_TR, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_classSize", fileSuffix = ".pdf")
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_TR, color = composition)) + geom_boxplot() + facet_grid( . ~ nTree) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize_perNTree", fileSuffix = ".pdf")



# H3.1 Le stress hydrique réduit la NBE growth, surtout pour les SAP-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_GPP, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_RU", fileSuffix = ".pdf")

# H3.2 le stress joue peu sur la NBE tr
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_TR, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
# --> comment rendre compte de l'effet du stress hydrique sur les dégats du peuplement ?
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_RU", fileSuffix = ".pdf")



# NBE_REW ~ RU
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_REWmin, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)


ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_RU_shortage_max, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_ru_shortage_max_RU", fileSuffix = ".pdf")




# PLOTS NBE SPECIES on Beech ----

# H1.1 et H1.2 : Evaluation of NBE on growth and transpiration

# RU 50 and 100
# Hydric stress reduces NBE on growth, not on transpiration
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_TR_beech, y = NBE_GPP_beech, color = composition, size = nTree, shape = as.factor(RU))) +
  geom_point() + expand_limits(y = 0, x = 0)  +
  theme(legend.position = "right") + 
  scale_shape_manual(values = c(2,1 )) +
  scale_size_continuous(range = c(3,7)) # +facet_wrap(. ~ RU)

saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_NBEtr_RU50RU100_beech", fileSuffix = ".pdf")





# H2.1.a La densité augmente le NBE-growth sur les hetre-dominant, et réduit le NBE-growth sur les sapin-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_GPP_beech, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_nTree_beech", fileSuffix = ".pdf")




# H2.1.b La densité augmente le NBE-transpiration sur tous les types de peuplements
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_TR_beech, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_nTree_beech", fileSuffix = ".pdf")



# H2.2.a L'aĝe a un effet positif sur la NBE-gpp sur les hetre-dominant et négatif sur la NBE-gpp sur les sapin-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_GPP_beech, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize_beech", fileSuffix = ".pdf")



# H2.2.b L'aĝe a un effet nul sur la NBE-water, sauf pour les petites densité (effet positif) et grosse densité (effet négatif chez les hetre-dominant)
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_TR_beech, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_classSize_beech", fileSuffix = ".pdf")
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_TR_beech, color = composition)) + geom_boxplot() + facet_grid( . ~ nTree) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize_perNTree_beech", fileSuffix = ".pdf")



# H3.1 Le stress hydrique réduit la NBE growth, surtout pour les SAP-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_GPP_beech, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_RU_beech", fileSuffix = ".pdf")

# H3.2 le stress joue peu sur la NBE tr
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_TR_beech, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
# --> comment rendre compte de l'effet du stress hydrique sur les dégats du peuplement ?
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_RU_beech", fileSuffix = ".pdf")



# NBE_REW ~ RU
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_REWmin, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)


ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_RU_shortage_max, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_ru_shortage_max_RU", fileSuffix = ".pdf")





# PLOTS NBE SPECIES on Fir ----
# H1.1 et H1.2 : Evaluation of NBE on growth and transpiration

# RU 50 and 100
# Hydric stress reduces NBE on growth, not on transpiration
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_TR_fir, y = NBE_GPP_fir, color = composition, size = nTree, shape = as.factor(RU))) +
  geom_point() + expand_limits(y = 0, x = 0) +
  theme(legend.position = "right") + 
  scale_shape_manual(values = c(2,1 )) +
  scale_size_continuous(range = c(3,7)) # +facet_wrap(. ~ RU)

saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_NBEtr_RU50RU100_fir", fileSuffix = ".pdf")





# H2.1.a La densité augmente le NBE-growth sur les hetre-dominant, et réduit le NBE-growth sur les sapin-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_GPP_fir, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_nTree_fir", fileSuffix = ".pdf")




# H2.1.b La densité augmente le NBE-transpiration sur tous les types de peuplements
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_TR_fir, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_nTree_fir", fileSuffix = ".pdf")



# H2.2.a L'aĝe a un effet positif sur la NBE-gpp sur les hetre-dominant et négatif sur la NBE-gpp sur les sapin-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_GPP_fir, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize_fir", fileSuffix = ".pdf")



# H2.2.b L'aĝe a un effet nul sur la NBE-water, sauf pour les petites densité (effet positif) et grosse densité (effet négatif chez les hetre-dominant)
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_TR_fir, color = composition)) + geom_boxplot() + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_classSize_fir", fileSuffix = ".pdf")
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_TR_fir, color = composition)) + geom_boxplot() + facet_grid( . ~ nTree) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_classSize_perNTree_fir", fileSuffix = ".pdf")



# H3.1 Le stress hydrique réduit la NBE growth, surtout pour les SAP-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_GPP_fir, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEg_RU_fir", fileSuffix = ".pdf")

# H3.2 le stress joue peu sur la NBE tr
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_TR_fir, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
# --> comment rendre compte de l'effet du stress hydrique sur les dégats du peuplement ?
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBEtr_RU_fir", fileSuffix = ".pdf")



# NBE_REW ~ RU
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_REWmin, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)


ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_RU_shortage_max, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + expand_limits(y = 0)
saveGgPlot(folderPlot, 
           plot_height = 480, plot_width = NULL,
           ratio = 3/2,
           scale = 1, fileName = "NBE_ru_shortage_max_RU", fileSuffix = ".pdf")





# 10.04.2024 check daily simulation ----

simuPlots(dsimuListRU50_st$RU50_01_2PB__1HETpur__Narbres_75__Nha_833, tableName = "d", xvar = "d",
          yvar = waterDailyVariables)




