# Camille Rouet 2024 CC-BY-ND

# CONFIGURATION ----
rm(list = ls()) ; gc() # clear

capsisPath = ""
varPath = paste0(capsisPath, "var/") 
workfilesPath = ""
PROGRAM_ON_SERVER = FALSE

source("scripts/define_folders.R")
source("scripts/analysis-of-simulations/CRMethodsForSimulations.R")



# IMPORT ----
# IMPORT SIMULATIONS
fmOutput = 1 # similar to fmSettings.output
logList = logListList[[fmOutput]]
keepFilter = ""

currentSimulation = "2024-04-03_ONF/"
simulationFolderGlobal = paste0("2024_simu_article/", currentSimulation)
folderPlot = paste0("local_plots/", currentSimulation, "bazardeplot/")

simuListRU50 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU50/"),
                                          logList = logList, keepFilter = keepFilter)
simuListRU50 = cleanSimuListNames(simuListRU50)


simuListRU100 = importSimuList(paste0(workfilesPath, "simulations_capsis/", simulationFolderGlobal, "RU100/"),
                              logList = logList, keepFilter = keepFilter)
simuListRU100 = cleanSimuListNames(simuListRU100)

# simuPlots(simuListRU50$RU50_01_2PB__1HETpur__Narbres_75__Nha_833, tableName = "y", xvar = "y", yvar = yearlyVariables)


# CONVERT ----

simuListRU50_st = getStandScaleSimuList(simuListRU50)
simuListRU100_st = getStandScaleSimuList(simuListRU100)


folderPlot = paste0("local_plots/", currentSimulation, "allSimuCheck_yearly/")

yvars = c("Gha", "GPP", "NPP", "LAImaxThisYear", "BiomassOfReserves", "BiomassOfReservesBeforeRefill")

for(i in 1:length(simuListRU50_st)){
  name_simu = names(simuListRU50_st)[i]
  aplot = simuPlots(simuListRU50_st[[i]], tableName = "y", xvar = "y", yvar = yvars)

  saveGgPlot(aplot, saveGgPlot,
             plot_height = 960, plot_width = NULL,
             ratio = 4/3,
             scale = 1, fileName = name_simu, fileSuffix = ".pdf")
  rm(aplot)
  rm(name_simu)
}

for(i in 1:length(simuListRU100_st)){
  name_simu = names(simuListRU100_st)[i]
  aplot = simuPlots(simuListRU100_st[[i]], tableName = "y", xvar = "y", yvar = yvars)

  saveGgPlot(aplot, saveGgPlot,
             plot_height = 960, plot_width = NULL,
             ratio = 4/3,
             scale = 1, fileName = name_simu, fileSuffix = ".pdf")
  rm(aplot)
  rm(name_simu)
}


yvarList = c("GPP", "TR", "REWmin", "LAImaxThisYear", "dbh")
standYearTable_RU50 = makeStandYearTableFromStandScale(simuListRU50_st, yvarList)
standYearTable_RU100 = makeStandYearTableFromStandScale(simuListRU100_st, yvarList)

standYearTable = rbind(standYearTable_RU50, standYearTable_RU100)
standYearTable$RU = getRUfromCode_site(standYearTable$code_site)
standYearTable$classSize = getClassSizefromCode_site(standYearTable$code_site)
standYearTable$composition = getCompositionfromCode_site(standYearTable$code_site)
standYearTable$compositionIndex = getCompositionIndexfromCode_site(standYearTable$code_site)
standYearTable$nTree = getNTreefromCode_site(standYearTable$code_site)
standYearTable$nha = getNhafromCode_site(standYearTable$code_site)
standYearTable$standArea_m2 = standYearTable$nTree / standYearTable$nha *10000
colnames(standYearTable)[colnames(standYearTable) == "LAImaxThisYear"] = "LAI"
standYearTable$gha = standYearTable$nha * pi * (standYearTable$dbh/100 / 2)**2
standYearTable$GPPperLAI = standYearTable$GPP / standYearTable$LAI





# Plot organization of data

# LAI depends solely on number of tree (not on classSize, nor composition) 
ggplot(subset(standYearTable, year == standYearTable$year[1] & RU == 50), aes(y = LAI, x = nTree)) + geom_point() + facet_grid(classSize ~ composition)

# gha is the same for each composition, gha depends on density (nTree) and class size
ggplot(subset(standYearTable, year == standYearTable$year[1] & RU == 50), aes(y = gha, x = nTree, size = classSize)) + geom_point() + facet_grid( ~ composition)



# Plots basic knowledge

# LAI / nTree drives GPP
ggplot(standYearTable, aes(x = LAI, y = GPP, color = as.factor(nTree))) + geom_point() + ylim(c(0, NA))


# Plot advanced results
ggplot(standYearTable, aes(x = composition, y = GPPperLAI, color = as.factor(RU))) + geom_point() + ylim(c(0, NA)) + facet_grid(classSize ~ nTree)

ggplot(standYearTable, aes(x = composition, y = TR, color = as.factor(RU))) + geom_point() + ylim(c(0, NA))  + facet_grid(classSize ~ nTree)
