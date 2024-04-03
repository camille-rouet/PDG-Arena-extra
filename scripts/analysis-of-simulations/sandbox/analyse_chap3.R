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

currentSimulation = "2024-03-29_ONF/"
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


# folderPlot = paste0("local_plots/", currentSimulation, "allSimuCheck_yearly/")
# 
# for(i in 1:length(simuListRU50_st)){
#   name_simu = names(simuListRU50_st)[i]
#   aplot = simuPlots(simuListRU50_st[[i]], tableName = "y", xvar = "y", yvar = yearlyVariables)
#   
#   saveGgPlot(aplot, saveGgPlot,  
#              plot_height = 960, plot_width = NULL, 
#              ratio = 4/3, 
#              scale = 1, fileName = name_simu, fileSuffix = ".pdf")
#   rm(aplot)
#   rm(name_simu)
# }
# 
# for(i in 1:length(simuListRU100_st)){
#   name_simu = names(simuListRU100_st)[i]
#   aplot = simuPlots(simuListRU100_st[[i]], tableName = "y", xvar = "y", yvar = yearlyVariables)
#   
#   saveGgPlot(aplot, saveGgPlot,  
#              plot_height = 960, plot_width = NULL, 
#              ratio = 4/3, 
#              scale = 1, fileName = name_simu, fileSuffix = ".pdf")
#   rm(aplot)
#   rm(name_simu)
# }
