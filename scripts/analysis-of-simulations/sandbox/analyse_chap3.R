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


yvarList = c("GPP", "TR", "REWmin", "LAImaxThisYear", "dbh", "RU_shortage_max")
standYearTable_RU50 = makeStandYearTableFromStandScale_universal(simuListRU50_st, yvarList)
standYearTable_RU100 = makeStandYearTableFromStandScale_universal(simuListRU100_st, yvarList)

standYearTable = rbind(standYearTable_RU50, standYearTable_RU100)
standYearTable$RU = getRUfromCode_site(standYearTable$code_site)
standYearTable$RU_shortage_max_relative = standYearTable$RU_shortage_max / standYearTable$RU
standYearTable$classSize = getClassSizefromCode_site(standYearTable$code_site)
standYearTable$composition = getCompositionfromCode_site(standYearTable$code_site)
standYearTable$compositionIndex = getCompositionIndexfromCode_site(standYearTable$code_site)
standYearTable$isMixed = standYearTable$compositionIndex %in% c(2,3)

standYearTable$beechCompositionRate = ifelse(standYearTable$composition == "HETpur", 1, 
                                             ifelse(standYearTable$composition == "HETsap", 0.70, 
                                                    ifelse(standYearTable$composition == "SAPhet", 0.30, 
                                                           ifelse(standYearTable$composition == "SAPpur", 0,NA))))


standYearTable$nTree = getNTreefromCode_site(standYearTable$code_site)
standYearTable$nha = getNhafromCode_site(standYearTable$code_site)
standYearTable$standArea_m2 = standYearTable$nTree / standYearTable$nha *10000
colnames(standYearTable)[colnames(standYearTable) == "LAImaxThisYear"] = "LAI"
standYearTable$gha = standYearTable$nha * pi * (standYearTable$dbh/100 / 2)**2
standYearTable$LAI = round(standYearTable$LAI * 100) * 0.01
standYearTable$GPPperLAI = standYearTable$GPP / standYearTable$LAI




# One table for all years
standYearTable_meanYear = aggregate(standYearTable, FUN = mean, by = list(code_site_bis = standYearTable$code_site))
standYearTable_meanYear$code_site = standYearTable_meanYear$code_site_bis
standYearTable_meanYear$code_site_bis = NULL
standYearTable_meanYear$RU = getRUfromCode_site(standYearTable_meanYear$code_site)
standYearTable_meanYear$RU_shortage_max_relative = standYearTable_meanYear$RU_shortage_max / standYearTable_meanYear$RU
standYearTable_meanYear$classSize = getClassSizefromCode_site(standYearTable_meanYear$code_site)
standYearTable_meanYear$composition = getCompositionfromCode_site(standYearTable_meanYear$code_site)
standYearTable_meanYear$compositionIndex = getCompositionIndexfromCode_site(standYearTable_meanYear$code_site)
standYearTable_meanYear$isMixed = standYearTable_meanYear$compositionIndex %in% c(2,3)
standYearTable_meanYear$nTree = getNTreefromCode_site(standYearTable_meanYear$code_site)
standYearTable_meanYear$nha = getNhafromCode_site(standYearTable_meanYear$code_site)
standYearTable_meanYear$GPPperLAI = standYearTable_meanYear$GPP / standYearTable_meanYear$LAI


# Plot organization of data

# VARIABLES = nTree (LAI) x classSize x RU x composition
# nTree (3) x classSize (3) x RU (2) x composition (4) = 72

# LAI depends solely on number of tree (not on classSize, nor composition) 
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = LAI, x = nTree, size = gha)) + geom_point() + facet_grid(classSize ~ composition)
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = LAI, x = gha, size = nha)) + geom_point() + facet_grid(classSize ~ composition)

# gha is the same for each composition, gha depends on density (nTree) and class size
ggplot(subset(standYearTable_meanYear, year == standYearTable_meanYear$year[1] & RU == 50), aes(y = gha, x = nTree, size = classSize)) + geom_point() + facet_grid( ~ composition)



# Plots basic knowledge

# LAI / nTree drives GPP and TR.
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = as.factor(nTree))) + geom_point() + ylim(c(0, NA))

# ClassSize has a negligeable effect
ggplot(standYearTable_meanYear, aes(x = LAI, y = GPP, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + ylim(c(0, NA)) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + ylim(c(0, NA)) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
ggplot(standYearTable_meanYear, aes(x = LAI, y = RU_shortage_max, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + ylim(c(0, NA)) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)
ggplot(standYearTable_meanYear, aes(x = LAI, y = REWmin, color = classSize, size = 4 -as.numeric(as.factor(classSize)))) + geom_point() + ylim(c(0, NA)) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) + guides(size = FALSE)

# Same plot, but with code_site shown and saved as file
myplot = ggplot(standYearTable_meanYear, aes(x = LAI, y = TR, color = classSize, size = 4 -as.numeric(as.factor(classSize)), label = code_site)) + geom_point() + ylim(c(0, NA)) + facet_grid(RU ~ composition) + scale_color_manual(values=c("black", "grey", "white")) +
  geom_text_repel(size=2, alpha=0.6, color = "black") + guides(size = FALSE)
saveGgPlot(myplot, plot_height = 1080, plot_width = 1920, scale = 1, plotfolderPlot = folderPlot, fileName = "TR")

# hydric stress (low RU) increase maximum water shortage
ggplot(standYearTable_meanYear, aes(x = as.factor(RU), y = RU_shortage_max_relative, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + ylim(c(0, NA))


# Plot advanced results
ggplot(standYearTable_meanYear, aes(x = composition, y = GPPperLAI, color = as.factor(RU))) + geom_point() + ylim(c(0, NA)) + facet_grid(classSize ~ nTree)

ggplot(standYearTable_meanYear, aes(x = composition, y = TR, color = as.factor(RU))) + geom_point() + ylim(c(0, NA))  + facet_grid(classSize ~ nTree)




# Compute Net Biodiversity Effects (NBE)
standYearTable_meanYear$NBE_gpp = 0
standYearTable_meanYear$NBE_tr = 0
standYearTable_meanYear$NBE_ru_shortage_max = 0
standYearTable_meanYear$NBE_ru_shortage_max_relative = 0
standYearTable_meanYear$NBE_rewMin = 0
for(i in 1:dim(standYearTable_meanYear)[1]){
  aSimu = standYearTable_meanYear[i, ]
  if(!grepl(aSimu$composition, pattern = "pur")){
    nTreeSimu = aSimu$nTree
    classSizeSimu = aSimu$classSize
    RUSimu = aSimu$RU
    
    beechMono = subset(standYearTable_meanYear, nTree == nTreeSimu & classSize == classSizeSimu & RU == RUSimu & composition == "HETpur")
    firMono = subset(standYearTable_meanYear, nTree == nTreeSimu & classSize == classSizeSimu & RU == RUSimu & composition == "SAPpur")
    
    proportion = aSimu$beechCompositionRate
    
    numericalColumns = colnames(aSimu)[vapply(aSimu[1,], FUN = is.numeric, FUN.VALUE = FALSE, USE.NAMES = F)]
    mixedAdditiveASimu = beechMono[, numericalColumns] * proportion + firMono[, numericalColumns] * (1-proportion) 
    
    standYearTable_meanYear[i, ]$NBE_gpp = (standYearTable_meanYear[i, ]$GPP - mixedAdditiveASimu$GPP) / mixedAdditiveASimu$GPP
    standYearTable_meanYear[i, ]$NBE_tr = (standYearTable_meanYear[i, ]$TR - mixedAdditiveASimu$TR) / mixedAdditiveASimu$TR
    standYearTable_meanYear[i, ]$NBE_ru_shortage_max = (standYearTable_meanYear[i, ]$RU_shortage_max - mixedAdditiveASimu$RU_shortage_max) / mixedAdditiveASimu$RU_shortage_max
    standYearTable_meanYear[i, ]$NBE_ru_shortage_max_relative = (standYearTable_meanYear[i, ]$RU_shortage_max_relative - mixedAdditiveASimu$RU_shortage_max_relative) / mixedAdditiveASimu$RU_shortage_max_relative
    standYearTable_meanYear[i, ]$NBE_rewMin = (standYearTable_meanYear[i, ]$REWmin - mixedAdditiveASimu$REWmin) / mixedAdditiveASimu$REWmin
    
  }
}

# H1.1 et H1.2 : Evaluation of NBE on growth and transpiration
# NBE is globally positive on growth and transpiration. NBE growth is more pronounced on Beech-dominant mixture, NBE transpiration is more pronounced on fir-dominant mixture
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_tr, y = NBE_gpp, color = composition)) + geom_point() + xlim(c(0, NA)) + ylim(c(0, NA))
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_tr, y = NBE_gpp, color = composition, size = as.factor(RU))) + geom_point() + xlim(c(0, NA)) + ylim(c(0, NA))


# NBE is not exactly identical on transpiration and RU shortage max
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = NBE_ru_shortage_max, y = NBE_tr, color = composition)) + geom_point() + xlim(c(0, NA)) + ylim(c(0, NA))


# H2.1.a La densité augmente le NBE-growth sur les fir-dominant, et réduit le NBE-growth sur les beech-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_gpp, color = composition)) + geom_boxplot() + ylim(c(0, NA))
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_gpp, color = composition)) + geom_boxplot() + facet_grid( . ~ classSize) + ylim(c(0, NA))


# H2.1.b La densité augmente le NBE-transpiration sur tous les types de peuplements
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_tr, color = composition)) + geom_boxplot() + ylim(c(0, NA))
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(nTree), y = NBE_tr, color = composition)) + geom_boxplot() + facet_grid( . ~ classSize) + ylim(c(0, NA))


# H2.2.a L'aĝe a un effet positif sur la NBE-gpp sur les beech-dominant et négatif sur la NBE-gpp sur les fir-dominant
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_gpp, color = composition)) + geom_boxplot() + ylim(c(0, NA))
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_gpp, color = composition)) + geom_boxplot() + facet_grid( . ~ nTree) + ylim(c(0, NA))


# H2.2.b L'aĝe a un effet nul sur la NBE-water, sauf pour les petites densité (effet positif) et grosse densité (effet négatif chez les beech-dominant)
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_tr, color = composition)) + geom_boxplot() + ylim(c(0, NA))
# bonus
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = classSize, y = NBE_tr, color = composition)) + geom_boxplot() + facet_grid( . ~ nTree) + ylim(c(0, NA))


# H3.1 Le stress hydrique réduit la NBE growth
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_gpp, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + ylim(c(0, NA))

# H3.2 le stress hydrique réduit la NBE tr (nivellement)
ggplot(subset(standYearTable_meanYear, isMixed), aes(x = as.factor(RU), y = NBE_tr, color = composition)) + geom_boxplot() + facet_grid( . ~ .) + ylim(c(0, NA))
# --> comment rendre compte de l'effet du stress hydrique sur les dégats du peuplement ?


