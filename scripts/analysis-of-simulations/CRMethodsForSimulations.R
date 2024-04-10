# This work © 2024 by Camille Rouet is licensed under [CC BY-NC 4.0](http://creativecommons.org/licenses/by-nc/4.0/).
# Author: Camille Rouet 2021-2024
# Methods of this file are dedicated to the analysis of castaneaOnly and physiodemogenetics simulations
# The input files of theses methods are the log of the castanea library found in var

if(! "PROGRAM_ON_SERVER" %in% ls()){
  PROGRAM_ON_SERVER = FALSE
}

# library(tidyverse) # Group of packages that contains ggplot2, readr...
library(readr)
library(data.table)
library(ggplot2)
library(cowplot)
library(tibble)
library(dplyr)
library(Metrics)
library(corrplot)
library(ggcorrplot)
library(ggrepel)
library(units)
library(tidyr)
library(stringr)
library(ggpubr)

if(!PROGRAM_ON_SERVER){
  library(tcltk2)
  if(system.file(package='beepr') != ""){ library(beepr) }
}

if(system.file(package='rmote') != ""){ library(rmote) }



# unités
WtoPhoton = 4.54 # 1 W = 4.54 micromol of photo per sec

logList1 = c("yearlyResults", "fractionMaps", "standFractionMaps", "SLyearlyLightResult", "SLenergyAttributionTrees", "SLenergyAttributionCells", "SLenergyBalanceCheck") ;
logList2 = c(logList1, "dailyResults", "soilWater", "woodGrowth", "radiationDaily", "treesDailyResults", "standDailyResults")
logList3 = c(logList2, "hourlyResults", "radiationAbsorbedCanopy", "thermic", "inputRadiation", "radiationSoil")
logList4 = c(logList3, "radiationAbsorbedLayer", "photosynthesis", "radiationAbsorbedLayerPerLine")
logListList = list(logList1, logList2, logList3, logList4) ;
rm(logList1) ; rm(logList2) ; rm(logList3) ; rm(logList4)

# variables on $dailyResults
initYearlyVariables = c("year", "idFmCell", "Xposition", "Yposition", "species", "age", "rw", "dbhInitSimulation", "BAI", "height", "LAImaxThisYear")
yearlyVariables = c("dbh", "GPP", "NPP", "LAImaxThisYear", "BiomassOfReserves", "BAI")

globalDailyVariables = c("GPP", "BSS", "dailyRemainingCarbone", "LAIday", "NEE", "REW",  "TR", "ETRsoil")
globalDailyVariables2 = c("GPP", "BSS", "dailyRemainingCarbone", "LAIday", "LMAmoy", "REW",  "TR", "ETRsoil")
carbonDailyVariables = c("GPP", "BSS", "dailyRemainingCarbone", "leafGrowth", "NEE", "RVM", "RVC", "Rheterotrophic")
carbonDailyVariables2 = c("GPP", "BSS", "reservesMortality", "reservesGrowth", "LAIday", "dailyRemainingCarbone", "leafGrowth", "woodGrowth", "RVC", "RVM")
waterDailyVariables = c("LAIday", "RU_currentLevel", "REW", "TR", "ETRsoil", "ETRcan", "drainage", "Rain")

# variables on $woodGrowth
woodDailyVariables =  c( "BF", "biomassOfReserves", "biomassOfFineRoot", "woodGrowth", "biomassOfTrunk", "RVM")
woodDailyVariables2 =  c( "BF", "biomassOfReserves", "biomassOfFineRoot", "woodGrowth", "biomassOfTrunk", "RVM", "dailyRemainingCarbone", "LAIday")
biomassVariables = c('biomassOfTrunk', 'biomassOfBranch', 'biomassOfCoarseRoot', 'biomassOfFineRoot', 'biomassOfReserves', 'BF')

# variables on radiationDaily
radiationDailyVariables = c("Rg", "rain", "Tmoy", "lightCompetitionIndex", "vegPAR", "vegPIR", "soilPAR", "soilPIR")
radiationDailyVariablesPARPIRdirdiff = c("vegPARdir","vegPARdif", "vegPIRdir",  "vegPIRdif", "soilPARdir","soilPARdif", "soilPIRdir",  "soilPIRdif") # MJ / m2
radiationDailyVariablesPAR = c("inciPARdir", "vegPARdir", "soilPARdir", "inciPARdif","vegPARdif","soilPARdif") # MJ / m2
radiationDailyVariablesPIR = c("inciPIRdir", "vegPIRdir", "soilPIRdir", "inciPIRdif","vegPIRdif","soilPIRdif") # MJ / m2
radiationDailyVariablesDir = c("inciPARdir", "vegPARdir", "soilPARdir", "inciPIRdir", "vegPIRdir", "soilPIRdir") # MJ / m2
radiationDailyVariablesDif = c("inciPARdif","vegPARdif","soilPARdif", "inciPIRdif","vegPIRdif","soilPIRdif") # MJ / m2
radiationDailyVariablesGlobDirDiff = c("inciTot", "vegTot", "soilTot", "inciDir", "vegDir", "soilDir", "inciDif", "vegDif", "soilDif") # MJ / m2
radiationDailyVariablesGlobPARPIR = c("inciTot", "vegTot", "soilTot", "inciPAR", "vegPAR", "soilPAR", "inciPIR", "vegPIR", "soilPIR") # MJ / m2

climateDailyVariables = c("Rg", "rain", "Tmoy", "Tmin", "Tmax") # to use with $dailyMeteo

commonBaseVariablesYearly = c("year", "idFmCell", "species")
commonBaseVariablesDaily = c(commonBaseVariablesYearly, "day")






# ..............................................................................
# IMPORT -----------------------------------------------------------------------


# this method join the table of a tmpListTable, a list of table
# tableName1 and tableName2 are two data.frame of tmpListTable that should share the same structure (eg, one row per cell per day )
joinTable = function(tmpListTable, tableName1, tableName2){
  joinColumns = colnames(tmpListTable[[tableName2]])[ colnames(tmpListTable[[tableName2]]) %in% colnames(tmpListTable[[tableName1]]) ]
  
  # check compatibility
  valid = FALSE
  tab1 = tmpListTable[[tableName1]][, joinColumns]
  tab2 = tmpListTable[[tableName2]][, joinColumns]
  
  if(sum(dim(tab1) - dim(tab2)) == 0){
    sumDiffJoinColumns = sum(tab1 - tab2)
    if(sumDiffJoinColumns == 0){
      valid = TRUE
    }
  }
  
  # check table compatibility
  if(valid){
    tmpListTable[[tableName1]] = inner_join(tmpListTable[[tableName1]], tmpListTable[[tableName2]], by = joinColumns)
  }else{
    print(paste0("Join columns data are not matching for ", tableName1, " and ", tableName2) )
  }
  
  return(tmpListTable)
}



# Import simulation files (e.g. from the var folder of CAPSIS)
# Example : to import castanea_dailyResults.log and castanea_hourlyResults.log
# which are located in "~capsis/var", the formula is :
# > CAST = importSimu(folderPath = "~capsis/var", prefix = "castanea_", logList = c("dailyResults", "hourlyResults"))
# The result is a list of table, accessed by, eg, "CAST$dailyResults" and "CAST$hourlyResults".
# logList : vector containing the desired table name to import
# prefix : prefix to all log files to import
importSimu = function(folderPath, logList = NULL, prefix = "", fileExtention = ".log", 
                      printMessage = TRUE, shortMessage = FALSE, years = NULL){
  
  if(!file.exists(folderPath)){
    stop(paste0("Folder ", folderPath, " does not exists."))
  }
  
  # yearlyResultsFile = paste0(folderPath, prefix, "yearlyResults", fileExtention)
  # if(!file.exists( yearlyResultsFile )){
  #   stop(paste0("yearlyResults not found at ", yearlyResultsFile))
  # }
  
  if(is.null(logList)){
    logList = logListList[[length(logListList)]]
  }
  
  if(printMessage){
    cat("\nImporting :\t", prefix, "\n")
    if(!shortMessage){
      cat("from :\t\t", folderPath, "\n")
    }
  }
  
  
  tmp = list()
  logExist = c()
  
  # import all log files
  for(i in 1:length(logList)){
    suffixName = logList[i]
    suffix = paste0(logList[i], fileExtention)
    path = paste0(folderPath, prefix, suffix)
    
    if(file.exists(path)){
      logExist = c(logExist, suffixName)
      tmp[[suffixName]] = read_delim(path, delim = ";", col_types = cols())
    }
  }
  
  
  if("SLenergyAttributionTrees" %in% logExist){
    tmp[["SLenergyAttributionTrees"]]$idFmCell = tmp[["SLenergyAttributionTrees"]]$idTree
  }
  
  # check if dailyResults is empty
  if("dailyResults" %in% logExist){
    if(dim(tmp[["dailyResults"]])[1] < 1){
      print("Warning : dailyResults does not contain any line")
    }
    
    # automaiccaly
    # tmp[["dailyResults"]] = pixellateVariableOnDailyTable(tmp[["dailyResults"]], variableName = "Rain", categorySize = 20)
  }
  
  # check if hourlyResults is empty
  if("hourlyResults" %in% logExist){
    if(dim(tmp[["hourlyResults"]])[1] < 1){
      print("Warning : hourlyResults does not contain any line")
    }
  }
  
  # # join daily variables
  # if("dailyResults" %in% logExist & "dailyMeteo" %in% logExist){
  #   tmp = joinTable(tmp, "dailyResults", "dailyMeteo")
  # }
  
  # join daily variables
  if("dailyResults" %in% logExist & "woodGrowth" %in% logExist){
    tmp = joinTable(tmp, "dailyResults", "woodGrowth")
  }
  
  # join daily variables
  if("dailyResults" %in% logExist & "radiationDaily" %in% logExist){
    
    tmp = joinTable(tmp, "dailyResults", "radiationDaily")
    
    # Extract last day results
    tmp[["dailyResultsLastDay"]] = subset(tmp$dailyResults, (year%%4 == 0 & day == 366) | (year%%4 != 0 & day == 365))
    
    # add BSSLastDay to yearlyResults
    # apply to each line the concatenation of $year and $idFmCell to have year_idFmCell
    year_id_dResultsLastDay = apply(tmp$dailyResultsLastDay[, c("year", "idFmCell")], MARGIN = 1, FUN = function(x) paste(x, collapse = "_"))
    year_id_yResults = apply(tmp$yearlyResults[, c("year", "idFmCell")], MARGIN = 1, FUN = function(x) paste(x, collapse = "_"))
    
    # use match() to match correct idFmcell x year lines
    match(year_id_yResults, year_id_dResultsLastDay)
    
    tmp[["yearlyResults"]]$BSSLastDay = tmp[["dailyResultsLastDay"]]$BSS[match(year_id_yResults, year_id_dResultsLastDay)]
  }
  
  
  # join daily variables
  if("dailyResults" %in% logExist & "treesDailyResults" %in% logExist){
    tmp = joinTable(tmp, "dailyResults", "treesDailyResults")
  }
  
  # join hourly variables
  if("hourlyResults" %in% logExist & "radiationAbsorbedCanopy" %in% logExist){
    tmp = joinTable(tmp, "hourlyResults", "radiationAbsorbedCanopy")
  }
  
  # join hourly variables
  if("hourlyResults" %in% logExist & "radiationAbsorbedLayerPerLine" %in% logExist){
    tmp = joinTable(tmp, "hourlyResults", "radiationAbsorbedLayerPerLine")
  }
  
  
  
  # import inventory
  tmp[["inventory"]] = suppressWarnings( importInventory(folderPath, suffix = "inventoryCopy", prefix = prefix))
  # add a climate file name entry
  tmp[["inventory"]]$climateFileName = unique(tmp[["yearlyResults"]]$climateFileName)
  
  
  # import settings log
  fmSettingsPath = paste0(folderPath, prefix, "fmSettings", fileExtention)
  if(file.exists(fmSettingsPath)){
    tmp[["fmSettings"]] = read_delim(file = fmSettingsPath, delim = ";", col_types = cols())
  }
  
  # PDG settings
  PDGInitialParametersPath = paste0(folderPath, prefix, "PDGInitialParameters", fileExtention)
  if(file.exists( PDGInitialParametersPath) ){
    tmp[["PDGInitialParameters"]] = read_delim(file = PDGInitialParametersPath, delim = ";", col_types = cols())
  }
  
  
  # import soil inventory (of the first cell only)
  soilLine = importSoilForFirstCell(folderPath, suffix = "inventoryCopy", prefix = prefix)
  if(!is.null(soilLine)){
    tmp[["soilFirstCell"]] = soilLine
  }
  
  # filter on year
  if(!is.null(years)){
    tmp = sfiltering(tmp, year = years)
  }
  
  if(printMessage){
    message = ""
    if("PDGInitialParameters" %in% names(tmp)){
      option_1 = "pdgArenaTagMode"
      option_2 = "pdgArenaDisableWaterInteraction"
      if("pdgLightTagMode" %in% names(tmp$PDGInitialParameters)){
        # old denomination
        option_1 = "pdgLightTagMode"
        option_2 = "pdgLightDisableWaterInteraction"
      }
      
      message1 = ifelse(tmp$PDGInitialParameters[[option_1]], "SL tag mode", "SL simple mode")
      message2 = ifelse(tmp$PDGInitialParameters[[option_2]], "without water interactions", "with water interactions")
      message = paste0("\t\tOptions: ", message1, ", ", message2)
    }
    cat("Simu date :\t", tmp$inventory$dateTime, message,  "\n")
    cat("Inventory :\t", tmp$inventory$originalinventory, "\n")
  }
  
  if(prefix != ""){
    tmp$name = prefix
  }
  
  return(tmp)
} # end importSimu



# import the inventory of a given simulation
importInventory = function(folderPath, suffix, prefix = "", fileExtention = ".log"){
  filePath = paste0(folderPath, prefix, suffix, fileExtention)
  tmp = read_delim(filePath, delim = "=", comment = "#", col_names = FALSE, trim_ws = TRUE, col_types = cols())
  
  # remove NAs
  linesToRemove = c()
  for(i in 1:length(tmp$X1)){
    if(is.na(tmp$X1[i]) | is.na(tmp$X2[i]) ){
      linesToRemove = c(linesToRemove, i)
    }
  }
  tmp = tmp[-linesToRemove,]
  
  # convert to table in row
  inventory = data.frame(matrix(nrow = 1, ncol = length(tmp$X1)))
  varNames = tolower(tmp$X1)
  names(inventory) = varNames
  inventory[1,] = tmp$X2
  
  # turn string to numeric
  for(i in 1:length(inventory)){
    if(!is.na(as.numeric(inventory[i]))){
      inventory[i] = as.numeric(inventory[i])
    }
  }
  
  # standArea for PDG Arena
  if("ncol" %in% colnames(inventory) & "nlin" %in% colnames(inventory) & "cellwidth" %in% colnames(inventory)){
    inventory$standArea_m2 = inventory$ncol * inventory$nlin * inventory$cellwidth**2 
  }
  
  
  colInv = colnames(inventory)
  inventorySplit1 = inventory[, colInv[1:which(colInv == "time")]]
  inventorySplit2 = inventory[, colInv[(which(colInv == "time") +1) : length(colInv)]]
  inventorySplit1$dateTime = paste0(inventory$date, " ", inventory$time)
  inventory = cbind(inventorySplit1, inventorySplit2)
  
  return(inventory)
}


# extract the soil information from the first cell found in inventory
importSoilForFirstCell = function(folderPath, suffix, prefix = "", fileExtention = ".log"){
  
  filePath = paste0(folderPath, prefix, suffix, fileExtention)
  
  tmp = read_delim(filePath, delim = "nodelim", col_names = FALSE, trim_ws = TRUE, col_types = cols())
  
  lenCellRecord = -1
  cellRecordFound = FALSE
  cellRecordPass = FALSE
  i = 1
  while(!cellRecordFound & !cellRecordPass){
    # read each line
    line = strsplit( paste0(tmp[i, ]), "\t")[[1]]
    
    # when cells are found (PDG)
    if(line[1] == "#cID"){
      soilVarNamesLine = line
      cellRecordFound = TRUE
      lenCellRecord = length(soilVarNamesLine)
      soilVarNames = soilVarNamesLine[4:lenCellRecord]
      
      soilVariables = strsplit( paste0(tmp[i+1, ]), "\t")[[1]][4:lenCellRecord]
      soilVariables = as.numeric(soilVariables)
    }
    
    # when cells are found (CASTANEA)
    if(line[1] == "#Type"){
      soilVarNamesLine = line
      cellRecordFound = TRUE
      lenCellRecord = length(soilVarNamesLine)
      soilVarNames = soilVarNamesLine[7:23]
      soilVarNames[soilVarNames == "stone"] = "stoneContent"
      soilVarNames[soilVarNames == "soilHeight"] = "solHeight"
      
      soilVariables = strsplit( paste0(tmp[i+1, ]), "\t")[[1]][7:23]
      soilVariables = as.numeric(soilVariables)
    }
    
    i = i +1
    if(i > dim(tmp)[1]){
      print("Cell line(s) were not found in inventory. Maybe the first column name has changed.")
      cellRecordPass = TRUE
    }
  }
  
  
  if(cellRecordFound){
    soilFirstCell = data.frame(t(soilVariables))
    colnames(soilFirstCell) = soilVarNames
    soilFirstCell$RU = (soilFirstCell$wfc - soilFirstCell$wilt) * (1 - soilFirstCell$stoneContent) * soilFirstCell$solHeight * soilFirstCell$bulk
    soilFirstCell$wfcMinusWilt = soilFirstCell$wfc - soilFirstCell$wilt
    
    # rearranging columns
    columns = colnames(soilFirstCell)
    posWilt = which(columns == "wilt")
    lastPos = length(columns)
    rearrangedColumnsPosition = c(1:posWilt, lastPos, (posWilt+1):(lastPos-1))
    
    
    soilFirstCell = soilFirstCell[, rearrangedColumnsPosition]
  }
  
  return(soilFirstCell)
  
} # end importSoilForFirstCell


### import climate in Safran
importClimateFromSimu = function(asimu){
  # climate file name
  climatePathOfSimu = unique(asimu$yearlyResults$climateFileName)
  climateFile = strsplit(climatePathOfSimu, split = "/")[[1]]
  climateFile = climateFile[length(climateFile)]
  
  # climate file folder
  climateFolder = paste0(capsisPath, "data/castaneaonly/climate/safran/")
  
  climatePath = paste0(climateFolder, climateFile)
  
  # check number of line to skip before headers
  climateComment = suppressWarnings(read_delim(climatePath, delim = "\t", n_max = 100)[[1]])
  nCommentLine = sum(sapply(climateComment, function(x) grepl(pattern = '#', x)))
  
  # import
  climateTable = suppressWarnings( read_delim(climatePath, delim = "\t", skip = nCommentLine+1) )
  colnames(climateTable)[colnames(climateTable) == "h"] = "hour"
  colnames(climateTable)[colnames(climateTable) == "m"] = "month"
  colnames(climateTable)[colnames(climateTable) == "d"] = "day"
  colnames(climateTable)[colnames(climateTable) == "# y"] = "year"
  
  print(paste0("Climate ", climateFile, " war imported from folder ", climateFolder))
  return(climateTable)
}


### import climate in Safran
importClimate = function(climatePath){
  
  climatePath = paste0(climateFolder, climateFile)
  
  # check number of line to skip before headers
  climateComment = suppressWarnings(read_delim(climatePath, delim = "\t", n_max = 100)[[1]])
  nCommentLine = sum(sapply(climateComment, function(x) grepl(pattern = '#', x)))
  
  # import
  climateTable = suppressWarnings( read_delim(climatePath, delim = "\t", skip = nCommentLine+1) )
  colnames(climateTable)[colnames(climateTable) == "h"] = "hour"
  colnames(climateTable)[colnames(climateTable) == "m"] = "month"
  colnames(climateTable)[colnames(climateTable) == "d"] = "day"
  colnames(climateTable)[colnames(climateTable) == "# y"] = "year"
  
  print(paste0("Climate ", climateFile, " war imported from folder ", climateFolder))
  return(climateTable)
}




# Convert the hourly climate into a daily climate
convertClimateHourlyToDaily = function(climateHourly){
  
  # old way
  # climateDaily = tibble(year = 0, month = 0, day = 0, gr = 0, tmoy = 0, tmin = 0, tmax = 0, rh = 0, rhmin = 0, ws = 0, p = 0, .rows = 0)
  # lastLine = 0
  # 
  # nLine = dim(climateHourly)[1]
  # 
  # pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
  #                      max = nLine, # Maximum value of the progress bar
  #                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
  #                      width = 50,   # Progress bar width. Defaults to getOption("width")
  #                      char = "=")   # Character used to create the bar
  # 
  # for(i in 1:nLine){
  #   ayear = climateHourly[i, ]$year
  #   amonth = climateHourly[i, ]$month
  #   aday = climateHourly[i, ]$day
  #   
  #   if(i > 1) {
  #     if(aday == climateDaily[lastLine, ]$day){
  #       next
  #     }
  #   }
  #   
  #   climateSubsetForDay = subset(climateHourly, day == aday & year == ayear)
  #   
  #   climateDaily = bind_rows(climateDaily, 
  #                            tibble(year = ayear, month = amonth, day = aday, 
  #                                   gr = mean(climateSubsetForDay$gr), 
  #                                   tmoy = mean(climateSubsetForDay$t), 
  #                                   tmin = min(climateSubsetForDay$t), 
  #                                   tmax = max(climateSubsetForDay$t), 
  #                                   rh = mean(climateSubsetForDay$rh), 
  #                                   rhmin = min(climateSubsetForDay$rh), 
  #                                   ws = mean(climateSubsetForDay$ws), 
  #                                   p = sum(climateSubsetForDay$p)
  #                                   ))
  #   lastLine = lastLine + 1
  #   setTxtProgressBar(pb, i) # change progress bar
  #   
  # }
  
  
  dayAndYear = unique(climatVl[, c("year", "month", "day")])
  
  climateDaily = tibble(year = 0, month = 0, day = 0, gr = 0, tmoy = 0, tmin = 0, tmax = 0, rh = 0, rhmin = 0, ws = 0, p = 0, .rows = dim(dayAndYear)[1])
  climateDaily$year = dayAndYear$year
  climateDaily$day = dayAndYear$day
  climateDaily$month = dayAndYear$month
  
  meanColumns = c("gr", "t", "rh", "ws")
  minColumns = c("t", "rh")
  maxColumns = c("t")
  sumColumns = c("p")
  
  meanColumnsNames = c("gr", "tmoy", "rh", "ws")
  minColumnsNames = c("tmin", "rhmin")
  maxColumnsNames = c("tmax")
  sumColumnsNames = c("p")
  
  climateDaily[, meanColumnsNames] = aggregate(x = climateHourly[, meanColumns], FUN = mean, by = list(days = climateHourly$day, years = climateHourly$year))[, meanColumns]
  climateDaily[, minColumnsNames] = aggregate(x = climateHourly[, minColumns], FUN = min, by = list(days = climateHourly$day, years = climateHourly$year))[, minColumns]
  climateDaily[, maxColumnsNames] = aggregate(x = climateHourly[, maxColumns], FUN = max, by = list(days = climateHourly$day, years = climateHourly$year))[, maxColumns]
  climateDaily[, sumColumnsNames] = aggregate(x = climateHourly[, sumColumns], FUN = sum, by = list(days = climateHourly$day, years = climateHourly$year))[, sumColumns]
  
  # close(pb) # Close progress bar
  
  return(climateDaily)
}




# Call a subset function on all sub datatable of dataSimu to extract only a selection of year and fmCell
# dataSimu should be a list of data table resulting from importSimu().
# It can also be directly a data table with fmCellId, species, year or k column
# It can also be a list of simulation (a list of lists of table)
# firstYear as TRUE in order to select only the first year of simulation
sfiltering = function(dataSimu, fmCellIDs = NULL, fmSpeciesIDs = NULL, years = NULL, layers = NULL, firstYear = FALSE){
  
  # if dataSimu is simply null..
  if(is.null(dataSimu)){
    return(dataSimu)
  }
  
  # if dataSimu is a list of simulation
  if("yearlyResults" %in% names(dataSimu[[1]])){
    for(simuName in names(dataSimu)){
      dataSimu[[simuName]] = sfiltering(dataSimu[[simuName]], fmCellIDs, fmSpeciesIDs, years, layers, firstYear)
    }
    return(dataSimu)
  }
  
  # if dataSimu is not a list of tab, but a tab directly
  if(is.data.frame(dataSimu)){
    
    if(is.null(years)){
      years = unique(dataSimu$year)
    }
    
    if(firstYear){
      years = unique(dataSimu$year)[1]
    }
    
    if(is.null(fmCellIDs)){
      fmCellIDs = unique(dataSimu$idFmCell)
    }
    if(is.null(fmSpeciesIDs)){
      fmSpeciesIDs = unique(dataSimu$species)
    }
    
    res = subset(dataSimu, year %in% years & idFmCell %in% fmCellIDs & species %in% fmSpeciesIDs )
    return(res)
  }
  
  # Species filter
  if(!is.null(fmSpeciesIDs)){
    
    if(is.null(fmCellIDs)){
      fmCellIDs = unique(dataSimu$yearlyResults$idFmCell)
    }
    
    fmCellIDs_speciesCompatible = unique(dataSimu$yearlyResults[dataSimu$yearlyResults$species %in% fmSpeciesIDs, ]$idFmCell)
    
    # filter for fmSpeciesIDs AND fmCellsIDs (if not given, not filter on fmCellIDs)
    fmCellIDs = fmCellIDs[fmCellIDs %in% fmCellIDs_speciesCompatible]
  }
  
  res = dataSimu # make a clone (in R, this makes a total clone, it is not just a new pointer)
  
  
  
  
  
  # reduce data set for each data table
  for(tabName in names(res)){
    if(! tabName %in% c("inventory", "fmSettings", "PDGInitialParameters", "soilFirstCell")){
      
      # idFmCell selection
      
      if(!is.null(fmCellIDs)){
        subdata = res[[tabName]]
        if("idFmCell" %in% names(subdata)){
          res[[tabName]] = subset(subdata, subdata$idFmCell %in% fmCellIDs)
        }
      }
      
      # if first year
      if(firstYear){
        years = unique(dataSimu$yearlyResults$year)[1]
      }
      
      # year selection
      if(!is.null(years)){
        subdata = res[[tabName]]
        if("year" %in% names(subdata)){
          res[[tabName]] = subset(subdata, subdata$year %in% years)
        }
      }
      
      # layer selection
      if(!is.null(layers)){
        subdata = res[[tabName]]
        if("k" %in% names(subdata)){
          res[[tabName]] = subset(subdata, subdata$k %in% layers)
        }
      }
    }
  }
  
  if(dim(res$yearlyResults)[1] == 0){
    res$isEmpty = TRUE
    warning(paste0("Filtering simulation ", res$name, " gives an empty simulation." ))
  }
  return(res)
}







# ..............................................................................
# MONOSP. COMBINAISON ----------------------------------------------------------

# return the proportion of basal area of a given species sp in the simulation data
# data can be a simulation list of a yearlyResults
# data should be 
getInitialSpeciesBasalAreaProportionInSimulation = function(data, targetSp){
  # verify data
  if(inherits(data, "list")){
    
    if(! "yearlyResults" %in% names(data)){
      stop("data is a list but does not contain yearlyResults")
    }
    yearlyResults = data$yearlyResults
    
  }else if(inherits(data, "data.frame")){
    
    if(!"year" %in% names(data)){
      stop("data is a data.frame but does not have a year columnn")
    }
    yearlyResults = data
    
  }else{
    stop("data is not a list of simulation nor a yearlyResults")
  }
  
  # verify the number of individuals
  if( length(unique(yearlyResults$idFmCell)) <= 1){
    stop("yearlyResults does not contains more than one individual.")
  }
  
  # verify targetSp
  if(! targetSp %in% yearlyResults$species){
    return(0)
  }
  
  
  firstYear = min(yearlyResults$year)
  yearlyResultsFirstYear = yearlyResults[yearlyResults$year == min(yearlyResults$year), ]
  yearlyResultsFirstYearTargetSp = yearlyResultsFirstYear[yearlyResultsFirstYear$species == targetSp, ]
  # yearlyResultsFirstYearOtherSp = yearlyResultsFirstYear[yearlyResultsFirstYear$species != targetSp, ]
  
  sumBasalArea = sum(pi * (yearlyResultsFirstYear$dbhInitSimulation/2)**2)
  sumBasalAreaTargetSp = sum(pi * (yearlyResultsFirstYearTargetSp$dbhInitSimulation/2)**2)
  
  proportionTargetSp = sumBasalAreaTargetSp / sumBasalArea
  
  return(proportionTargetSp)
}


# return the proportion of number of tree of a given species sp in the simulation data
# data can be a simulation list of a yearlyResults
# data should be 
getInitialSpeciesTreeProportionInSimulation = function(data, targetSp){
  # verify data
  if(inherits(data, "list")){
    
    if(! "yearlyResults" %in% names(data)){
      stop("data is a list but does not contain yearlyResults")
    }
    yearlyResults = data$yearlyResults
    
  }else if(inherits(data, "data.frame")){
    
    if(!"year" %in% names(data)){
      stop("data is a data.frame but does not have a year columnn")
    }
    yearlyResults = data
    
  }else{
    stop("data is not a list of simulation nor a yearlyResults")
  }
  
  # verify the number of individuals
  if( length(unique(yearlyResults$idFmCell)) <= 1){
    stop("yearlyResults does not contains more than one individual.")
  }
  
  # verify targetSp
  if(! targetSp %in% yearlyResults$species){
    return(0)
  }
  
  
  firstYear = min(yearlyResults$year)
  yearlyResultsFirstYear = yearlyResults[yearlyResults$year == min(yearlyResults$year), ]
  yearlyResultsFirstYearTargetSp = yearlyResultsFirstYear[yearlyResultsFirstYear$species == targetSp, ]
  # yearlyResultsFirstYearOtherSp = yearlyResultsFirstYear[yearlyResultsFirstYear$species != targetSp, ]
  
  nTree = length(unique(yearlyResultsFirstYear$idFmCell))
  nTreeSp = length(unique(yearlyResultsFirstYearTargetSp$idFmCell))
  
  proportionTargetSp = nTreeSp / nTree
  
  return(proportionTargetSp)
}




# Combine the value of regular monospecific (stand scale) simulation
# Result is a regular plurispecific (stand scale) simulation with no species interaction
# dataIrregPluriSp is the simulation from which data1 and data2 are monospecific adaptations
# proportionBA1 can be directly given instear of dataIrregPluriSp
combineStandScaleSimulations = function(data1, data2, proportionBA1 = NULL, proportionN1 = NULL, dataIrregPluriSp = NULL, outputDepth = 1){
  
  # if data1 is a list of indivScaleSimu
  if(!"yearlyResults" %in% names(data1) & "yearlyResults" %in% names(data1[[1]])){
    return( combineStandScaleSimulationsList(data1, data2, proportionBA1, proportionNtree1, dataIrregPluriSp, outputDepth) )
  }
  
  # if data are not stand scale
  if(length(unique(data1$yearlyResults$idFmCell)) > 1){
    stop("data1 is not stand scale.")
  }
  
  if(length(unique(data2$yearlyResults$idFmCell)) > 1){
    stop("data2 is not stand scale.")
  }
  
  if((is.null(proportionBA1) | is.null(proportionN1)) & is.null(dataIrregPluriSp)){
    stop("proportionBA1 and proportionN1 OR dataIrregPluriSp should be given.")
  } 
  
  if((!is.null(proportionBA1) | !is.null(proportionN1)) & !is.null(dataIrregPluriSp)){
    stop("proportionBA1 and proportionN1 OR dataIrregPluriSp should be given, not both.")
  } 
  
  if(is.null(proportionBA1)){
    proportionBA1 = getInitialSpeciesBasalAreaProportionInSimulation(dataIrregPluriSp, unique(data1$yearlyResults$species))
    proportionN1 = getInitialSpeciesTreeProportionInSimulation(dataIrregPluriSp, unique(data1$yearlyResults$species))
  }
  
  
  
  # Create a list of tibble (= the simulation object) with rows of 0
  dataRes = list()
  
  allTibNames = names(data1)
  
  
  # Create tabs for non-empty tabs and fill them with weighted averaged value
  for(i in 1:length(allTibNames)){
    tibName = allTibNames[i]
    tib = data1[[i]]
    tibColNames = colnames(tib)
    tibNVar = length(tibColNames)
    
    # skip some tables
    if (tibName %in% c("inventory", "name", "fmSettings", "PDGInitialParameters", "soilFirstCell")) {
      next
    }
    
    dataRes[[i]] = tibble(data.frame(t(rep(0, tibNVar))), .rows = 0)
    colnames(dataRes[[i]]) = tibColNames
    names(dataRes)[i] = allTibNames[i]
    
    # skip if data frame is an empty table
    isProperTibble = FALSE
    if(inherits(tib, "data.frame")){
      nRows = dim(tib)[1]
      isProperTibble = nRows > 0
    }
    
    if(!isProperTibble){
      next
    }
    
    # creates more rows
    dataRes[[i]][1:nRows, ] = 0
    
    tib_isNumericColumns = vapply(tib[1, ], FUN = is.numeric, FUN.VALUE = T, USE.NAMES = F)
    
    # set default values for character and boolean columns
    tib_isCharacterColumns = vapply(tib[1, ], FUN = is.character, FUN.VALUE = T, USE.NAMES = F)
    tib_isBooleanColumns = !tib_isNumericColumns & !tib_isCharacterColumns
    dataRes[[i]][, tib_isCharacterColumns] = ""
    dataRes[[i]][, tib_isBooleanColumns] = FALSE
    
    if(tibName == "yearlyResults"){
      
      identic_columns = c("year", "idFmCell", "longitude", "latitude", "Tmoy", "Tmax", "Tmin", "Rg", "PRI", "RH", "RU", "ca", "campbellCoef", "basePredawn", "wsat", "wilt",
                          "wfc", "outputDepth", "thinningNumber", "vetoFrost", "vetoDrought", "fit2018FileName", 
                          "incident_yearlyMJm2",
                          "incidentPAR_yearlyMJm2", "incidentPIR_yearlyMJm2",
                          "incidentDir_yearlyMJm2", "incidentDiff_yearlyMJm2",
                          "incidentPARdir_yearlyMJm2", "incidentPIRdir_yearlyMJm2",
                          "incidentPARdiff_yearlyMJm2", "incidentPIRdiff_yearlyMJm2")
      
      # dependant du nombre d'arbre de l'espèce
      BA_relative_columns= c("GPP", "altitude", "species", "Xposition", "Yposition", "LAImaxNextYear", "LAImaxBeforeFrost", "LAImaxThisYear", "LAIloss",
                             "NEE", "NPP", "Reco", "GPP", "RN", "ETR", "ETRveg", "ETRsol", "TR", "ETP", "CWD", "DBBV", "DBF", "DBRF", "Rauto", "Rmaintenance", "Rgrowth",
                             "Rwood", "RfineRoots", "Rcanopy", "Rhetero", "BSSmin", "BiomassOfLeaves", "BiomassOfReserves", "BiomassOfReservesBeforeRefill", "Ctop", "Csol", "ratioReserves",
                             "BiomassOfFineRoot", "BiomassOfTrunk", " BiomassOfCoarseRoot", "BiomassOfBranch", "AliveWoodBiomass", "BiomassNitrogen", "DBSS", "fineRootsMortality",
                             "drainage", "Psoilmin", "RDI", "Nha", "Gha", "Vha", "VhaInit", "reservetoBud", "numberOfSeeds", "numberOfFruits", "mortality", "thinnedVha", "NSC",
                             "PLCsoil", "rateOfMortality", "yearlyLeafGrowth", "yearlyLeafRespiration", "yearlyInputNitrogen", 
                             "RU_level_min", "REWmin", "RU_shortage_max",
                             "veg_yearlyMJm2", 
                             "vegPAR_yearlyMJm2", "vegPIR_yearlyMJm2", 
                             "vegDir_yearlyMJm2", "vegDiff_yearlyMJm2", 
                             "vegPARdir_yearlyMJm2", "vegPIRdir_yearlyMJm2", "vegPARdiff_yearlyMJm2", "vegPIRdiff_yearlyMJm2"
                             )
      
      # independant du nombre d'arbre de l'espèce (si on ajoute des arbres identiques sur le même espace)
      tree_relative_columns= c("age", "rw", "dbh", "dbhInitYear", "dbhInitSimulation", "BAI", "WVI", "WVI_cylinder", "WVI_real", "height", "leafSurfaceThisYear", "StressLevel", "istressBiljou", "gel", "kplantini",
                               "LeafDormancyBreakDate", "BBday", "maxFrostEffect", " FlowerDormancyBreakDate", "FloweringDate", "FruitGrowthInitDate", "FruitMaturationDate",
                               "endOfLeafAreaGrowth", "endOfLeafGrowth", "LeafSenescenceDate", " endLeaf", "beginFall", "dayOfWoodStop", "crownProjectionNextYear",
                               "crownProjectionThisYear", "Pleafmin", "TSUMBB", "Fcrit", "FHminfe", "g1", "nc", "nf", "LMA", "sureauSlope", "coefrac", "GBVmin", "woodStop", "CRBV",
                               "potsoiltowood", "treevolume", "tronviv", "treeNumberOfSeeds", "treeNumberOfFruits", "PLC", "Ksmax", 
                               "timberValue", "dayOptim")
      
      # dont know variables : Delta13C, rateOfSeeds, SeedMassMu, SeedMassSpecies, SeedMortality, ReservesToReproduce, AuxiliaryReproStructure, rateOfEmptySeeds, 
      # maximumWUE, maximumLeafRespirationRate, maximumWoodRespirationRate,  allocLeaves, DroughtEffectOnRespirationCoef
      # not concerned variables : 
      
      
      # For each column, the value is the weighted mean of the monospecific simulation's values
      for(j_res in 1:tibNVar){
        aVariable = tibColNames[j_res]
        
        if(tib_isNumericColumns[j_res]){
          if(aVariable %in% identic_columns){
            dataRes[[i]][ , j_res] = data1[[i]][ , j_res]
          }else if(aVariable %in% BA_relative_columns){
            dataRes[[i]][ , j_res] = proportionBA1 * data1[[i]][ , j_res] + (1 - proportionBA1) * data2[[i]][ , j_res]
          }else if(aVariable %in% tree_relative_columns){
            dataRes[[i]][ , j_res] = proportionN1 * data1[[i]][ , j_res] + (1 - proportionN1) * data2[[i]][ , j_res]
          }else{
            # do nothing
          }
        }else{ # For non numeric, values of data1 are used
          dataRes[[i]][ , j_res] = data1[[i]][ , j_res]
        }
      } # end loop on columns
    } # end yearlyResults
    
    
    
    if(tibName == "dailyResults" & outputDepth >= 2){
      
      identic_columns = c("idFmCell", "year", "day", "month", "Tmax", "Tmin", "Rg", "Rain")
      
      # dependant du nombre d'arbre de l'espèce
      BA_relative_columns= c("species", "GPP", "TR", "CanopyConductance", "ETRcan", "ETP", "LAIday", "dailyInputNitrogen", "DBBV", "DBRF", "BSS", "BF", "dailyRemainingCarbone", "MBF", "LAImax", "RVM", "RVC",
                             "maxHourlyTranspiration", "Rheterotrophic", "Ctop", "Csol", "LAIloss", "NEE", "Reco", "biomassOfFineRoot", "biomassOfReserves", "BiomassOfReservesBeforeRefill", "biomassOfTrunk", "biomassOfBranch",
                             "biomassOfCoarseRoot", "woodGrowth", "reservesGrowth", "fineRootsGrowth", "reservesMortality", "leafGrowth", "RCF", "woodMortality", "fineRootsMortality",
                             "coarseRootsMortality", "biomassOfNitrogen", "vegTot", "vegPAR", "vegPIR", "vegDir", "vegDif", "vegPARdir", "vegPARdif", "vegPIRdir", "vegPIRdif", "dailyPotentialPARdir",
                             "dailyAbsorbedPARdir", "lightCompetitionIndex", "inciTot", "inciPAR", "inciPIR", "inciDir", "inciDif", "inciPARdir", "inciPARdif", "inciPIRdir", "inciPIRdif", "soilTot",
                             "soilPAR", "soilPIR", "soilDir", "soilDif", "soilPARdir", "soilPARdif", "soilPIRdir", "soilPIRdif", "ETRsoil", "ETR", "reproGrowth", "RCRepro", "REW", "RU_currentLevel", "CWD",
                             "rsol", "rtop", "drainage", "drainageDeep", "rlit", "potsoil", "RU", "AsolPAR", "AsolPIR", "psd", "egt", "ec", "transpiration")
      
      # independant du nombre d'arbre de l'espèce (si on ajoute des arbres identiques sur le même espace)
      tree_relative_columns= c("lateFrost", "potlefmin", "TSUMBB", "SlopePotGs", "SlopeVcmax", "PLCmemo", "dayOptim", "crownProjection", "stomatalControl", "stressLevel", "istressBiljouDay",
                               "droughtIndexForFruitMaturation", "LMAmoy", "LMA5", "LMA4", "LMA3", "LMA2", "LMA1", "GBV", "GSS", "GRF", "GRG", "coefrac", "WoodStop", "ConcBSS", "GBVmemory", "RS", 
                               "PotentialNumberOfBuds", "PotentialNumberOfSeeds", "SeedNumber", "SeedMass")
      
      # dont know variables : Rnj d13C droughtEffectOnRespirationCoef frec reci QDIXAcclim 
      # not concerned variables : LAImaxAboveTreeBottomLine LAImaxAboveLayer1 LAImaxAboveLayer2 LAImaxAboveLayer3 LAImaxAboveLayer4 LAImaxAboveLayer5 
      # LAIaboveTreeBottomLine LAIaboveLayer1 LAIaboveLayer2 LAIaboveLayer3 LAIaboveLayer4 LAIaboveLayer5 species    
      
      
      # For each column, the value is the weighted mean of the monospecific simulation's values
      for(j_res in 1:tibNVar){
        aVariable = tibColNames[j_res]
        
        if(tib_isNumericColumns[j_res]){
          if(aVariable %in% identic_columns){
            dataRes[[i]][ , j_res] = data1[[i]][ , j_res]
          }else if(aVariable %in% BA_relative_columns){
            dataRes[[i]][ , j_res] = proportionBA1 * data1[[i]][ , j_res] + (1 - proportionBA1) * data2[[i]][ , j_res]
          }else if(aVariable %in% tree_relative_columns){
            dataRes[[i]][ , j_res] = proportionN1 * data1[[i]][ , j_res] + (1 - proportionN1) * data2[[i]][ , j_res]
          }else{
            # do nothing
          }
        }else{ # For non numeric, values of data1 are used
          dataRes[[i]][ , j_res] = data1[[i]][ , j_res]
        }
      } # end loop on columns
    } # end dailyResults
    
    
    
  } # end loop on simulation tables
  
  # 0. Copy paste tables that are identical between simulations
  
  for(tableName in c("soilFirstCell", "fmSettings", "PDGInitialParameters")){
    if(tableName %in% allTibNames){
      dataRes[[tableName]] = data1[[tableName]]
    }
  }
  
  # some meta informations
  dataRes$inventoryParent1 = data1$inventory
  dataRes$inventoryParent2 = data2$inventory
  dataRes$nameParent1 = data1$name
  dataRes$nameParent2 = data2$name
  
  dataRes$inventory = data1$inventory
  dataRes$inventory$originalinventory = paste0(data1$name, "__AND__", data2$name)
  dataRes$inventory$date = Sys.Date()
  dataRes$inventory$time = format(Sys.time(), "%H:%M:%OS")
  dataRes$inventory$dateTime = paste(dataRes$inventory$date, dataRes$inventory$time)
  
  if(!is.null(data1$name) & !is.null(data2$name)){
    dataRes$name = paste0(data1$name, "__COMBINEDWITH__", data1$name)
  }
  
  dataRes$combinationFromMonoSpMarker = paste0("This simulation is a combination from two monospecific simulation ", format(as.POSIXct(Sys.time()), format = "%d/%m/%Y, at %H:%M:%S"))
  return(dataRes)
} # end combineStandScaleSimulations




# Combine the tree of two irregular monospecific simulation
# Result is a irregular plurispecific simulation with no species interaction
# Trees that are conserved are those whose idFmCell < 1000000
combineIndividualScaleSimulation = function(data1, data2, name = NULL){
  
  # if data1 is a list of indivScaleSimu
  if(!"yearlyResults" %in% names(data1) & "yearlyResults" %in% names(data1[[1]])){
    return( combineIndividualScaleSimulationList(data1, data2, name) )
  }
  
  # if data are stand scale
  if(length(unique(data1$yearlyResults$idFmCell)) == 1){
    stop("data1 has only one individual.")
  }
  
  if(length(unique(data2$yearlyResults$idFmCell)) == 1){
    stop("data2 has only one individual.")
  }
  
  
  # 1. Create a list of tibble (= the simulation object) with rows of 0
  dataRes = list()
  
  allTibNames = names(data1)
  
  # Create new hollow tabs for non-empty tabs 
  for(i in 1:length(allTibNames)){
    tibName = allTibNames[i]
    tib = data1[[i]]
    tibColNames = colnames(tib)
    tibNVar = length(tibColNames)
    
    # Skip every tab that dont have idFmCell and some others
    if ( (! "idFmCell" %in% tibColNames) | (tibName %in% c("inventory", "fmSettings", "PDGInitialParameters", "soilFirstCell") )) {
      next
    }
    
    dataRes[[i]] = tibble(data.frame(t(rep(0, tibNVar))), .rows = 0)
    colnames(dataRes[[i]]) = tibColNames
    names(dataRes)[i] = allTibNames[i]
    
    # skip if data frame is an table
    isProperTibble = FALSE
    if(inherits(tib, "data.frame")){
      nRows = dim(tib)[1]
      isProperTibble = nRows > 0
    }
    
    if(!isProperTibble){
      next
    }
    
    # copy the tibble of data1
    dataRes[[i]] = data1[[i]]
    
    
    # separate original and species-converted trees
    legitimateIdFmcell1 = data1[[i]]$idFmCell[data1[[i]]$idFmCell < 1000000]
    legitimateIdFmcell2 = data2[[i]]$idFmCell[data2[[i]]$idFmCell < 1000000]
    
    illegitimateIdFmcell1 = data1[[i]]$idFmCell[data1[[i]]$idFmCell >= 1000000]
    illegitimateIdFmcell2 = data2[[i]]$idFmCell[data2[[i]]$idFmCell >= 1000000]
    
    # Check if legitimate trees of both table does not match
    if(sum(match(legitimateIdFmcell1, legitimateIdFmcell2, nomatch = 0)) != 0){
      stop(paste0("data1 and data2 have trees idFmCells in common for table ", tibName))
    }
    
    
    # Check if legitimate trees of data2 matches all with illegitimate trees of data1
    # double négation : le nombre d'arbres de legitimateIdFmcell2 + 100000 qui ne matche pas avec illegitimateIdFmcell1 est nul
    if(sum(match(legitimateIdFmcell2 + 1000000, illegitimateIdFmcell1, nomatch = 0) == 0)){
      stop(paste0("All illegitimate trees of data2 are not found in data 1 for table ", tibName))
    }
    
    
    # loop in legitimate trees of data2
    for(legitimateData2_idFmCell in unique(legitimateIdFmcell2)){
      dataRes[[i]][dataRes[[i]]$idFmCell == legitimateData2_idFmCell+1000000, ] = data2[[i]][data2[[i]]$idFmCell == legitimateData2_idFmCell, ]
    }
  
  } # end loop on simulation tables
  
  # 0. Copy paste tables that are identical between simulations
  for(tableName in c("soilFirstCell", "fmSettings", "PDGInitialParameters")){
    if(tableName %in% allTibNames){
      dataRes[[tableName]] = data1[[tableName]]
    }
  }
  
  # Some meta-information
  dataRes$inventoryParent1 = data1$inventory
  dataRes$inventoryParent2 = data2$inventory
  dataRes$nameParent1 = data1$name
  dataRes$nameParent2 = data2$name
  
  dataRes$inventory = data1$inventory
  dataRes$inventory$originalinventory = paste0(data1$name, "__AND__", data2$name)
  dataRes$inventory$date = Sys.Date()
  dataRes$inventory$time = format(Sys.time(), "%H:%M:%OS")
  dataRes$inventory$dateTime = paste(dataRes$inventory$date, dataRes$inventory$time)
  
  if(is.null(name)){
    dataRes$name = paste0(data1$name, "__COMBINEDWITH__", data1$name)
  }
  
  dataRes$combinationFromMonoSpMarker = paste0("This simulation is a combination from two monospecific simulation ", format(as.POSIXct(Sys.time()), format = "%d/%m/%Y, at %H:%M:%S"))
  return(dataRes)
} # combineIndividualScaleSimulation








# ..............................................................................
# STAND SCALE ------------------------------------------------------------------


# Compute stand variable based on individual variable
# for each variable : sum ( variable(tree) * crownProjection(tree) ) / stand area 
# variableSelection is a selection of variables in the unit "X / m2 of soil under crown projection"
computeStandCrownProjectionVariables = function(tableToFill, indivScaleTable, variableSelection, crownProjectionTable, standArea, 
                                                    years = NULL, days = NULL, hours = NULL, layers = NULL){
  
  # If given, years, days and hours should has the same number of rows than indivScaleTable. 
  # These elements indicate the temporality of each measures. 
  # The variables are summed in the scope of one unique temporality, eg. : Year 2003, Day 1 (if years and days are given) 
  # Year 2005, Day 15, Hour 23 (if years, days and hours are given)
  
  # filter Variables if not existing...
  variableSelection = variableSelection[variableSelection %in% names(indivScaleTable)]
  
  # sum elements times crown projection
  if(is.null(years)){ # no indication on years, days, hours..
    aggregateOnStandTab = t(apply( X = indivScaleTable[ , variableSelection] * crownProjectionTable, 
                                   MARGIN = 2, 
                                   FUN = sum) )
    
  }else if(is.null(days)){ # indication on years
    aggregateOnStandTab = aggregate(indivScaleTable[ , variableSelection] * crownProjectionTable, 
                                    by = list(year = years),
                                    FUN = sum)
    
    aggregateOnStandTab$element = aggregateOnStandTab$year
    tableToFill$element = tableToFill$year
    
  }else if(is.null(hours)){ # indication on years and days
    aggregateOnStandTab = aggregate(indivScaleTable[ , variableSelection] * crownProjectionTable, 
                                    by = list(year = years, day = days),
                                    FUN = sum)
    
    aggregateOnStandTab$element = paste0(aggregateOnStandTab$year, "_", aggregateOnStandTab$day)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day)
    
  }else if(is.null(layers)){ # indication on years, days and hours
    aggregateOnStandTab = aggregate(indivScaleTable[ , variableSelection] * crownProjectionTable, 
                                    by = list(year = years, day = days, hour = hours),
                                    FUN = sum) 
    
    aggregateOnStandTab$element = paste0(aggregateOnStandTab$year, "_", aggregateOnStandTab$day, "_", aggregateOnStandTab$hour)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day, "_", tableToFill$hour)
    
  }else{ # indication on years, days, hours and layers
    
    aggregateOnStandTab = aggregate(indivScaleTable[ , variableSelection] * crownProjectionTable, 
                                    by = list(year = years, day = days, hour = hours, k = layers),
                                    FUN = sum)
    
    aggregateOnStandTab$element = paste0(aggregateOnStandTab$year, "_", aggregateOnStandTab$day, "_", aggregateOnStandTab$hour, "_", aggregateOnStandTab$k)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day, "_", tableToFill$hour, "_", tableToFill$k)
  }
  
  # divide by total area to get back in X / m2
  aggregateOnStandTab[, variableSelection] = aggregateOnStandTab[, variableSelection] / standArea
  
  
  # replace values of the computed table into the original table
  if(!is.null(years)){
    tableToFill[, variableSelection] = aggregateOnStandTab[match(tableToFill$element, aggregateOnStandTab$element), variableSelection]
  }else{
    tableToFill[, variableSelection] = aggregateOnStandTab[ , variableSelection]
  }
  
  tableToFill$element = NULL
  
  return(tableToFill)
}




# Compute stand variable based on individual variable by averaging variables in meanVar
computeStandMeanVariables = function(tableToFill, indivScaleTable, variableSelection, years = NULL, days = NULL, hours = NULL, layers = NULL){
  # If given, years, days and hours should has the same number of rows than indivScaleTable. 
  # These elements indicate the temporality of each measures. 
  # The variables are averaged in the scope of one unique temporality, eg. : year 2003, day 1 (if years and days are given) 
  # Year 2005, Day 15, hours 23 (if years, days and hours are given)
  
  # filter Variables if not existing...
  variableSelection = variableSelection[variableSelection %in% names(indivScaleTable)]
  
  if(is.null(years)){ # no indication on years, days, hours..
    averagedTab = t( apply(indivScaleTable[, variableSelection], MARGIN = 2, FUN = mean) )
    
  }else if(is.null(days)){ # indication on years
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years), FUN = mean)
    
    averagedTab$element = averagedTab$year
    tableToFill$element = tableToFill$year
    
  }else if(is.null(hours)){ # indication on years and days
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years, day = days), FUN = mean)
    
    averagedTab$element = paste0(averagedTab$year, "_", averagedTab$day)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day)
    
  }else if(is.null(layers)){ # indication on years, days and hours
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years, day = days, hour = hours), FUN = mean)
    
    averagedTab$element = paste0(averagedTab$year, "_", averagedTab$day, "_", averagedTab$hour)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day, "_", tableToFill$hour)
    
  }else{ # indication on years, days, hours and layers
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years, day = days, hour = hours, k = layers), FUN = mean)
    
    averagedTab$element = paste0(averagedTab$year, "_", averagedTab$day, "_", averagedTab$hour, "_", averagedTab$k)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day, "_", tableToFill$hour, "_", tableToFill$k)
  }
  
  # replace values of the computed table into the original table
  if(!is.null(years)){
    tableToFill[, variableSelection] = averagedTab[match(tableToFill$element, averagedTab$element), variableSelection]
  }
  else{
    tableToFill[, variableSelection] = averagedTab[, variableSelection]
  }
  
  tableToFill$element = NULL
  
  return(tableToFill)
}





# Compute stand variable based by doing quadratic mean on individual variables
computeStandQuadraticMeanVariables = function(tableToFill, indivScaleTable,  variableSelection, years = NULL, days = NULL, hours = NULL, layers = NULL){
  
  quadraticMean = function(x){
    return( sqrt(mean(x**2)) )
  }
  
  # If given, years, days and hours should has the same number of rows than indivScaleTable. 
  # These elements indicate the temporality of each measures. 
  # The variables are averaged in the scope of one unique temporality, eg. : year 2003, day 1 (if years and days are given) 
  # Year 2005, Day 15, hours 23 (if years, days and hours are given)
  
  # filter Variables if not existing...
  variableSelection = variableSelection[variableSelection %in% names(indivScaleTable)]
  
  if(is.null(years)){ # no indication on years, days, hours..
    averagedTab = t( apply(indivScaleTable[, variableSelection], MARGIN = 2, FUN = quadraticMean) )
    
  }else if(is.null(days)){ # indication on years
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years), FUN = quadraticMean)
    
    averagedTab$element = averagedTab$year
    tableToFill$element = tableToFill$year
    
  }else if(is.null(hours)){ # indication on years and days
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years, day = days), FUN = quadraticMean)
    
    averagedTab$element = paste0(averagedTab$year, "_", averagedTab$day)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day)
    
  }else if(is.null(layers)){ # indication on years, days and hours
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years, day = days, hour = hours), FUN = quadraticMean)
    
    averagedTab$element = paste0(averagedTab$year, "_", averagedTab$day, "_", averagedTab$hour)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day, "_", tableToFill$hour)
    
  }else{ # indication on years, days, hours and layers
    averagedTab = aggregate(indivScaleTable[, variableSelection], by = list(year = years, day = days, hour = hours, k = layers), FUN = quadraticMean)
    
    averagedTab$element = paste0(averagedTab$year, "_", averagedTab$day, "_", averagedTab$hour, "_", averagedTab$k)
    tableToFill$element = paste0(tableToFill$year, "_", tableToFill$day, "_", tableToFill$hour, "_", tableToFill$k)
  }
  
  # replace values of the computed table into the original table
  if(!is.null(years)){
    tableToFill[, variableSelection] = averagedTab[match(tableToFill$element, averagedTab$element), variableSelection]
  }
  else{
    tableToFill[, variableSelection] = averagedTab[, variableSelection]
  }
  
  tableToFill$element = NULL
  
  return(tableToFill)
}






# Create a stand scale simulation object base on a PDG Arena simulation
# indivScaleSimu is the output of importSimu (list of tables). It can also be a list of simulation
# Simulation can be of one year or several year
# output can be reduced to make computation up to yearly (1), daily (2), hourly (3) or hourly x layer (4) level
# fmSpeciesIDs can be given to filter the tree belonging to a species or a set of species
getStandScaleSimu = function(indivScaleSimu, output = 4, fmSpeciesIDs = NULL){
  
  if(!is.null(indivScaleSimu$isEmpty)){
    indivScaleSimu$standScaleMarker = paste0("This tab was stand scaled the ", format(as.POSIXct(Sys.time()), format = "%d/%m/%Y, at %H:%M:%S"))
    return(indivScaleSimu)
  }
  
  if(!is.null(fmSpeciesIDs)){
    indivScaleSimu = sfiltering(indivScaleSimu, fmSpeciesIDs = fmSpeciesIDs)
  }
  
  # indivScaleSimu is a list of indivScaleSimu
  if(!"yearlyResults" %in% names(indivScaleSimu) & "yearlyResults" %in% names(indivScaleSimu[[1]])){
    return( getStandScaleSimuList(indivScaleSimu, output) )
  }
  
  # 1. Create a list of tibble (= the simulation object) with rows of 0
  standSimu = list()
  
  allTibNames = names(indivScaleSimu)
  
  standArea = indivScaleSimu$inventory$ncol * indivScaleSimu$inventory$nlin * indivScaleSimu$inventory$cellwidth**2 # in m2
  
  nFmCell = length(unique(indivScaleSimu$yearlyResults$idFmCell))
  oneFmCellId = unique(indivScaleSimu$yearlyResults$idFmCell)[1]
  years = unique(indivScaleSimu$yearlyResults$year)
  oneTreeSimu = sfiltering(indivScaleSimu, fmCellIDs = c(oneFmCellId))
  
  option = "pdgArenaDisableWaterInteraction"
  if("pdgLightTagMode" %in% names(indivScaleSimu$PDGInitialParameters)){
    # old denomination
    option = "pdgLightDisableWaterInteraction"
  }
  
  isWaterInteractionsDisabled = indivScaleSimu$PDGInitialParameters[[option]]
  
  
  # Create new hollow tabs for tree-level tabs 
  for(i in 1:length(allTibNames)){
    tibName = allTibNames[i]
    tib = indivScaleSimu[[i]]
    tibColNames = colnames(tib)
    tibNVar = length(tibColNames)
    
    
    if ("idFmCell" %in% tibColNames & ! tibName %in% c("standDailyResults", "soilWater", "inventory", "fmSettings", "PDGInitialParameters", "soilFirstCell")) {
      standSimu[[i]] = tibble(data.frame(t(rep(0, tibNVar))))
      colnames(standSimu[[i]]) = tibColNames
      names(standSimu)[i] = allTibNames[i]
      
      nRows = length(which(tib$idFmCell == oneFmCellId))
      
      if(nRows > 0)
        standSimu[[i]][1:nRows, ] = 0
    }
  }
  
  # 0. Copy paste all table that are identical between stand-scale and individual-scale simulation
  
  for(tableName in c("inventory", "soilFirstCell", "fmSettings", "PDGInitialParameters",
                     "standDailyResults", "soilWater", "fractionMaps", "standFractionMaps", "name")){
    if(tableName %in% allTibNames){
      standSimu[[tableName]] = indivScaleSimu[[tableName]]
    }
  }
  
  
  # 2. Compute some stand scale variables
  
  # yearly level
  if("yearlyResults" %in% allTibNames){
    
    
    # variable that are identitcal in all cells
    identicalVar = c("year", "longitude", "latitude", "altitude", "Tmoy", "Tmax", "Tmin", "Rg", 
                     "PRI", "RH", "ca", "climateFileName", "inventoryFileName", "outputFile", 
                     "incident_yearlyMJm2", "incidentPAR_yearlyMJm2", "incidentPIR_yearlyMJm2", 
                     "incidentDir_yearlyMJm2", "incidentDiff_yearlyMJm2", 
                     "incidentPARdir_yearlyMJm2", "incidentPIRdir_yearlyMJm2", 
                     "incidentPARdiff_yearlyMJm2", "incidentPIRdiff_yearlyMJm2")
    
    if(length(unique(indivScaleSimu$yearlyResults$species)) == 1){
      identicalVar = c(identicalVar, "species")
    }
    
    # variables to average linearly over trees to get stand level variables
    meanVar = c("age", "rw", "BAI", "WVI", "WVI_cylinder", "WVI_real", "height", "crownProjectionNextYear", "crownProjectionThisYear", "StressLevel",
                "TSUMBB", "Fcrit", "g1", "nc", "nf", "allocLeaves", "maximumWUE", "maximumLeafRespirationRate", "maximumWoodRespirationRate", 
                "PLC", "NSC", "Pleafmin", "Psoilmin", "campbellCoef", "basePredawn", "wsat", "wfc", "wilt", "RDI", "treevolume", 
                "LeafDormancyBreakDate", "BBday", "FlowerDormancyBreakDate", "FloweringDate", "FruitGrowthInitDate", "FruitMaturationDate", 
                "endOfLeafAreaGrowth", "endOfLeafGrowth", "LeafSenescenceDate", "endLeaf", "beginFall", "dayOfWoodStop",
                "LMA", "sureauSlope", "coefrac", "woodStop", "CRBV", "potsoiltowood")
    
    
    quadraticMeanVar = c("dbh", "dbhInitYear", "dbhInitSimulation")
    
    # variables to multiply by crown projection and then divide by stand area to get stand level variables
    crownProjectionRatioVar = c("GPP", "NPP", "Reco", "LAImaxThisYear", "LAImaxNextYear", "LAImaxBeforeFrost", "LAIloss", "NEE", "ETRveg", 
                                "TR", "ETP", "DBSS",
                                "DBBV", "DBF", "DBRF", "Rauto", "Rmaintenance", "Rgrowth", "Rwood", "RfineRoots", "Rcanopy", "Rhetero",
                                "BSSmin", "BiomassOfLeaves", "BiomassOfReserves", "BiomassOfReservesBeforeRefill", "BiomassOfFineRoot", "BiomassOfTrunk",
                                "BiomassOfCoarseRoot", "BiomassOfBranch", "fineRootsMortality", "AliveWoodBiomass", "BiomassNitrogen", "yearlyLeafGrowth", "yearlyLeafRespiration",
                                "veg_yearlyMJm2", "vegPAR_yearlyMJm2", "vegPIR_yearlyMJm2", "vegDir_yearlyMJm2", "vegDiff_yearlyMJm2", 
                                "vegPARdir_yearlyMJm2", "vegPIRdir_yearlyMJm2", "vegPARdiff_yearlyMJm2", "vegPIRdiff_yearlyMJm2")
    
    
    
    if(isWaterInteractionsDisabled){
      # warning : original crownProjection should be used, not the actual one
      variableSelection = c(crownProjectionRatioVar, "RU")
    }else{
      # Water variables are computed at stand scales and are then identical between trees
      identicalVar = c(identicalVar, "RU")
      meanVar = c(meanVar, "REWmin", "RU_level_min", "RU_shortage_max") # these variables should be identical, but we do an average just in case
    }
    
    # variables that are tree level and soil level : see instead the daily soilWater table generated at stand scale
    # c("drainage", "yearlyInputNitrogen", "ETR", "ETRsol")
    
    crownProjections = indivScaleSimu$yearlyResults$crownProjectionThisYear
    
    standSimu$yearlyResults[, identicalVar] = oneTreeSimu$yearlyResults[, identicalVar]
    
    standSimu$yearlyResults = computeStandMeanVariables(tableToFill = standSimu$yearlyResults, indivScaleTable = indivScaleSimu$yearlyResults,
                                                            variableSelection = meanVar,
                                                            years = indivScaleSimu$yearlyResults$year)
    
    standSimu$yearlyResults = computeStandQuadraticMeanVariables(tableToFill = standSimu$yearlyResults, indivScaleTable = indivScaleSimu$yearlyResults,
                                                            variableSelection = quadraticMeanVar,
                                                            years = indivScaleSimu$yearlyResults$year)
    
    
    
    standSimu$yearlyResults = computeStandCrownProjectionVariables(tableToFill = standSimu$yearlyResults, indivScaleTable = indivScaleSimu$yearlyResults,
                                                                       variableSelection = crownProjectionRatioVar,
                                                                       crownProjectionTable = crownProjections,
                                                                       standArea = standArea,
                                                                       years = indivScaleSimu$yearlyResults$year)
    
    
    # variable to recompute at stand scale
    # c("Nha", "Gha", "Vha", "VhaInit")
    for(ayear in unique(indivScaleSimu$yearlyResults$year)){
      indivScaleSimuYear = indivScaleSimu$yearlyResults[indivScaleSimu$yearlyResults$year == ayear, ]
      standArea_ha = standArea / 1e4
      standSimu$yearlyResults[standSimu$yearlyResults$year == ayear, ]$Nha = dim(indivScaleSimuYear)[1] / standArea_ha
      standSimu$yearlyResults[standSimu$yearlyResults$year == ayear, ]$Vha = sum(indivScaleSimuYear$treevolume) / standArea_ha
      standSimu$yearlyResults[standSimu$yearlyResults$year == ayear, ]$Gha = sum((indivScaleSimuYear$dbh/100/2)**2*pi) / standArea_ha
    }
    
    
  }
  
  # Daily level
  if(output >=2 & mean(c("dailyResults", "woodGrowth", "radiationDaily", "treesDailyResults") %in% allTibNames) == 1){
    
    # variable that are identitcal in all cells
    identicalVar = c("year", "day", "month", "Tmax", "Tmin", "Rg", "Rain")
    
    if(length(unique(indivScaleSimu$dailyResults$species)) == 1){
      identicalVar = c(identicalVar, "species")
    }
    
    # variables to average linearly over trees to get stand level variables
    meanVar = c("lateFrost", "potlefmin", "potsoil", "TSUMBB", "SlopePotGs",
                "SlopeVcmax", "crownProjection", "stomatalControl", "stressLevel", 
                "istressBiljouDay", "droughtIndexForFruitMaturation", "LMAmoy", "LMA1", "LMA2", "LMA3", "LMA4", "LMA5",
                
                # woodGrowth variables
                "GBV", "GSS", "GRF", "GRG", "coefrac", "WoodStop", "ConcBSS", 
                "RS", "PotentialNumberOfBuds", "PotentialNumberOfSeeds", "SeedNumber", "SeedMass",
                
                # treesDailyResults variables
                "lightCompetitionIndex",
                "LAImaxAboveTreeBottomLine", "LAImaxAboveLayer1", "LAImaxAboveLayer2", "LAImaxAboveLayer3", "LAImaxAboveLayer4", "LAImaxAboveLayer5",
                "LAIaboveTreeBottomLine", "LAIaboveLayer1", "LAIaboveLayer2", "LAIaboveLayer3", "LAIaboveLayer4", "LAIaboveLayer5")
    
    # variables to multiply by crown projection and then divide by stand area to get stand level variables
    crownProjectionRatioVar = c( "GPP", "TR", "CanopyConductance", "ETRcan", "ETP", "LAIday",
                                 "dailyInputNitrogen", "DBBV", "BSS", "MBF", "LAImax", 
                                 "RVM", "RVC", "maxHourlyTranspiration", "LAIloss", "NEE", "Rheterotrophic",
                                 
                                 # woodGrowth variables
                                 "BF", "biomassOfFineRoot", "biomassOfReserves", "BiomassOfReservesBeforeRefill", "biomassOfTrunk",
                                 "biomassOfBranch", "biomassOfCoarseRoot", "woodGrowth",
                                 "reservesGrowth", "fineRootsGrowth", "reservesMortality",
                                 "leafGrowth", "RCF", "woodMortality",
                                 "fineRootsMortality", "coarseRootsMortality", "reproGrowth",
                                 "RCRepro", "biomassOfNitrogen", "dailyRemainingCarbone",
                                 
                                 #radiationDaily variables
                                 "vegTot", "vegPAR", "vegPIR", "vegDir", "vegDif", 
                                 "vegPARdir", "vegPARdif", "vegPIRdir", "vegPIRdif")
    
    # Dont know variable : "d13C"  "droughtEffectOnRespirationCoef" "frec"  "reci" "QDIXAcclim"   "dayOptim"  "Reco"  "ETR", 
    
    # 16.01.2023
    # variables "REW", "rsol", "RU_currentLevel", "drainage" are now found in standDailyResults
    
    standSimu$dailyResults[, identicalVar] = oneTreeSimu$dailyResults[, identicalVar]
    
    # new way
    standSimu$dailyResults = computeStandMeanVariables(tableToFill = standSimu$dailyResults,
                                                           indivScaleTable = indivScaleSimu$dailyResults, 
                                                           variableSelection = meanVar, 
                                                           years = indivScaleSimu$dailyResults$year, 
                                                           days = indivScaleSimu$dailyResults$day)
    
    
    
    standSimu$dailyResults = computeStandCrownProjectionVariables(tableToFill = standSimu$dailyResults,
                                                                      indivScaleTable = indivScaleSimu$dailyResults,
                                                                      crownProjectionTable = indivScaleSimu$dailyResults$crownProjection,
                                                                      variableSelection = crownProjectionRatioVar,
                                                                      standArea = standArea,
                                                                      years = indivScaleSimu$dailyResults$year,
                                                                      days = indivScaleSimu$dailyResults$day)
    
    
    # re-computation of REW
    if(isWaterInteractionsDisabled){
      
      yearOne = min(indivScaleSimu$yearlyResults$year)
      indivScaleSimuYearlyResultsYearOne = sfiltering(indivScaleSimu$yearlyResults, years = yearOne)
      
      standRU = sum(indivScaleSimuYearlyResultsYearOne$RU * indivScaleSimuYearlyResultsYearOne$crownProjectionThisYear) / standArea
      standSimu$dailyResults$REW = standSimu$dailyResults$RU_currentLevel / standRU
    }
    
    
    
  } # end daily level variables
  
  
  
  # Daily level
  if(output >=2 & ("dailyMeteo" %in% allTibNames)){
    
    # variable that are identitcal in all cells
    identicalVar = colnames(standSimu$dailyMeteo)
    
    standSimu$dailyMeteo[, identicalVar] = oneTreeSimu$dailyMeteo[, identicalVar]
    
  } # end daily level variables
  
  
  
  
  
  # add stand daily results into dailyResults
  if("standDailyResults" %in% allTibNames){
    # standDailyResults overwrite common variables (eg, tree-level soil variables)
    duplicateColNames = colnames(standSimu$dailyResults) %in% colnames(standSimu$standDailyResults) & (! colnames(standSimu$dailyResults) %in% c("year", "day"))
    standSimu$dailyResults[duplicateColNames] = NULL
    
    standSimu = joinTable(standSimu, "dailyResults", "standDailyResults")
  }
  
  
  
  
  # # add soilWater into dailyResults, only if water competition if held at stand scale
  # # if water interactions are disable, not needed
  # if(!isWaterInteractionsDisabled){
  #   if("soilWater" %in% allTibNames){
  #     # soilWater overwrite common variables (eg, tree-level soil variables)
  #     duplicateColNames = colnames(standSimu$dailyResults) %in% colnames(standSimu$soilWater) & (! colnames(standSimu$dailyResults) %in% c("year", "day"))
  #     standSimu$dailyResults[duplicateColNames] = NULL
  #     
  #     standSimu = joinTable(standSimu, "dailyResults", "soilWater")
  #   }
  # }
  
  
  # hourly level
  if(output >= 3 & "hourlyResults" %in% allTibNames){
    # variable that are identitcal in all cells
    identicalVar = c("year", "day", "hour", "PARh", "PARdiro", "PARdifo", "Sgh", "skyl", "th", "RHh", "Soh", "frach", "beta", "lat", "long")
    
    if(length(unique(indivScaleSimu$hourlyResults$species)) == 1){
      identicalVar = c(identicalVar, "species")
    }
    
    # variables to average linearly over trees to get stand level variables
    meanVar = c("potsoil", "Tveg", "crownProjection")
    
    # variables to multiply by crown projection and then divide by stand area to get stand level variables
    crownProjectionRatioVar = c("canopyPhotosynthesis", "canopyRespiration", "transpiration", "evapotranspiration", "canopyConductance", "LAItot", "strat")
    
    # Dont know variable : rbl rnhsol;rnhveg;ras;hourlySoilEvaporation ; canopyDelta13C ; NetRadiationFromLayer ; gcsun;coefSun
    
    
    standSimu$hourlyResults[, identicalVar] = oneTreeSimu$hourlyResults[, identicalVar]
    
    # new way
    # print("V2")
    standSimu$hourlyResults = computeStandMeanVariables(tableToFill = standSimu$hourlyResults,
                                                            indivScaleTable = indivScaleSimu$hourlyResults,
                                                            variableSelection = meanVar,
                                                            years = indivScaleSimu$hourlyResults$year, 
                                                            days  = indivScaleSimu$hourlyResults$day, 
                                                            hours = indivScaleSimu$hourlyResults$hour)
    
    standSimu$hourlyResults = computeStandCrownProjectionVariables(tableToFill = standSimu$hourlyResults,
                                                                       indivScaleTable = indivScaleSimu$hourlyResults,
                                                                       crownProjectionTable = indivScaleSimu$hourlyResults$crownProjection,
                                                                       variableSelection = crownProjectionRatioVar,
                                                                       standArea = standArea,
                                                                       years = indivScaleSimu$hourlyResults$year,
                                                                       days = indivScaleSimu$hourlyResults$day,
                                                                       hours = indivScaleSimu$hourlyResults$hour)
    
    
    
  } # end hourlyresults table
  
  
  if("radiationAbsorbedCanopy" %in% allTibNames){
    
  }
  
  if("radiationAbsorbedLayer" %in% allTibNames){
    
  }
  
  if("inputRadiation" %in% allTibNames){
    
  }
  
  if("radiationSoil" %in% allTibNames){
    
  }
  
  if("thermic" %in% allTibNames){
    
  }
  
  # hourly x layer level
  if(output >= 4 & "photosynthesis" %in% allTibNames){
    # variable that are identitcal in all cells
    identicalVar = c("year", "day", "hour", "k", "RH", "ws", "sureau")
    
    # variables to average linearly over trees to get stand "level" variables
    meanVar = c("propsun", "propshade", "potsoil",  "NLAI", "LMAlayer", "crownProjection", "PARdif", "PARsun","vcmax", "vjmax", "tleaf", "nc", "g1", "QDIX", "netRadiation",
                "grossPhotosynthesis", "waterConductance", "respiration")
    
    # variables to multiply by crown projection and then divide by stand area to get stand level variables
    crownProjectionRatioVar = c("strat", "aboveLAI")
    
    
    
    # Dont know variable : rbl
    # questioning : NLAI, vcmax, vjmax, tleaf, nc, g1, QDIX, netRadiation
    
    standSimu$photosynthesis[, identicalVar] = oneTreeSimu$photosynthesis[, identicalVar]
    
    # new way
    # print("V2")
    standSimu$photosynthesis = computeStandMeanVariables(tableToFill = standSimu$photosynthesis,
                                                             indivScaleTable = indivScaleSimu$photosynthesis,
                                                             variableSelection = meanVar,
                                                             years = indivScaleSimu$photosynthesis$year,
                                                             days  = indivScaleSimu$photosynthesis$day,
                                                             hours = indivScaleSimu$photosynthesis$hour,
                                                             layers =indivScaleSimu$photosynthesis$k)
    
    standSimu$photosynthesis = computeStandCrownProjectionVariables(tableToFill = standSimu$photosynthesis,
                                                                        indivScaleTable = indivScaleSimu$photosynthesis,
                                                                        crownProjectionTable = indivScaleSimu$photosynthesis$crownProjection,
                                                                        variableSelection = crownProjectionRatioVar,
                                                                        standArea = standArea,
                                                                        years = indivScaleSimu$photosynthesis$year,
                                                                        days  = indivScaleSimu$photosynthesis$day,
                                                                        hours = indivScaleSimu$photosynthesis$hour,
                                                                        layers = indivScaleSimu$photosynthesis$k)
    
    
    
  } # end photosynthesis table
  
  standSimu$standScaleMarker = paste0("This tab was stand scaled the ", format(as.POSIXct(Sys.time()), format = "%d/%m/%Y, at %H:%M:%S"))
  return(standSimu)
}



# take a list of simulation or one simulation
# returns a table with information about the simulation
simuExamination = function(simuList, standSimuList = NULL){
  
  # if one simu, transform the simu into a list of simu
  if("yearlyResults" %in% names(simuList)){
    simuList = list(simu = simuList)
    
    
    if(!is.null(standSimuList)){
      standSimuList = list(simu = standSimuList)
    }
  }
  
  # progress bar
  n_iter <- length(names(simuList)) # Number of iterations of the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  table = NULL
  for(simulationName in names(simuList)){
    
    oneSimu = simuList[[simulationName]]
    
    if(!is.null(standSimuList)){
      oneSimu_st = standSimuList[[simulationName]]
    }else{
      oneSimu_st = getStandScaleSimu(oneSimu)
    }
    
    
    
    nYears = length(unique(oneSimu$yearlyResults$year))
    
    
    firstYear = min(oneSimu$yearlyResults$year)
    firstYearOneSimu = sfiltering(oneSimu, years = firstYear)
    firstYearOneSimu_st = sfiltering(oneSimu_st, years = firstYear)
    
    
    lastYear = max(oneSimu$yearlyResults$year)
    lastYearOneSimu = sfiltering(oneSimu, years = lastYear)
    lastYearOneSimu_st = sfiltering(oneSimu_st, years = lastYear)
    
    lastDayLastYear = max(lastYearOneSimu$dailyResults$day)
    lastDayLastYearOneSimu = lastYearOneSimu$dailyResults[lastYearOneSimu$dailyResults$day == lastDayLastYear, ]
    lastDayLastYearOneSimu_st = lastYearOneSimu_st$dailyResults[lastYearOneSimu_st$dailyResults$day == lastDayLastYear, ]
    
    perNaYearly = mean(is.na(oneSimu$yearlyResults))
    perNaDaily = mean(is.na(oneSimu$dailyResults))
    perNegBSS = mean(lastDayLastYearOneSimu$BSS < 0)
    
    meanDbh = mean(lastYearOneSimu$yearlyResults$dbh)
    sdDbh = sd(lastYearOneSimu$yearlyResults$dbh)
    
    whichNegBss = lastDayLastYearOneSimu[lastDayLastYearOneSimu$BSS < 0, ]$idFmCell
    meanDbhNegBSS = NA
    sdDbhNegBSS = NA
    
    perBeechNegBSS = NA
    if(length(whichNegBss) > 0){
      meanDbhNegBSS = mean(lastYearOneSimu$yearlyResults[whichNegBss,]$dbh)
      sdDbhNegBSS = sd(lastYearOneSimu$yearlyResults[whichNegBss,]$dbh)
      
      tab = table(lastYearOneSimu$yearlyResults[whichNegBss,]$speciesName)
      perBeechNegBSS = 0
      if("Fagus sylvatica" %in% names(tab)){
        perBeechNegBSS = tab[["Fagus sylvatica"]] / sum(tab)
      }
    }
    
    table_line = tibble( simulationName = simulationName, 
                         simulationName2 = oneSimu$name,
                         meanDbh = meanDbh,
                         sdDbh = sdDbh,
                         nYears = nYears,
                         lastYear = lastYear,
                         perNaYearly = mean(is.na(oneSimu$yearlyResults)),
                         perNaDaily = mean(is.na(oneSimu$dailyResults)),
                         perNegBSS = perNegBSS,
                         perBeechNegBSS = perBeechNegBSS,
                         meanDbhNegBSS = meanDbhNegBSS,
                         sdDbhNegBSS = sdDbhNegBSS,
                         lastDayTreesMeanBSS = mean(lastDayLastYearOneSimu$BSS),
                         lastDayTreeMinBSS = min(lastDayLastYearOneSimu$BSS),
                         lastDayStandBSS = lastDayLastYearOneSimu_st$BSS,
                         lastYearStandMaxBSS = max(lastYearOneSimu_st$dailyResults$BSS),
                         firstYearLAImax = firstYearOneSimu_st$yearlyResults$LAImaxThisYear,
                         standLAImaxInventory = oneSimu$inventory$standlai)
    
    if(simulationName == names(simuList)[1]){
      table = table_line
    }else{
      table = rbind(table, table_line)
    }
    
    setTxtProgressBar(pb, which(names(simuList) == simulationName)) # change progress bar
  }
  close(pb) # Close progress bar
  return(table)
} # end method simuExamination


# Find the prefix of log files in a given folder
findLogPrefixList = function(folderPath){
  fileList = list.files(folderPath)
  filteredFileList = fileList[grepl(fileList, pattern = "yearlyResults")]
  logPrefixList = sapply(filteredFileList, FUN = function(x) strsplit(x, split = "yearlyResults.log")[[1]])
  names(logPrefixList) = NULL
  logPrefixList = logPrefixList[!grepl(logPrefixList, pattern = "prelimCastSim")]
  return(logPrefixList)
}






# ..............................................................................
# MULTIPLE SIMULATIONS ---------------------------------------------------------


# Import all simulations found in logPrefixList.log (created in var)
# simuSetPrefix : string to find the logPrefixList_simuSetPrefix.txt which contain all log prefixes
# skipFilter : select simulation whose prefix does not correspond with skipFilter
# keepFilter : select simulation whose prefix correspond to keepFilter
# preliminarySimulation : TRUE to import corresponding preliminary simulations
importSimuList = function(folderPath, years = NULL, simuSetPrefix = "", skipFilter = NULL, keepFilter = NULL, logList = NULL, preliminarySimulation = FALSE){
  
  # find logPrefix list
  logPrefixList = findLogPrefixList(folderPath);
  
  # preliminary CASTANEA simulation
  if(preliminarySimulation){
    logPrefixList = paste0("prelimCastSim_", logPrefixList)
  }
  
  # remove logPrefix that do not contain keepFilter
  if (!is.null(keepFilter)) {
    keepSimu = grepl(logPrefixList, pattern = keepFilter)
    logPrefixList = logPrefixList[keepSimu]
  }
  
  # remove logPrefix that do contain skipFilter
  if (!is.null(skipFilter)) {
    skip = grepl(logPrefixList, pattern = skipFilter)
    logPrefixList = logPrefixList[!skip]
  }
  
  if(length(logPrefixList) == 0){
    stop(paste0("No simulation has been found at ", folderPath))
  }
  
  cat(paste0("Importing ", length(logPrefixList), " simulations...\n"))
  cat("list : \t", logPrefixList)
  cat("\n")
  
  
  cat("from :")
  cat("\t", folderPath)
  cat("\n")
  
  if(is.null(logList)){
    logList = logListList[[length(logListList)]]
  }
  
  # remove last underscore
  logPrefixList2 = sapply(logPrefixList, FUN = function(x) substr(x, 1, nchar(x)-1), simplify = TRUE, USE.NAMES = FALSE)
  
  simuList = list()
  
  # progress bar
  n_iter <- length(logPrefixList) # Number of iterations of the progress bar
  # pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
  # max = n_iter, # Maximum value of the progress bar
  # style = 3,    # Progress bar style (also available style = 1 and style = 2)
  # width = 50,   # Progress bar width. Defaults to getOption("width")
  # char = "=")   # Character used to create the bar
  
  
  if(PROGRAM_ON_SERVER){
    pb_serv <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                         max = 100, # Maximum value of the progress bar
                         style = 3,    # Progress bar style (also available style = 1 and style = 2)
                         width = 50,   # Progress bar width. Defaults to getOption("width")
                         char = "=")   # Character used to create the bar
  }else{
    # linux only ! Windows : use winProgressBar (or comment this line)
    pb <- tkProgressBar(min = 0, max = 100, initial = 0, width = 500)
  }
  
  
  for(i in 1:length(logPrefixList)){
    logPrefix = logPrefixList[i]
    logPrefix2 = logPrefixList2[i] # without last underscore
    
    progressBarValue = round(100*(i-1)/n_iter)
    
    if(PROGRAM_ON_SERVER){
      setTxtProgressBar(pb_serv, progressBarValue) # change progress bar
    }else{
      # linux only ! Windows : use winProgressBar (or comment this line)
      setTkProgressBar(pb, value = progressBarValue, title = "Importing simulations from a list", label = paste0(progressBarValue, " %", " - Importing ", logPrefix2))
    }
    
   
    if(!file.exists(folderPath)){
      if(PROGRAM_ON_SERVER){
        close(pb_serv)
      }else{
        close(pb)
      }
      stop(paste0("Folder ", folderPath, " does not exists."))
    }
    simuList[[logPrefix2]] = importSimu(folderPath, prefix = logPrefix, logList = logList, printMessage = TRUE, shortMessage = TRUE, years = years)
    
    if(preliminarySimulation){
      simuList[[logPrefix2]]$standScaleMarker = "This is a CASTANEA preliminary simulation."
    }
    
    
    if(simuSetPrefix != ""){
      simuList[[logPrefix2]]$name = sub(x = simuList[[logPrefix2]]$name, pattern = simuSetPrefix, replacement = "")
    }
    
  }
  
  if(PROGRAM_ON_SERVER){
    close(pb_serv) # Close progress bar
  }else{
    close(pb) # Close progress bar
  }
  
  cat(paste0("\n", length(logPrefixList), " simulations imported :\n"))
  cat(paste(logPrefixList, sep = "    "))
  cat("\n")
  
  cat("\nDone\n")
  
  if(length(simuList) == 1){
    cat("\nOnly one simulation has been found. The simulation is return directly (not inside a list).\n")
    simuList = simuList[[1]]
  }
  
  return(simuList)
} # end importSimuList


# get stand scale simulation for all simulation of simuList
getStandScaleSimuList = function(simuList, output = 4){
  
  # progress bar
  n_iter <- length(names(simuList)) # Number of iterations of the progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  
  standSimuList = list()
  
  # for each simulation of the list
  for(simulationName in names(simuList)){
    standSimuList[[simulationName]] = getStandScaleSimu(simuList[[simulationName]])
    setTxtProgressBar(pb, which(names(simuList) == simulationName)) # change progress bar
  }
  
  close(pb) # Close progress bar
  
  return(standSimuList)
}



# specific to GMAP
# add composition and triplet information in yearlyResults and dailyResults
addCompositionInfo = function(simuList){
  for(site in names(simuList)){
    
    
    for(composition_prop in c("_m", "_sp", "_ph")){
      if(grepl(site, pattern = composition_prop) ){
        composition_code = composition_prop
      }
    }
    
    triplet = paste0(strsplit(site, composition_code)[[1]], collapse = "")
    
    composition = substr(composition_code, 2, nchar(composition_prop))
    
    simuList[[site]]$yearlyResults$triplet = triplet
    simuList[[site]]$yearlyResults$composition = composition
    
    
    simuList[[site]]$dailyResults$triplet = triplet
    simuList[[site]]$dailyResults$composition = composition
    
  }
  
  return(simuList)
}


# For a list of simulation, combine the value of two regular monospecific (stand scale) simulation
# Result is a regular plurispecific (stand scale) simulation with no species interaction
# --
# @data1 and 2 should be lists of regular monospecific stand scale simulation, with one simulation per inventory and representing different species of the same plots
# @proportionBA1 should be a vector containing the initial proportion of species 1 in the mixture for every inventory
# @dataIrregPluriSp (possibly provided instead of proportionBA1) should be a list of irregular plurispecific simulation, with one simulation per inventory  
# Order of every input variable should match !
combineStandScaleSimulationsList = function(data1, data2, 
                                           proportionBA1 = NULL, proportionN1 = NULL, dataIrregPluriSp = NULL, outputDepth = 1){
  
  # check data
  if(!inherits(data1, 'list')){
    stop("data1 is not a list. It should be a list of simulation.")
  }else if(!inherits(data1[[1]], 'list') ){
    stop("data1[[1]] is not a list. It should be a simulation's list.")
  }else if(!"yearlyResults" %in% names(data1[[1]]) ){
    stop("data1[[1]] does not contains yearlyResults")
  }
  # check data
  if(!inherits(data2, 'list')){
    stop("data2 is not a list. It should be a list of simulation.")
  }else if(!inherits(data2[[1]], 'list') ){
    stop("data2[[1]] is not a list. It should be a simulation's list.")
  }else if(!"yearlyResults" %in% names(data2[[1]]) ){
    stop("data2[[1]] does not contains yearlyResults")
  }
  # check data
  if(!is.null(dataIrregPluriSp)){
    if(!inherits(dataIrregPluriSp, 'list')){
      stop("dataIrregPluriSp is not a list. It should be a list of simulation.")
    }else if(!inherits(dataIrregPluriSp[[1]], 'list') ){
      stop("dataIrregPluriSp[[1]] is not a list. It should be a simulation's list.")
    }else if(!"yearlyResults" %in% names(dataIrregPluriSp[[1]]) ){
      stop("dataIrregPluriSp[[1]] does not contains yearlyResults")
    }
  }
  # check data
  if((!is.null(proportionBA1) | !is.null(proportionN1)) & !is.null(dataIrregPluriSp)){
    stop("proportionBA1 and proportionN1 OR dataIrregPluriSp should be given, not both")
  } 
  
  nInventories = length(data1)
  
  # check length of list ..
  if(nInventories != length(data2)){
    stop("The number of simulation is different between data1 and data2")
  }
  # check length of list ..
  if(!is.null(proportionN1)){
    if(nInventories != length(proportionN1)){
      stop("The length of data1 and proportionN1 is different")
    }
  }
  # check length of list ..
  if(!is.null(proportionBA1)){
    if(nInventories != length(proportionBA1)){
      stop("The length of data1 and proportionBA1 is different")
    }
  }
  
  # check length of list ..
  if(!is.null(dataIrregPluriSp)){
    if(nInventories != length(dataIrregPluriSp)){
      stop("The number of simulation is different between data1 and dataIrregPluriSp")
    }
  }
  
  
  resList = list()
  
  for(i in 1:nInventories){
    aName = getCodeSiteFromSimulationName(names(data1)[[i]])
    resList[[aName]] = combineStandScaleSimulations(data1[[i]], 
                                                  data2[[i]],
                                                  proportionBA1[[i]],
                                                  proportionN1[[i]],
                                                  dataIrregPluriSp[[i]],
                                                  outputDepth)
  }
  
  return(resList)
} # end combineStandScaleSimulationsList



# For a list of simulation, combine the value of two irregular monospecific simulation
# Result is a irregular plurispecific simulation with no species interaction
# --
# @data1 and 2 should be lists of irregular monospecific simulation, with one simulation per inventory and representing different species of the same plots
# @name is optionnal vector to give a name to each result
# Order of every input variable should match !
combineIndividualScaleSimulationList = function(data1, data2, name = NULL){
  
  # check data
  if(!inherits(data1, 'list')){
    stop("data1 is not a list. It should be a list of simulation.")
  }else if(!inherits(data1[[1]], 'list') ){
    stop("data1[[1]] is not a list. It should be a simulation's list.")
  }else if(!"yearlyResults" %in% names(data1[[1]]) ){
    stop("data1[[1]] does not contains yearlyResults")
  }
  # check data
  if(!inherits(data2, 'list')){
    stop("data2 is not a list. It should be a list of simulation.")
  }else if(!inherits(data2[[1]], 'list') ){
    stop("data2[[1]] is not a list. It should be a simulation's list.")
  }else if(!"yearlyResults" %in% names(data2[[1]]) ){
    stop("data2[[1]] does not contains yearlyResults")
  }
  
  nInventories = length(data1)
  
  # check length of list ..
  if(nInventories != length(data2)){
    stop("The number of simulation is different between data1 and data2")
  }
  # check length of list ..
  if(!is.null(name)){
    if(nInventories != length(name)){
      stop("The length of data1 and name is different")
    }
  }
  
  
  resList = list()
  
  for(i in 1:nInventories){
    name = getCodeSiteFromSimulationName(names(data1)[i])
    resList[[name]] = combineIndividualScaleSimulation(data1[[i]], 
                                                       data2[[i]],
                                                  name[i])
  }
  
  return(resList)
} # end combineIndividualScaleSimulationList 


# return a result in form of a list of getInitialSpeciesBasalAreaProportionInSimulation
# for every simulation of simuList
getInitialSpeciesBasalAreaProportionInSimulationList = function(simuList, targetSp){
  res = list()
  
  for(aname in names(simuList)){
    res[[aname]] = getInitialSpeciesBasalAreaProportionInSimulation(simuList[[aname]], targetSp)
  }
  return(res)
}

# return a result in form of a list of getInitialSpeciesTreeProportionInSimulation
# for every simulation of simuList
getInitialSpeciesTreeProportionInSimulationList = function(simuList, targetSp){
  res = list()
  
  for(aname in names(simuList)){
    res[[aname]] = getInitialSpeciesTreeProportionInSimulation(simuList[[aname]], targetSp)
  }
  return(res)
}


cleanSimuListNames = function(simuList){
  
  for(i in 1:length(names(simuList))){
    # replace all brackets and dash by nothing
    names(simuList)[i] = gsub(x = names(simuList)[i], pattern = "(\\[|\\]|\\)|\\(|-)", replacement = "")
  }
  
  return(simuList)
}


# AREAS COMPARISON -------------------------------------------------------------

# SEVERAL GRAPHS IN ONE PLOT, 
# Area representation,perfect for following daily or hourly variation along year / day
# Can show multiple variable if yVar is a list
#
# (prefered) if tableName is provided, data1 and data2 must be a list of tables
# (old entry) if not, data1 and data2 are datatables (tibbles..)
# dataNames is a list of the names of data1 and data2, used in the legend
# xvar is a string, eg. "day" (must be a column name of data1 and data2)
# yvar is a string or a list of string, eg. c('GPP','REW','LAIday')  (must be column names of data1 and data2)
# Facultative parameter :
# coord_cartesian to zoom on a part of the graph, eg. coord_cartesian(xlim = c(100,200), ylim = c(5,10)) 
# type = "cloud" to represent data with points
# alphaValue to adjust the transparency of the areas or points, from 0 to 1
simuPlots = function(data1_, data2_ = NULL, tableName = NULL, dataNames = NULL, xvar, yvar, 
                     title = "", coord_cartesian = NULL, type = "area", alphaValue = 0.5, 
                     years = NULL, fmSpeciesIDs = NULL, fmCellIDs = NULL, layers = NULL,
                     nColMultiPlot = NULL, dataReductionFactor = NULL,
                     treeLegend = FALSE, sortYvar = FALSE){
  
  # EDIT title
  if(title == "" & !is.null(suppressWarnings(data1_$name))){
    title = data1_$name
    if(!is.null(data2_$name)){
      if(data2_$name != data1_$name){
        title = paste0(title, " ", data2_$name)
      }
    }
  }
  
  annotation = ""
  
  if(!is.list(data1_) & !is.data.frame(data1_)){
    stop("Error : data1_ is not a list neither a data frame ")
  }
  
  # if tableName is not provided and data1_ is not a table, problem
  if(is.null(tableName) & !is.data.frame(data1_)){
    stop("Error : data1_ is not a data frame and tableName was not provided ")
  }
  
  if(!is.null(tableName) & is.data.frame(data1_)){
    warning("Warning : tableName was provided but data1_ is already a table. tableName will be ignored.")
  }
  
  isStandScale = FALSE
  
  if(!is.null(tableName)){
    
    # Identify table name if a diminuted name of provided
    tableName = identifyTableName(tableName, names(data1_))
    
    
    if(tableName %in% names(data1_)){
      data1 = data1_[[tableName]]
    }else{
      stop(paste0("tableName ", tableName, " was not recognized"))
    }
    
    isData1StandScale = !is.null(data1_$standScaleMarker)
    isStandScale = isData1StandScale
    
    if(!is.null(data2_)){
      data2 = data2_[[tableName]]
      isData2StandScale = !is.null(data2_$standScaleMarker)
      
      if(isData1StandScale != isData2StandScale){
        warning("Warning: data1_ and data2_ represent respectively stand scale and individual scale simulations.")
      }
    }else{
      data2 = data2_
    }
    
  }else{
    data1 = data1_
    data2 = data2_
  }
  
  # check if multiple years
  # if this is a daily simulation but no year was provided
  if("year" %in% names(data1) & "day" %in% names(data1) & is.null(years)){
    if(length(unique(data1$year)) > 1){
      years = unique(data1$year)[1]
      cat(paste0("\nWarning: several years were provided in the table. Year ", years, " was chosen."))
    }
  }
  
  
  if(!is.null(years)){
    data1 = sfiltering(data1, years = years)
    data2 = sfiltering(data2, years = years)
  }
  
  if(!is.null(fmSpeciesIDs)){
    data1 = sfiltering(data1, fmSpeciesIDs = fmSpeciesIDs)
    data2 = sfiltering(data2, fmSpeciesIDs = fmSpeciesIDs)
  }
  if(!is.null(fmCellIDs)){
    data1 = sfiltering(data1, fmCellIDs = fmCellIDs)
    data2 = sfiltering(data2, fmCellIDs = fmCellIDs)
  }
  if(!is.null(layers)){
    data1 = sfiltering(data1, layers = layers)
    data2 = sfiltering(data2, layers = layers)
  }
  
  
  # verification of nullity of data tables
  if(dim(data1)[1] <= 0){
    stop("Error: data1 is empty !")
  }
  
  if(!is.null(data2)){
    if(dim(data2)[1] <= 0){
      stop("Error: data2 is empty !")
    }
  }
  
  # verification of colnames presence
  if(mean(yvar %in% colnames(data1)) != 1){
    varNotInData = which(!yvar %in% colnames(data1))
    varNotInDataString = paste0(yvar[varNotInData], collapse = ", ")
    warning(paste0( "Warning: these yvar are not in colnames of data1 and are ignored:\n", varNotInDataString))
    
    yvar = yvar[-varNotInData]
  }
  
  if(!is.null(data2)){
    if(mean(yvar %in% colnames(data2)) != 1){
      varNotInData = which(!yvar %in% colnames(data2))
      varNotInDataString = paste0(yvar[varNotInData], collapse = ", ")
      warning(paste0( "Warning: These yvar are not in colnames of data2 and are ignored:\n", varNotInDataString))
      
      yvar = yvar[-varNotInData]
    }
  }
  
  
  # rename xvar automatically
  if(!xvar %in% colnames(data1) & "year" %in% colnames(data1) ){
    if( grepl("year", pattern = xvar) ){
      xvar = "year"
    }
  }
  
  if(!xvar %in% colnames(data1) & "day" %in% colnames(data1) ){
    if(grepl("day", pattern = xvar) ){
      xvar = "day"
    }
  }
  
  
  if(!xvar %in% colnames(data1) & "hour" %in% colnames(data1) ){
    if(grepl("hour", pattern = xvar) ){
      xvar = "hour"
    }
  }
  
  # verification of xvar presence
  if(!xvar %in% colnames(data1)){
    additionalMsg = ""
    if(!is.null(tableName)){
      additionalMsg = paste0(" (table: ", tableName, ")")
    }
    stop(paste0( "Warning: xvar ", xvar, " is not in colnames of data1", additionalMsg))
  }
  
  if(!is.null(data2)){
    if(!xvar %in% colnames(data2)){
      additionalMsg = ""
      if(!is.null(tableName)){
        additionalMsg = paste0(" (table: ", tableName, ")")
      }
      stop(paste0( "Warning: xvar ", xvar, " is not in colnames of data2", additionalMsg))
    }
  }
  
  
  # verification that not more than 1 value is given for a given value of xvar
  minXvar = min(unique(data1[[xvar]]))
  maxXvar = max(unique(data1[[xvar]]))
  
  for(i in unique(round(seq( minXvar, maxXvar, length.out = 20 )) )){
    len1 = dim(data1[data1[[xvar]] == i, ])[1]
    
    if(len1 > 1){
      cat("\nWarning: data1 has more than one variable for a given xvar value. It may contain multiple trees")
      break
    }
    
    
    if(!is.null(data2)){
      len2 = dim(data2[data2[[xvar]] == i, ])[1]
      
      if(len2 > 1){
        stop("Warning: data2 has more than one variable for a given xvar value")
      }
    }
  }
  
  # This is not necessary, since number of value for a given xvar is checked to be <= 1
  # The case when a table has value for only a certain range of xvar should not be rejected
  # # verification of correspondance
  # dimData1 = dim(data1[yvar])
  # dimData2 = dim(data2[yvar])
  # if(!is.null(data2) &  mean( dimData1 == dimData2 ) != 1 ){
  #   msg = paste0("Error: Data frames dimensions do not correspond.",
  #                "\ndata1: ", paste0(dimData1, collapse = ", "), 
  #                "\ndata2: ", paste0(dimData2, collapse = ", "))
  #   
  #   stop(msg)
  # }
  
  
  # DATA NAMES
  if(is.null(dataNames)){
    dataNames = c("data1", "data2")
    
    if(!is.null(suppressWarnings(data1_$name))){
      dataNames[1] = data1_$name
    }
    
    if(!is.null(suppressWarnings(data2_$name))){
      dataNames[2] = data2_$name
    }
  }
  
  if(isStandScale){
    annotation = paste0(annotation, "stand\n")
  }
  
  if(!is.null(years)){
    annotation = paste0(annotation, years, "\n")
  }
  
  
  nVar = length(yvar)
  if(nVar > 1){ # SEVERAL VARIABLES
    if(nVar < 4){
      ncol = 1
    }else{
      # nrow = ceiling(sqrt(nVar))
      ncol = 2
    }
    
    if(!is.null(nColMultiPlot)){
      ncol = nColMultiPlot
      # nrow = ceiling(nVar / ncol)
    }
    
    # a plot isnt showed when there is 7 plot with 3 columns..
    if(nVar == 7 & ncol == 3){
      ncol = 2
    }
    
    
    aPlotList = list()
    i = 1
    for(yvar_i in yvar){
      if(i == 1){ # first plot, with title
        aPlotList[[i]] = getOneSimuPlot(data1, data2, xvar, yvar_i, dataNames, 
                                        title = title, hideLabel = FALSE, annotation = annotation,
                                        type = type, alphaValue = alphaValue, isStandScale = isStandScale, dataReductionFactor = dataReductionFactor, treeLegend = treeLegend, sortYvar = sortYvar)
      }else{ # next plot
        aPlotList[[i]] = getOneSimuPlot(data1, data2, xvar, yvar_i, dataNames, 
                                        hideLabel = TRUE, type = type, alphaValue = alphaValue, isStandScale = isStandScale, dataReductionFactor = dataReductionFactor, sortYvar = sortYvar)
      }
      aPlotList[[i]] = aPlotList[[i]] + coord_cartesian
      i = i + 1
    }
    # do.call("grid.arrange", args = plotList)
    return(plot_grid(plotlist = aPlotList, ncol = ncol, align = "v", byrow = FALSE))
  }else{ # ONE VARIABLE TO PLOT
    return(getOneSimuPlot(data1, data2, xvar, yvar, dataNames, title, hideLabel = FALSE, type, alphaValue = alphaValue, isStandScale = isStandScale, dataReductionFactor = dataReductionFactor, treeLegend = treeLegend, sortYvar = sortYvar))
  }
}


# ..............................................................................

# new method
# ONE GRAPH (WARNING : for the user, use simuPlots instead)
# Show a simulation, area representation, perfect for following daily or hourly variation along year / day
# can be used for one dataset or two
getOneSimuPlot = function(data1, data2 = NULL, xvar, yvar, dataNames = NULL, title = "", 
                          hideLabel = FALSE, type = "area", alphaValue = 0.5, 
                          isStandScale = FALSE, dataReductionFactor = NULL, treeLegend = FALSE, sortYvar = FALSE, annotation = NULL){
  
  # if dataReductionFactor is given and greater than 1
  if(!is.null(dataReductionFactor)){
    if(dataReductionFactor > 1){
      
      data1 = dataReduction(data1, dataReductionFactor)
      if(!is.null(data2)){
        data2 = dataReduction(data2, dataReductionFactor)
      }
    }
  }
  
  name1 = "data1"
  name2 = "data2"
  
  if(!is.null(dataNames)){
    name1 = dataNames[1]
    name2 = dataNames[2]
  }
  
  
  
  multipleLineMode = FALSE
  if(!isStandScale & type != "bars"){
    if("idFmCell" %in% colnames(data1)){
      if(length(unique(data1$idFmCell)) == 1){
        name1 = paste0(name1, ifelse(name1 == "", "", " "), unique(data1$idFmCell) )
      }else{
        # in multipleLine mode if multiple tree are provided in one data.frame
        multipleLineMode = TRUE
        # in this mode, data2 is not represented
        data2 = NULL
        # force line in multipleLine mode
        type = "line"
      }
    }
    
    if("idFmCell" %in% colnames(data2)){
      if(length(unique(data2$idFmCell)) == 1)
        name2 = paste0(name2, ifelse(name2 == "", "", " "), unique(data2$idFmCell) )
    }
  }
  
  
  
  if(treeLegend && type != "line"){
    cat("treeLegend parameter is only for line plot. It won't be used.")
  }
  
  if(sortYvar && type != "bars"){
    cat("sortYvar parameter is only for bars plot. It won't be used.")
  }
  
  if(type == "area"){ 
    aes1 = aes_string(x = xvar, y = yvar, fill = paste0("'", name1, "'"))
    myPlot = ggplot(data1, aes1) + geom_area()
    
    if(is.null(data2)){
      myPlot = myPlot+ scale_fill_manual("",
                                         breaks = c(name1), 
                                         values = c(hsv(0.4,0.5,0.25,alphaValue + 0.25*(1 - alphaValue)))) # couleur verte, alpha moin prononcé
    }else{
      aes2 = aes_string(x = xvar, y = yvar, fill = paste0("'", name2, "'"))
      myPlot = myPlot + geom_area(data = data2, aes2) +
        scale_fill_manual("",
                          breaks = c(name1, name2), 
                          values = c(hsv(0.05,0.8,0.75,alphaValue), hsv(0.55,1,0.45,alphaValue * 0.85))) 
                          # values = c(hsv(0.55,1,0.7,alphaValue), hsv(0.05,1,0.5,alphaValue)))
    }
    
    myPlot = myPlot + guides(colour = guide_legend(override.aes = list(alpha = alphaValue)))
    
  }else if(type == "cloud"){
    dotSize = 1.5
    dotShape = 16 # see https://www.datanovia.com/en/blog/ggplot-point-shapes-best-tips/
    aes1 = aes_string(x = xvar, y = yvar, color = paste0("'", name1, "'"))
    myPlot =  ggplot(data1, aes1) + geom_point(size = dotSize, shape = dotShape)+
      guides(colour = guide_legend(override.aes = list(alpha = 1))) # the points in legend wont have transparency
    if(is.null(data2)){
      myPlot = myPlot +
        scale_color_manual("",
                           breaks = c(name1), 
                           values =  c(hsv(0.4,0.6,0.3,0.8*alphaValue + 0.8*0.5*(1 - alphaValue))))
    }else{
      aes2 = aes_string(x = xvar, y = yvar, color = paste0("'", name2, "'"))
      myPlot = myPlot + 
        geom_point(size = dotSize, shape = dotShape, data = data2, aes2)+
        scale_color_manual("",
                           breaks = c(name1, name2), 
                           values = c(hsv(0.55,1,0.7,alphaValue), hsv(0.05,1,0.6,0.8*alphaValue)))
    }
  }else if(type == "line"){
    lineWidth = 0.5
    
    if(multipleLineMode){
      aes1 = aes_string(x = xvar, y = yvar, color = "as.factor(idFmCell)")
      lineWidth = 0.5
    }else{
      aes1 = aes_string(x = xvar, y = yvar, color = paste0("'", name1, "'"))
    }
    
    myPlot = ggplot(data1, aes1) + geom_line(linewidth = lineWidth)
    
    if(is.null(data2)){
      # # no color modification
      # myPlot = myPlot+ scale_color_manual("",
      #                                    breaks = c(name1),
      #                                    values = c(hsv(0.4,0.5,0.25,alphaValue + 0.25*(1 - alphaValue)))) # couleur verte, alpha moin prononcé
    }else{
      aes2 = aes_string(x = xvar, y = yvar, color = paste0("'", name2, "'"))
      myPlot = myPlot + geom_line(data = data2, aes2) +
        scale_color_manual("",
                          breaks = c(name1, name2), 
                          values = c(hsv(0.05,1,0.8,alphaValue), hsv(0.55,1,0.5,alphaValue))) 
    }
    if(treeLegend){
      myPlot = myPlot + guides(color = guide_legend(override.aes = list(linewidth = 1)))
    }else{
      myPlot = myPlot + guides(color = FALSE) # remove legend for lines
    }
    
  }else if(type == "bars"){
    # data1 only : a bar plot with one colored bar per xvar
    # data1 and data2 : a sorted bar plot with bars colored according to the dataset 
    
    data1$dataSet = name1
    if(!is.null(data2)){
      data2$dataSet = name2
      data1[[xvar]] = paste0(data1[[xvar]], "(1)")
      data2[[xvar]] = paste0(data2[[xvar]], "(2)")
      data1 = rbind(data1, data2)
      
      sortYvar = TRUE
    }
    
    if(sortYvar){
      sortedIndex = sort(data1[[yvar]], index.return = T)$ix 
      data1[[xvar]] = factor(data1[[xvar]], levels = data1[[xvar]][sortedIndex])
      
    }else{
      data1[[xvar]] = factor(data1[[xvar]])
    }
    
    
    
    if(is.null(data2)){ # fill with white, color with xVar value
      xVarColors = gg_color_hue(length( data1[[xvar]] ))
      
      myPlot = ggplot(data1, aes_string(x = xvar, color = xvar, y = yvar, fill = "dataSet"))
      myPlot = myPlot + geom_bar(stat="identity", fill = gray(1,0.5), linewidth = 1, width = 0.75) 
      myPlot = myPlot + scale_color_manual(breaks = data1[[xvar]],
                                           values = xVarColors)
      
      myPlot = myPlot + guides(color = FALSE) # remove legend for lines
      
    }else{ # fill with dataSet, dont color
      myPlot = ggplot(data1, aes_string(x = xvar, y = yvar, fill = "dataSet"))
      myPlot = myPlot + geom_bar(stat = "identity") +
        scale_fill_hue(c=45, l=80)
    }
    
    # add value as text for each bar
    myPlot = myPlot +
      geom_text(aes_string(label=paste0("round(", yvar, ", 2)" )), vjust=1.6, color="black", size=3) +
      xlab(xvar) 
  } # end creation of plot
  
  # legend position = position of the legend anchor point in the figure from 0 (left / bottom) to 1(right / top)
  # legend.justification = position of this anchor point in the legend panel from 0 (left / bottom) to 1(right / top)
  myPlot = myPlot + theme(legend.position = c(1,1),
                          legend.justification = c(1,1),
                          legend.background = element_rect(fill = hsv(0,0,1,0.4)), # legend.key = element_rect(fill = hsv(0,0,1,0.4)),
                          legend.margin = margin(-6, 0, 0, 0),
                          # legend.box.margin = margin(-5, 5, 0, 0),
                          legend.title = element_blank(),
                          legend.key.size = unit(0.15, 'in'),
                          legend.text = element_text(margin = margin(l = -5, unit = "pt")))
  
  if(title != "")
    myPlot = myPlot + labs(title = title)
  if(hideLabel){
    myPlot = myPlot + theme(legend.position= "none")
  }else{
    # xpos = c(-Inf,-Inf,Inf,Inf),
    # ypos =  c(-Inf, Inf,-Inf,Inf),
    # annotateText = c("Text","tExt","teXt","texT"),
    # hjustvar = c(0,0,1,1) ,
    # vjustvar = c(0,1.0,0,1))
    annotations <- data.frame(
      xpos = c(-Inf),
      ypos =  c(Inf),
      annotateText = c(annotation),
      hjustvar = c(-0.1) ,
      vjustvar = c(1))
    myPlot = myPlot + geom_text(data = annotations, aes(x=xpos,y=ypos,hjust=hjustvar,
                                                        vjust=vjustvar,label=annotateText), color = hsv(0,0,0.35), size = 3)
  }
  
  # AXIS
  # move the y axis title to the top of the graph
  labelOpacity = 0.75
  labelSize = 10
  
  thisTitle = yvar
  if(title != ""){
    thisTitle = paste0(yvar, " / / ", title)
  }
  
  myPlot = myPlot + ggtitle(thisTitle) +
    theme(plot.title = element_text(hjust = 0, vjust = -1,  size = labelSize, color = gray(0, labelOpacity)),
          plot.margin = rep(grid::unit(0,"in"),4),
          axis.title.y = element_blank(),
          axis.title.x = element_text(color = gray(0, labelOpacity), size = labelSize))
  
  myPlot = myPlot + geom_abline(intercept = 0, slope = 0, linewidth = 0.25, color = hsv(0,0,0, 0.5))
  return(myPlot)
}




# ..............................................................................
# SINGLE VALUES COMPARISON -----------------------------------------------------

# show comparison of two simulations for one value of one or several attributes
# parameters similar to simuPlots
showSimuValue_comparison = function(data1, data2 = NULL, yvar, data1name, data2name = NULL, title = ""){
  nVar = length(yvar)
  if(nVar > 1){
    if(nVar < 4){
      nrow = nVar
    }else{
      nrow = ceiling(sqrt(nVar))
    }
    aPlotList = list()
    i = 1
    for(yvar_i in yvar){
      if(i == 1){
        aPlotList[[i]] = getSimuValue_comparison(data1, data2, yvar_i, data1name, data2name, title, hideLabel = FALSE)
      }else{
        aPlotList[[i]] = getSimuValue_comparison(data1, data2, yvar_i, data1name, data2name, hideLabel = TRUE)
      }
      i = i + 1
    }
    # do.call("grid.arrange", args = plotList)
    return(plot_grid(plotlist = aPlotList, nrow = nrow, align = "v", byrow = FALSE))
  }else{
    return(getSimuValue_comparison(data1, data2, yvar, data1name, data2name, title))
  }
}





# (WARNING : method showSimuValue_comparison is preferred for the user)
# show comparison of two simulations for one value
getSimuValue_comparison = function(data1, data2 = NULL, yvar, data1name, data2name = NULL, title = "", hideLabel = FALSE){
  # palette
  colValues = (0.45+0.32*seq(0,1, length.out = 2) )%%1 ; palette1 = hsv(colValues, 0.5, 0.8, 0.75)
  bindedData = bind_rows(name1 = data1[,yvar], name2 = data2[,yvar], .id = "dataName")
  bindedData[bindedData$dataName == "name1", ]$dataName = data1name
  bindedData[bindedData$dataName == "name2", ]$dataName = data2name
  xvar = "dataName"
  aesSingle = aes_string(x = xvar, y = yvar, fill = xvar)
  myPlot = ggplot(bindedData, aesSingle) + geom_col()+
    scale_fill_manual("n layer", values = palette1) +
    guides(fill = FALSE)+ # remove legend for fill
    xlab('Simulation')
  if(title != "")
    myPlot = myPlot + labs(title = title)
  if(hideLabel)
    myPlot = myPlot + theme(legend.position= "none")
  return(myPlot)
}






# Diverging bar plot, one variable per simulation
# see http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html#Diverging%20Bars
divergingPlot = function(simuTable, xvar, yvar){
  
  simuTable$yvar_sign = ifelse(simuTable[[yvar]] > 0, "above", "below")
  
  simuTable <- simuTable[order(simuTable[[yvar]]), ]  # sort
  simuTable[[xvar]] <- factor(simuTable[[xvar]], levels = simuTable[[xvar]])  # convert to factor to retain sorted order in plot.
  
  plot = ggplot(simuTable, aes_string(x=xvar, y=yvar, label=yvar)) + 
    geom_bar(stat='identity', aes(fill=yvar_sign), width=.5)  +
    scale_fill_manual(name="Mileage", 
                      labels = c("Above Average", "Below Average"), 
                      values = c("above"="#00ba38", "below"="#f8766d")) + 
    coord_flip()+ theme(legend.position= "none")
  
  return(plot)
}







# ..............................................................................
# HOURLY PARAMETERS -----------------------------------------------------------

# Variation throughout one day of a hourly variables.
# The data must have a "jour" column.
# Default layer is 5.
showVariationInAday = function(data1_, data2_ = NULL, tableName = NULL, dataNames = NULL, aDay, aLayer = 5, xvar, yvar, type = "cloud", alphaValue = 0.5, coord_cartesian = NULL){
  
  if(!is.null(tableName)){
    data1 = data1_[[tableName]]
    data2 = data2_[[tableName]]
  }else{
    data1 = data1_
    data2 = data2_
  }
  
  if("k" %in% names(data1)){
    data1 = subset(data1, k == aLayer)
    if(!is.null(data2)){
      data2 = subset(data2, k == aLayer)
    }
  }else if("layer" %in% names(data1)){
    data1 = subset(data1, layer == aLayer)
    if(!is.null(data2)){
      data2 = subset(data2, layer == aLayer)
    }
  }
  
  # different translation allowed
  if("jour" %in% names(data1)){
    data1 = subset(data1, jour == aDay)
    if(!is.null(data2)){
      data2 = subset(data2, jour == aDay)
    }
  }else if("day" %in% names(data1)){
    data1 = subset(data1, day == aDay)
    if(!is.null(data2)){
      data2 = subset(data2, day == aDay)
    }
  }
  
  # reaffect modified table after filtering
  if(!is.null(tableName)){
    data1_[[tableName]] = data1
    data2_[[tableName]] = data2
  }else{
    data1_ = data1
    data2_ = data2
  }
  
  simuPlots(data1_, data2_, tableName = tableName,
            dataNames = dataNames,
            xvar = xvar, yvar = yvar, type = type, alphaValue = alphaValue, coord_cartesian = coord_cartesian)+ 
    annotate("label", x = 0.025, y = 0.96, label = paste0("k: ", aLayer), color = rgb(0,0.3,0,0.75), alpha = 0.5)
}

# Variation throughout the year of hourly variables. 
# One hour of interest must be defined.
# The data must have a "heure" column.
# aLayer is ignored if table does not contain a layer column ("k"), Default layer is 5.
showOneHourEveryDay = function(data1_, data2_ = NULL, tableName = NULL, dataNames = NULL, anHour, aLayer = 5, xvar, yvar, type = "cloud",
                               alphaValue = 0.5, coord_cartesian = NULL){
  
  if(!is.null(tableName)){
    data1 = data1_[[tableName]]
    data2 = data2_[[tableName]]
  }else{
    data1 = data1_
    data2 = data2_
  }
  
  # reduction to one layer
  if("k" %in% names(data1)){
    data1 = subset(data1, k == aLayer)
    if(!is.null(data2)){
      data2 = subset(data2, k == aLayer)
    }
  }else if("layer" %in% names(data1)){
    data1 = subset(data1, layer == aLayer)
    if(!is.null(data2)){
      data2 = subset(data2, layer == aLayer)
    }
  }
  
  # different translation allowed
  if("hour" %in% names(data1)){
    data1 = subset(data1, hour == anHour)
    if(!is.null(data2)){
      data2 = subset(data2, hour == anHour)
    }
  }else if("heure" %in% names(data1)){
    data1 = subset(data1, heure == anHour)
    if(!is.null(data2)){
      data2 = subset(data2, heure == anHour)
    }
  }
  
  # reaffect modified table after filtering
  if(!is.null(tableName)){
    data1_[[tableName]] = data1
    data2_[[tableName]] = data2
  }else{
    data1_ = data1
    data2_ = data2
  }
  
  
  if(is.null(coord_cartesian))
    coord_cartesian = NULL
  
  
  simuPlots(data1_, data2_, tableName = tableName,
            dataNames = dataNames,
            xvar = xvar, yvar = yvar, type = type, alphaValue = alphaValue, coord_cartesian = coord_cartesian) + 
    annotate("label", x = 0.025, y = 0.96, label = paste0("k: ", aLayer), color = rgb(0,0.3,0,0.75), alpha = 0.5)
}




# ..............................................................................
# OTHER PLOTS ------------------------------------------------------------------


# Box plot for all years of all sites for a defined yearly variable (yvar)
# composition_align
plotYearlySites = function(simuList, yvar = "GPP", composition_align = FALSE){
  
  df = tibble(year = 0, composition = "", triplet = "", site = "a", yvariable = 0, option = "", .rows = 0)
  
  for(site in names(simuList)){
    simu = simuList[[site]]
    for(i in 1:dim(simu$yearlyResults)[1]){
      line = simu$yearlyResults[i, ]
      
      if(composition_align){
        composition = line$composition
        triplet = line$triplet
      }else{
        composition = ""
        triplet = ""
      }
      
      df = add_row(df, year = line$year, site = site, yvariable = line[[yvar]], composition = composition, triplet = triplet)
    }
  }
  
  years = unique(simu$yearlyResults$year)
  years_bis = seq(min(years), max(years), by = 2)
  
  if(composition_align){
    facet = facet_grid( vars(triplet), vars(composition) )
  }else{
    facet = facet_wrap(. ~ site)
  }
  
  plot = ggplot(df, aes(x = as.factor(year), y = yvariable)) + geom_boxplot() + facet +
    ylab(yvar) + 
    theme(axis.text.x = element_text(angle=30, hjust = 1)) +
    scale_x_discrete(breaks = years_bis )
  return(plot)
}




# ..............................................................................
# SHOW TABLES ------------------------------------------------------------------

# show two vectors (of same nature) side by side
showVectors = function(v1, v2 = NULL, names = NULL){
  bindTable = data.frame(v1) # if only one vector is proposed
  if(!is.null(v2)){ # if two vectors
    bindTable[2,] = v2
    if("date" %in% names(bindTable)){
      bindTable$date = format(bindTable$date, format="%B %d %Y")
    }
    bindTable[3,] = v1 == v2
    if(is.null(names)){ # change names
      rownames(bindTable) = c("v1", "v2", "isEqual")
    }else{
      rownames(bindTable) = c(names[1], names[2], "isEqual")
    }
  }
  # several identical path to have different variable names
  randomValue = runif(1)
  prop = 1/20
  if(randomValue < prop * 1){
    bindTable1 = t(bindTable)
    View(bindTable1) # view transpose of the binded inventories
  }else if (randomValue < prop * 2){
    bindTable2 = t(bindTable)
    View(bindTable2) # view transpose of the binded inventories
  }else if (randomValue < prop * 3){
    bindTable3 = t(bindTable)
    View(bindTable3) # view transpose of the binded inventories
  }else if (randomValue < prop * 4){
    bindTable4 = t(bindTable)
    View(bindTable4) # view transpose of the binded inventories
  }else if (randomValue < prop * 5){
    bindTable5 = t(bindTable)
    View(bindTable5) # view transpose of the binded inventories
  }else if (randomValue < prop * 6){
    bindTable6 = t(bindTable)
    View(bindTable6) # view transpose of the binded inventories
  }else if (randomValue < prop * 7){
    bindTable7 = t(bindTable)
    View(bindTable7) # view transpose of the binded inventories
  }else if (randomValue < prop * 8){
    bindTable8 = t(bindTable)
    View(bindTable8) # view transpose of the binded inventories
  }else if (randomValue < prop * 9){
    bindTable9 = t(bindTable)
    View(bindTable9) # view transpose of the binded inventories
  }else if (randomValue < prop * 10){
    bindTable10 = t(bindTable)
    View(bindTable10) # view transpose of the binded inventories
  }else if (randomValue < prop * 11){
    bindTable11 = t(bindTable)
    View(bindTable11) # view transpose of the binded inventories
  }else if (randomValue < prop * 12){
    bindTable12 = t(bindTable)
    View(bindTable12) # view transpose of the binded inventories
  }else if (randomValue < prop * 13){
    bindTable13 = t(bindTable)
    View(bindTable13) # view transpose of the binded inventories
  }else if (randomValue < prop * 14){
    bindTable14 = t(bindTable)
    View(bindTable14) # view transpose of the binded inventories
  }else if (randomValue < prop * 15){
    bindTable15 = t(bindTable)
    View(bindTable15) # view transpose of the binded inventories
  }else if (randomValue < prop * 16){
    bindTable16 = t(bindTable)
    View(bindTable16) # view transpose of the binded inventories
  }else if (randomValue < prop * 17){
    bindTable17 = t(bindTable)
    View(bindTable17) # view transpose of the binded inventories
  }else if (randomValue < prop * 18){
    bindTable18 = t(bindTable)
    View(bindTable18) # view transpose of the binded inventories
  }else if (randomValue < prop * 19){
    bindTable19 = t(bindTable)
    View(bindTable19) # view transpose of the binded inventories
  }else if (randomValue < prop * 20){
    bindTable20 = t(bindTable)
    View(bindTable20) # view transpose of the binded inventories
  }
  
}

# if vector are of different nature / size
showVectors2 = function(v1, v2, names = NULL){
  lmax = max(length(v1), length(v2))
  if(length(v1) < lmax){
    v1 = c(v1, rep(NA, lmax - length(v1)))
    bindTable = data.frame(v2)
  }
  if(length(v2) < lmax){
    v2 = c(v2, rep(NA, lmax - length(v2)))
    bindTable = data.frame(v1)
  }
  bindTable[1,] = names(v1)
  bindTable[2,] = v1
  bindTable[3,] = names(v2)
  bindTable[4,] = v2
  bindTable = t(bindTable)
  if(is.null(names)){
    colnames(bindTable) = c("Nv1", "v1", "Nv2", "v2")
  }else{
    colnames(bindTable) = c(names[1], names[2], "isEqual")
  }
  rownames(bindTable) = NA
  View(bindTable) # view transpose of the binded inventories
}

# ..............................................................................

# Show one inventory or two equivalent inventories side by side
showInv = function(simu1, simu2 = NULL, names = NULL){
  showVectors(simu1$inventory, simu2$inventory, names)
}

# Show inventories (used that function if they are different in length)
showInv2 = function(simu1, simu2, names = NULL){
  showVectors2(simu1$inventory, simu2$inventory, names)
}

# ..............................................................................

showYearlyResults = function(simu1, simu2 = NULL, names = NULL){
  showVectors(simu1$yearlyResults, simu2$yearlyResults, names)
}

# ..............................................................................


# Show one or two soils side by side
showSoilFirstcell = function(simu1, simu2 = NULL, names = NULL){
  showVectors(simu1$soilFirstCell, simu2$soilFirstCell, names)
}

# Show fmSettings for one or two table side by side
showFmSettings = function(simu1, simu2 = NULL, names = NULL){
  showVectors(simu1$fmSettings, simu2$fmSettings, names)
}

# Show PDGInitialParameters for one or two table side by side
showPDGInitialParameters = function(simu1, simu2 = NULL, names = NULL){
  showVectors(simu1$PDGInitialParameters, simu2$PDGInitialParameters, names)
}


# Show the daily mean value for essential hydric fluxes
# for one or two simulation
showWaterBudget = function(simu1, simu2 = NULL, names = NULL){
  
  table1 = getWaterBudgetTable(simu1)
  
  table2 = NULL
  if(!is.null(simu2)){
    table2 = getWaterBudgetTable(simu2)
  }
  
  
  showVectors(table1, table2, names)
}


# return a table with the daily mean value for essential hydric fluxes
getWaterBudgetTable = function(aSimu, meanDay = NULL, varDay = NULL){
  
  aSimuDailyResults = aSimu$dailyResults
  
  if(!is.null(meanDay)){
    aSimuDailyResults = subset(simu1$dailyResults, abs(day - meanDay) < varDay)
  }
  
  table = tibble(TR_mean = mean(aSimuDailyResults$TR), 
                 ETRsoil_mean = mean(aSimuDailyResults$ETRsoil),
                 evapoCan_mean = mean(aSimuDailyResults$ETRcan - aSimuDailyResults$TR),
                 drainage_mean = mean(aSimuDailyResults$drainage),
                 Rain_mean = mean(aSimuDailyResults$Rain) ) 
  
  return(table)
}







# ..............................................................................
# SECONDARY METHODS ------------------------------------------------------------


# extract GPP of first day, last day, and mean GPP of vegetation period
# /!\ work with dailyResults !! need fmOutput >= 2
getBSSautoLine = function(simulation, targetYear = NULL){
  
  # check
  if(is.null(targetYear)){
    targetYear = unique(simulation$dailyResults$year)
    
    if(length(targetYear) > 1){
      stop("getBSSautoLine() : targetYear was not defined and dailyResults has more than one year")
    }
  }
  
  invSplit = strsplit(simulation$inventory$originalinventory, "/" )[[1]]
  inventoryName = invSplit[length(invSplit)]
  inventoryName = strsplit( inventoryName, "[.]" ) [[1]][1]
  
  LAImaxThisYear = subset(simulation$yearlyResults, year == targetYear)$LAImaxThisYear
  species = subset(simulation$yearlyResults, year == targetYear)$species
  
  simulationDaily = subset(simulation$dailyResults, year == targetYear)
  
  firstDay = min(simulationDaily$day)
  lastDay = max(simulationDaily$day)
  
  vegDays = simulationDaily$day[ abs(simulationDaily$day - 140) < 30 ] # simulationDaily$day > simulation$yearlyResults$endOfLeafAreaGrowth & simulationDaily$day < simulation$yearlyResults$endLeaf 
  
  line = list(inventory = inventoryName,
              year = targetYear, 
              BSSday1 = simulationDaily$BSS[simulationDaily$day == firstDay],
              BSSlastDay = simulationDaily$BSS[simulationDaily$day == lastDay],
              BSSmeanVegSeason = mean( simulationDaily$BSS[simulationDaily$day %in% vegDays]),
              LAImax = LAImaxThisYear,
              species = species,
              RU = simulation$yearlyResults$RU)
  
  return(line)
}


getLastYearOfSimulation = function(simul){
  line = tibble(lastYear = max(simul$yearlyResults$year))
  return(line)
}




## Export data of interest (defined by myfunc, hwich coulb be eg, getBSSautoLine) from a list of simulation
# one line per simulation
# /!\ work with dailyResults !! need fmOutput >= 2
getAutoTableOnSimuList = function(simuList, myfunc){
  
  autoAnalysis = tibble(myfunc(simuList[[1]]), .rows = 0)
  
  for(inventory in names(simuList)){
    simu = simuList[[inventory]]
    
    autoAnalysis = bind_rows(autoAnalysis, myfunc(simu))
    
  } # end loop on inventories
  
  rownames(autoAnalysis) = names(simuList)
  
  if(is.numeric(autoAnalysis$species[1])){
    autoAnalysis$species = ifelse(autoAnalysis$species == 0, "stand", ifelse(autoAnalysis$species == 1, "Abies alba", ifelse(autoAnalysis$species == 3, "Fagus sylvatica", NA)))
  }
  
  
  
  return(autoAnalysis)
}













# For each simulation of a simuList, create a new variable that represent the mean value of a given variable inside a given class
# tableName     table that will be affected by the "pixellation"
# variableName  variable that will be affected by the "pixellation". For the moment the method is designed for daily table only
# categorySize  size of the calss over which a mean will be done
pixellateVariableSimuList = function(standSimuList, tableName, variableName, categorySize = 5){
  
  
  
  for(inventory in names(standSimuList)){
    asimu = standSimuList[[inventory]]
    
    tableName = identifyTableName(tableName, names(asimu))
    
    asimuTable = asimu[[tableName]]
    
    asimuTable = pixellateVariableOnDailyTable(asimuTable, variableName, categorySize)
    
    standSimuList[[inventory]][[tableName]] = asimuTable
    
  }
  return(standSimuList)
}


# Create a new variable that represent the mean value of a given variable inside a given class
# variableName  variable that will be affected by the "pixellation". For the moment the method is designed for daily table only
# categorySize  size of the calss over which a mean will be done
pixellateVariableOnDailyTable = function(asimuTable, variableName_, categorySize){
  
  for(variableName in variableName_){
    variableNameBis = paste0(variableName, "_pix", categorySize)
    
    
    yearsSimu = unique(asimuTable$year)
    idsSimu = unique(asimuTable$idFmCell)
    
    # new column
    asimuTable[[variableNameBis]] = -1
    
    for(ayear in yearsSimu){
      for(anId in idsSimu){
        asimuTableOneYear = subset(asimuTable, year == ayear & idFmCell == anId)
        
        alldays = asimuTableOneYear$day
        nCategories = length(alldays)%/%categorySize +1
        alldaysCategorized = alldays[ cut(alldays, nCategories) ]
        categories = unique(alldaysCategorized)
        
        for(daycat in categories){
          meanVariable = mean(asimuTableOneYear[[variableName]][alldaysCategorized == daycat])
          asimuTableOneYear[[variableNameBis]][alldaysCategorized == daycat] = meanVariable
        }
        
        
        testForError = TRUE
        if(testForError){
          error = mean(asimuTableOneYear[[variableNameBis]]) - mean(asimuTableOneYear[[variableName]])
          
          if(error > 1e-10){
            stop("Pixellation did not worked")
          }
        }
        
        asimuTable[asimuTable$year == ayear & asimuTable$idFmCell == anId, ] = asimuTableOneYear
      } # end loop on ids
    } # end loop on years
  } # end loop on variableName_
  
  return(asimuTable)
}




# reduce the number of point of a table by a given factor
# useful for plot
dataReduction = function(dataTable, dataReductionFactor){
  
  lengthDataTable = dim(dataTable)[1]
  xkept = seq(1, lengthDataTable, by = dataReductionFactor)
  
  # add the last row to keep table boundaries
  if(!lengthDataTable %in% xkept){
    xkept = c(xkept, lengthDataTable)
  }
  
  dataTable = dataTable[xkept, ]
  return(dataTable)
}


# find the real value of tableName if only a prefix is given
# eg, give "dailyResults" for "dai"
identifyTableName = function(tableName, possibleNames = NULL){
  
  if(!tableName %in% possibleNames){
    firstLetter = substr(tableName, 1, 1)
    searchPattern = grepl(possibleNames, pattern = tableName)
    
    # check for yearly, daily and hourly tag
    if(firstLetter == "d"){
      tableName = "dailyResults"
    }else if(firstLetter == "h"){
      tableName = "hourlyResults"
    }else if(firstLetter == "y"){
      tableName = "yearlyResults"
    }else if(grepl(pattern = "day", tableName) | grepl(pattern = "dai", tableName) | grepl(pattern = "jour", tableName)){
      tableName = "dailyResults"
    }else if(grepl(pattern = "hour", tableName) | grepl(pattern = "heur", tableName) | grepl(pattern = "hor", tableName)){
      tableName = "hourlyResults"
    }else if(grepl(pattern = "year", tableName) | grepl(pattern = "ann", tableName)){
      tableName = "yearlyResults"
    }else if(TRUE %in% searchPattern){
      index = which(searchPattern)[1]
      tableName = possibleNames[index]
    }
  }
  
  return(tableName)
}





# reproduce ggplot basic color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}




saveGgPlot = function(plot = NULL, plotfolderPlot, plot_height = 960, plot_width = NULL, ratio = 4/3, scale = 1, fileName = NULL, fileSuffix = ".pdf"){
  
  
  if(!dir.exists(folderPlot)){
    dir.create(folderPlot, recursive = T)
  }
  
  if(is.null(fileName)){
    index_filePath = paste0(folderPlot,"plot_index.txt")
    
    # read index file or set to 1
    if(file.exists(index_filePath)){
      increment_plot = as.integer(read_file(index_filePath))
    }else{
      increment_plot = 1
    }
    
    fileName = paste0("plot_", increment_plot)
    increment_plot = increment_plot + 1
    
    write_file(x = paste0(increment_plot), file = index_filePath, append = F)
  }
  
  imagePath = paste0(folderPlot, fileName, fileSuffix)
  
  if(is.null(plot_width)){
    plot_width = plot_height * ratio
  }else{
    plot_height = plot_width / ratio
  }
  
  ggsave(imagePath, height = plot_height, width = plot_width, 
         dpi = 100, units = "px", scale = scale)
}

# plot the last ggplot, see saveGgPlot
saveLastGgPlot = function(plotfolderPlot, plot_height = 960, plot_width = NULL, ratio = 4/3, scale = 1, fileName = NULL, fileSuffix = ".pdf"){
  saveGgPlot(plotfolderPlot = plotfolderPlot, plot_height = plot_height, plot_width = plot_width, ratio = ratio, scale = scale, fileName = fileName, fileSuffix = fileSuffix)
}







# ..............................................................................
# Simulation analysis ----------------------------------------------------------



# function: returns a numerical index based on num_arbre and num_tige from GMAP data (same method to define idFmCell in GMAPtoPDG.R)
getNumIndex = function(numArbre, numTige){
  
  # -> as.numeric("11e") gives numeric 11 ..
  # Then, isCoppice = suppressWarnings( is.na(as.numeric(numTige)) ) # with coppice, numTige as a letter
  # is replaced by : 
  numbers_only <- function(x) !grepl("\\D", x)
  isCoppice = !numbers_only(numTige)
  
  # define a numerical index for each stem (tige)
  num_index = rep(-1, length(numTige))
  for(i in 1:length(numTige)){
    aNumTige = numTige[i]
    aNumArbre = numArbre[i]
    
    if(isCoppice[i]){
      tige_letter = substr(aNumTige, start = nchar(aNumArbre) + 1, stop = nchar(aNumTige))
      
      # remove space if found in tige_letter
      if(grepl(tige_letter, pattern = " ")){
        tige_letter = gsub(x = tige_letter, pattern = " ", replacement = "")
      }
      
      tige_letter_as_index = utf8ToInt(tige_letter) - utf8ToInt('a') + 1 # from letter to integer : index from 1 to XX
      
      num_index[i] = aNumArbre + tige_letter_as_index * 1000
    }else{
      num_index[i] = aNumArbre
    }
  }
  
  return(num_index)
}



# MAKE A TREE-YEAR TABLE FROM SIMULATION
# simulationList is a list of simulation from importSimuList at individual scale
makeTreeYearTable = function(simulationList){
  
  # find table dimension 
  nYearTreeTot = 0
  nYearSimu = 0
  for(code_site in names(simulationList)){
    yRes = simulationList[[code_site]]$yearlyResults
    nYearTreeCodeSite = length(unique(yRes$idFmCell)) * length(unique(yRes$year)) 
    
    if(code_site != names(simulationList)[1] & nYearSimu != length(unique(yRes$year))){
      stop("not the same number of year for each tree.")
    }
    
    nYearSimu = length(unique(yRes$year)) 
    nYearTreeTot = nYearTreeTot + nYearTreeCodeSite
  }
  
  # create table 
  treeYearTable = tibble( treeGlobalId = "",
                       year = 0,
                       GPPm2_sim = 0, # gC/m2/yr
                       GPPabs_sim = 0, # gC/yr
                       NPPabs_sim = 0, # gC/yr
                       RautoAbs_sim = 0, # gC/yr
                       absorbedPAR_m2_sim = 0,
                       absorbedPAR_abs_sm = 0,
                       BAI_sim = 0, # BAI from mm2 to cm2
                       WVI_sim = 0, # in m3
                       .rows = nYearTreeTot) 
  
  
  index_line = 1
  nIter = nYearTreeTot
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  
  for(i in 1:length(names(simulationList))){ #  i = 1
    nameSimu = names(simulationList)[i]
    
    # information about simulation stand
    code_site = getCodeSiteFromSimulationName(nameSimu)
    
    # ONE LINE PER TREE-YEAR
    yearlyResultsSimu = simulationList[[nameSimu]]$yearlyResults
    nTreeYear = dim(yearlyResultsSimu)[1]
    for(j in 1:nTreeYear){ # j = 1
      yearlyResultsSimuLine = yearlyResultsSimu[j, ]
      
      treeYearTable[index_line, 1] = paste0(code_site, "_", yearlyResultsSimuLine$idFmCell)
      treeYearTable[index_line, 2] = yearlyResultsSimuLine$year
      treeYearTable[index_line, 3] = yearlyResultsSimuLine$GPP # gC/m2/yr
      treeYearTable[index_line, 4] = yearlyResultsSimuLine$GPP * yearlyResultsSimuLine$crownProjectionThisYear
      treeYearTable[index_line, 5] = yearlyResultsSimuLine$NPP * yearlyResultsSimuLine$crownProjectionThisYear
      treeYearTable[index_line, 6] = yearlyResultsSimuLine$Rauto * yearlyResultsSimuLine$crownProjectionThisYear
      treeYearTable[index_line, 7] = yearlyResultsSimuLine$vegPAR_yearlyMJm2
      treeYearTable[index_line, 8] = yearlyResultsSimuLine$vegPAR_yearlyMJm2 * yearlyResultsSimuLine$crownProjectionThisYear
      treeYearTable[index_line, 9] = yearlyResultsSimuLine$BAI * 1e-2 # BAI from mm2 to cm2
      treeYearTable[index_line, 10] = yearlyResultsSimuLine$WVI # in m3
      
      index_line = index_line + 1
      setTxtProgressBar(pb, index_line) # change progress bar
    }
  }
  close(pb)
  
  
  
  
  # Complete the table with stand attributes
  treeYearTable$site = vapply(treeYearTable$treeGlobalId, function(x) strsplit(x, "_")[[1]][1], FUN.VALUE = "example")
  treeYearTable$composition = vapply(treeYearTable$treeGlobalId, function(x) strsplit(x, "_")[[1]][3], FUN.VALUE = "example")
  treeYearTable$triplet = vapply(treeYearTable$treeGlobalId, function(x){ split = strsplit(x, "_")[[1]] ; return(paste(split[1], split[2], split[4],sep = "_") )}, FUN.VALUE = "example")
  treeYearTable$code_site = vapply(treeYearTable$treeGlobalId, function(x){ split = strsplit(x, "_")[[1]] ; return(paste(split[1:4], collapse = "_"))}, FUN.VALUE = "example")
  
  
  ## check the number of line of the table
  nTree = length(unique(treeYearTable$treeGlobalId))
  nYear = length(unique(treeYearTable$year))
  dim(treeYearTable)
  nTree * nYear == dim(treeYearTable)[1]
  
  # find which tree have a different number of years
  treeid_simulated = unique(treeYearTable$treeGlobalId)
  treeidtoremoveList = c()
  for(treeid in treeid_simulated){
    treeYearTableSubset = subset(treeYearTable, treeGlobalId == treeid)
    nYear = length(treeYearTableSubset$year)
    if( nYear != nYearSimu){
      treeidtoremoveList = c(treeidtoremoveList, treeid)
      print(treeid)
      print(nYear)
    }
  }
  rm(treeYearTableSubset, treeid, nYear)
  
  # remove tree with too much or less years from table
  if(length(treeidtoremoveList) > 0){
    treeYearTable = subset(treeYearTable, !treeGlobalId %in% treeidtoremoveList)
  }
  
  return(treeYearTable)
}# end makeTreeYearTable








# Make a table with tree-year lines containing BAI measured (dGMAP_dendro) and ADD INFO ON TREES (horsprod)
# one line per year of measure in dGMAP_dendro and per stem found in dGMAP_horsprod
makeGMAP_TreeYearTable = function(dGMAP_dendro, dGMAP_horsprod){
  
  # find table dimension 
  yearsLoop = 1993:2013
  nTreeYearTot = length(dGMAP_horsprod$treeGlobalId) * length(yearsLoop)
  
  # create table
  treeYearTableGMAP = tibble( treeGlobalId = "",
          code_site = "",
          year = 0,
          BAI_mes = 0,
          WVI_mes = 0,
          species = "",
          dbh = 0,
          hauteur = 0,
          hcb = 0,
          span = 0,
          goodCarrots = FALSE,
          .rows = nTreeYearTot)
  
  index_line = 1
  nIter = length(dGMAP_horsprod$treeGlobalId)
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  
  # for all line of dGMAP_dendro
  for(i in 1:nIter){ # loop on trees
    dGMAP_horsprod_tree = dGMAP_horsprod[i,]
    treeGlobalId = dGMAP_horsprod_tree$treeGlobalId
    
    code_site = dGMAP_horsprod_tree$code_site
    species = dGMAP_horsprod_tree$essence
    dbh = dGMAP_horsprod_tree$circonference / pi # in cm
    hauteur = dGMAP_horsprod_tree$htot
    hcb = dGMAP_horsprod_tree$h_base_houppier
    span = dGMAP_horsprod_tree$htot - dGMAP_horsprod_tree$h_base_houppier
    
    # has the tree been carroted ?
    goodCarrots = FALSE
    if(treeGlobalId %in% dGMAP_dendro$treeGlobalId){
      selection = which(dGMAP_dendro$treeGlobalId == treeGlobalId)[1] # in case of doublon
      dGMAP_dendro_tree = dGMAP_dendro[selection, ]
      goodCarrots = TRUE
    }
    
    # for all growth years of a line of dGMAP_dendro
    for(year in yearsLoop){
      
      BAI_mes = -999
      WVI_mes = -999 # wood volume increment in m3
      
      if(goodCarrots){
        colName = paste0("X", year)
        BAI_mes = dGMAP_dendro_tree[[colName]]
        WVI_mes = BAI_mes * 0.01**2 * hauteur # BAI from cm2 to m2
      }
      
      treeYearTableGMAP[index_line, 1] = treeGlobalId
      treeYearTableGMAP[index_line, 2] = code_site
      treeYearTableGMAP[index_line, 3] = year
      treeYearTableGMAP[index_line, 4] = BAI_mes
      treeYearTableGMAP[index_line, 5] = WVI_mes
      treeYearTableGMAP[index_line, 6] = species
      treeYearTableGMAP[index_line, 7] = dbh
      treeYearTableGMAP[index_line, 8] = hauteur
      treeYearTableGMAP[index_line, 9] = hcb
      treeYearTableGMAP[index_line, 10] = span
      treeYearTableGMAP[index_line, 11] = goodCarrots
      index_line = index_line + 1
      
      
    } # end loop on years
    setTxtProgressBar(pb, i) # change progress bar
  } # end loop on trees
  close(pb)
  
  return(treeYearTableGMAP)
}



# EXTRAPOLATE measured BAI of NON-MEASURED TREES
# For all groups of same year, stand and species, estimate mean(BAI/BA), then compute mean(BAI/BA) * BA for trees without BAI measure
extrapoleMissingMeasuredBAI = function(treeYearTableGMAP){
  treeYearTableGMAP$basalArea = pi * (treeYearTableGMAP$dbh/2)**2 # in cm2
  
  code_site_list = unique(treeYearTableGMAP$code_site)
  nIter = length(code_site_list)
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=")
  
  # for all stand
  for(i in 1:nIter){
    a_code_site = code_site_list[i]
    treeYearTable_stand = subset(treeYearTableGMAP, code_site == a_code_site)
    
    # for all year
    for(ayear in unique(treeYearTable_stand$year)){
      treeYearTable_stand_year = subset(treeYearTable_stand, year == ayear)
      
      # for all species
      for(aSpecies in unique(treeYearTable_stand_year$species)){
        treeYearTable_stand_year_species = subset(treeYearTable_stand_year, species == aSpecies)
        treeYearTable_stand_year_species_measured = subset(treeYearTable_stand_year, species == aSpecies & goodCarrots == TRUE)
        
        # si aucun individu de l'espece n'est present, on utilise le rapport BAI/BA des individus du peuplement cette anéne
        if(dim(treeYearTable_stand_year_species_measured)[1] == 0){
          treeYearTable_stand_year_species_measured = subset(treeYearTable_stand_year, goodCarrots == TRUE)
        }
        
        meanBAIonBA = mean(treeYearTable_stand_year_species_measured$BAI_mes / treeYearTable_stand_year_species_measured$basalArea)
        
        selection = treeYearTableGMAP$code_site == a_code_site & treeYearTableGMAP$year == ayear & treeYearTableGMAP$species == aSpecies & treeYearTableGMAP$goodCarrots == FALSE
        treeYearTableGMAP[selection, ]$BAI_mes = treeYearTableGMAP[selection, ]$basalArea * meanBAIonBA
        treeYearTableGMAP[selection, ]$WVI_mes = treeYearTableGMAP[selection, ]$BAI_mes * 0.01**2 * treeYearTableGMAP[selection, ]$hauteur # BAI from cm2 to m2
      }
    }
    setTxtProgressBar(pb, i) # change progress bar
  }
  close(pb)
  rm(treeYearTable_stand, treeYearTable_stand_year, treeYearTable_stand_year_species, treeYearTable_stand_year_species_measured)
  return(treeYearTableGMAP)
}







# Based on a treeYearTable, makes a table with one line per stand-year
# Sums up individual BAI for each stand-year
makeStandYearTable = function(treeYearTable){
  code_site_list = unique(treeYearTable$code_site)
  
  nStandYearTot = length(code_site_list) * length(unique(treeYearTable$year))
  standYearTable = tibble( code_site = "",
                           year = 0,
                           site = "",
                           composition = "",
                           triplet = "",
                           GPPabs_sim = 0, # gC/yr
                           NPPabs_sim = 0, # gC/yr
                           RautoAbs_sim = 0, # gC/yr
                           BAI_sim = 0,
                           BAI_mes = 0,
                           WVI_sim = 0,
                           WVI_mes = 0,
                           .rows = nStandYearTot)
  
  index_line = 1
  nIter = nStandYearTot
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  for(i in 1:length(code_site_list)){ # i = 1
    a_code_site = code_site_list[i]
    
    treeYearTable_codesite = subset(treeYearTable, code_site == a_code_site)
    years = unique(treeYearTable_codesite$year)
    for(ayear in years){
      
      treeYearTable_codesite_year = subset(treeYearTable_codesite, year == ayear)
      
      standYearTable[index_line, 1] = a_code_site
      standYearTable[index_line, 2] = ayear
      standYearTable[index_line, 3] = unique(treeYearTable_codesite_year$site)
      standYearTable[index_line, 4] = unique(treeYearTable_codesite_year$composition)
      standYearTable[index_line, 5] = unique(treeYearTable_codesite_year$triplet)
      standYearTable[index_line, 6] = sum(treeYearTable_codesite_year$GPPabs_sim)
      standYearTable[index_line, 7] = sum(treeYearTable_codesite_year$NPPabs_sim)
      standYearTable[index_line, 8] = sum(treeYearTable_codesite_year$RautoAbs_sim)
      standYearTable[index_line, 9] = sum(treeYearTable_codesite_year$BAI_sim)
      standYearTable[index_line, 10] = sum(treeYearTable_codesite_year$BAI_mes)
      standYearTable[index_line, 11] = sum(treeYearTable_codesite_year$WVI_sim)
      standYearTable[index_line, 12] = sum(treeYearTable_codesite_year$WVI_mes)
      index_line = index_line + 1
      
      setTxtProgressBar(pb, index_line) # change progress bar
    }
  }
  close(pb)
  
  rm(treeYearTable_codesite, treeYearTable_codesite_year, years, a_code_site)
  
  return(standYearTable)
}



# Based on a stand scale simulation list, makes a table with one line per stand-year
# 
# meanBAI*nTree/m2 * standArea for each stand-year
makeStandYearTableFromStandScale = function(standSimulist, standYearTable_GMAP){
  code_site_list = names(standSimulist)
  
  nYears = length(unique(standSimulist[[1]]$yearlyResults$year))
  nStandYearTot = length(code_site_list) * nYears
  standYearTable = tibble( code_site = "",
                            year = 0,
                            site = "",
                            composition = "",
                            triplet = "",
                            vegAbso_MJm2 = 0,
                            vegAbsorbance = 0,
                            vegAbsoPAR_MJm2 = 0,
                            vegAbsorbancePAR = 0,
                            GPPabs_sim = 0, # gC/yr
                            NPPabs_sim = 0, # gC/yr
                            RautoAbs_sim = 0, # gC/yr
                            BAI_sim = 0,# BAI from mm2 to cm2
                            WVI_sim = 0, # WVI in unit as in yearlyResults (normally, m3)
                            REWmin = 0, # %
                            RU_level_min = 0,
                            RU_shortage_max = 0,
                            transpiration = 0,
                            BAI_mes = 0,
                            WVI_mes = 0,
                           .rows = nStandYearTot)
  index_line = 1
  
  nIter = nStandYearTot
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  
  for(i in 1:length(code_site_list)){ # i = 1
    a_code_site = code_site_list[i]
    
    codesiteYearlyResults = standSimulist[[a_code_site]]$yearlyResults
    codesiteYearlyResults$year = round(codesiteYearlyResults$year)
    codesiteYearTable_GMAP = subset(standYearTable_GMAP, code_site == a_code_site)
    
    years = unique(codesiteYearlyResults$year)
    
    standArea_m2 = unique(codesiteYearTable_GMAP$standArea_m2)
    
    # if WVI is not in standSimulist variable, set it to -1
    if(!"WVI" %in% names(codesiteYearlyResults)){
      codesiteYearlyResults$WVI = -1
    }
      
    for(ayear in years){
      
      
      codesiteOneYearResults = subset(codesiteYearlyResults, year == ayear)
      codesiteOneYearTable_GMAP = subset(codesiteYearTable_GMAP, year == ayear)
      
      # hauteurMoyenne = codesiteOneYearResults$height
      
      nTreePerStand = codesiteOneYearResults$Nha * (standArea_m2 / 10000)
      
      standYearTable[index_line, 1] = a_code_site
      standYearTable[index_line, 2] = ayear
      standYearTable[index_line, 3] = getSiteFromCodeSite(a_code_site)
      standYearTable[index_line, 4] = getCompositionFromCodeSite(a_code_site)
      standYearTable[index_line, 5] = getTripletFromCodeSite(a_code_site)
      standYearTable[index_line, 6] = codesiteOneYearResults$veg_yearlyMJm2
      standYearTable[index_line, 7] = codesiteOneYearResults$veg_yearlyMJm2 / codesiteOneYearResults$incident_yearlyMJm2
      standYearTable[index_line, 8] = codesiteOneYearResults$vegPAR_yearlyMJm2
      standYearTable[index_line, 9] = codesiteOneYearResults$vegPAR_yearlyMJm2 / codesiteOneYearResults$incidentPAR_yearlyMJm2
      standYearTable[index_line, 10] = codesiteOneYearResults$GPP * standArea_m2
      standYearTable[index_line, 11] = codesiteOneYearResults$NPP * standArea_m2
      standYearTable[index_line, 12] = codesiteOneYearResults$Rauto * standArea_m2
      standYearTable[index_line, 13] = codesiteOneYearResults$BAI * nTreePerStand * 1e-2 # BAI from mm2 to cm2
      standYearTable[index_line, 14] = codesiteOneYearResults$WVI * nTreePerStand
      standYearTable[index_line, 15] = codesiteOneYearResults$REWmin
      standYearTable[index_line, 16] = codesiteOneYearResults$RU_level_min
      standYearTable[index_line, 17] = codesiteOneYearResults$RU_shortage_max
      standYearTable[index_line, 18] = codesiteOneYearResults$TR
      
      standYearTable[index_line, 19] = codesiteOneYearTable_GMAP$BAI_mes
      standYearTable[index_line, 20] = codesiteOneYearTable_GMAP$WVI_mes
      
      
      index_line = index_line + 1
      setTxtProgressBar(pb, index_line) # change progress bar
    } # /end loop years
  }# /end loop sites
  close(pb)
  
  return(standYearTable)
}






# Based on a stand scale simulation list, makes a table with one line per stand-year
# yvarList is the list of yearlyResults variabels to store
# This method is not dependent on a GMAP context
makeStandYearTableFromStandScale_universal = function(standSimulist, yvarList){
  code_site_list = names(standSimulist)
  
  nYears = length(unique(standSimulist[[1]]$yearlyResults$year))
  nStandYearTot = length(code_site_list) * nYears
  
  # create base result table
  standYearTable = tibble(code_site = "", year = 0, standSimulist[[1]]$yearlyResults[1, yvarList], .rows = nStandYearTot)
  
  # set all numerics to 0, "" or FALSE
  standYearTable[sapply(standYearTable[1, ], is.numeric)] = 0
  standYearTable[sapply(standYearTable[1, ], is.character)] = ""
  standYearTable[sapply(standYearTable[1, ], is.logical)] = FALSE
  
  
  index_line = 1
  
  nIter = nStandYearTot
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  
  for(i in 1:length(code_site_list)){ # i = 1
    a_code_site = code_site_list[i]
    
    codesiteYearlyResults = standSimulist[[a_code_site]]$yearlyResults
    codesiteYearlyResults$year = round(codesiteYearlyResults$year)
    
    years = unique(codesiteYearlyResults$year)
    
    for(ayear in years){
      standYearTable[index_line, ]$code_site = a_code_site
      standYearTable[index_line, ]$year = ayear
      codesiteOneYearResults = subset(codesiteYearlyResults, year == ayear)
      standYearTable[index_line, yvarList] = codesiteOneYearResults[yvarList]
      
      index_line = index_line + 1
      setTxtProgressBar(pb, index_line) # change progress bar
    } # /end loop years
  }# /end loop sites
  close(pb)
  
  return(standYearTable)
}







# Create a table with stand-year entries and simulated BAI from CASTANEA as column
# For multispecific simulations, BAI is the aggregate of two CASTANEA monospecifics simulations
makeCASTANEAstandYearTable = function(simulationList, simulationList_preli){
  standYearCASTANEATable = NULL
  
  # for all stand
  nIter = length(names(simulationList_preli))
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  for(i in 1:nIter){ # i = 1
    nameSimuPreli = names(simulationList_preli)[i]
    
    nameSimuPreliSplit = strsplit(nameSimuPreli, "_")[[1]]
    
    # information about simulation stand
    simulationSetPreli = paste0(nameSimuPreliSplit[1:2], collapse = "_")
    simulationSet = paste0(nameSimuPreliSplit[2], collapse = "_")
    
    nameSimu = paste0(nameSimuPreliSplit[2:length(nameSimuPreliSplit)], collapse = "_")
    code_site = paste0(nameSimuPreliSplit[-c(1,2)], collapse = "_")
    site = nameSimuPreliSplit[3]
    composition = nameSimuPreliSplit[5]
    triplet = paste(nameSimuPreliSplit[3], nameSimuPreliSplit[4], nameSimuPreliSplit[6], sep = "_")
    
    standArea_m2 = simulationList_preli[[i]]$inventory$standArea_m2
    
    yearlyResultsSimu = simulationList[[nameSimu]]$yearlyResults
    yearlyResultsSimuPreli = simulationList_preli[[nameSimuPreli]]$yearlyResults
    
    speciesList = unique(yearlyResultsSimu$species)
    if( TRUE %in% (unique(sort(yearlyResultsSimu$species)) != unique(sort(yearlyResultsSimuPreli$species))) ){
      stop("Not the same species in preliminary and simulation")
    }
    
    if( TRUE %in% (unique(sort(yearlyResultsSimu$year)) != unique(sort(yearlyResultsSimuPreli$year))) ){
      stop("Not the same years in preliminary and simulation")
    }
    
    # for all years
    for(ayear in unique(yearlyResultsSimuPreli$year)){
      yearlyResultsSimuPreliYear = subset(yearlyResultsSimuPreli, year == ayear)
      yearlyResultsSimuYear = subset(yearlyResultsSimu, year == ayear)
      
      GPP_aggregated = 0
      BAI_aggregated = 0
      BAI_aggregated2 = 0
      sumTreesBasalArea = sum(pi * (yearlyResultsSimuYear$dbh/2)**2) 
      
      for(aSpeciesId in yearlyResultsSimuPreliYear$species){
        yearlyResultsSimuYearSpecies = subset(yearlyResultsSimuYear, species == aSpeciesId)
        ratioSpecies = sum(pi * (yearlyResultsSimuYearSpecies$dbh/2)**2) / sumTreesBasalArea 
        
        yearlyResultsSimuPreliYearSpecies = subset(yearlyResultsSimuPreliYear, species == aSpeciesId)
        GPP_aggregated = GPP_aggregated + yearlyResultsSimuPreliYearSpecies$GPP * ratioSpecies
        
        # nTreeSpecies_Plot est le nombre d'arbre de l'espèce considéré comptés sur la placette (meilleure methode)
        # nTreeSpecies_est est le nombre d'arbre de l'espèce considéré recalculé à partir du Nha de la simu preli et du ratio species
        nTreeSpecies_Plot = length(unique(yearlyResultsSimuYearSpecies$idFmCell))
        species_stand_area_m2 = standArea_m2 * ratioSpecies
        nTreeSpecies_est = yearlyResultsSimuPreliYearSpecies$Nha / 10000 * species_stand_area_m2
        
        BAI_aggregated = BAI_aggregated + yearlyResultsSimuPreliYearSpecies$BAI * nTreeSpecies_Plot
        BAI_aggregated2 = BAI_aggregated2 + yearlyResultsSimuPreliYearSpecies$BAI * nTreeSpecies_est
      }
      
      
      table_line = tibble( code_site = code_site,
                           year = ayear,
                           GPPabs_simCAST = GPP_aggregated * standArea_m2,
                           BAI_simCAST = BAI_aggregated * 1e-2,
                           BAI_simCAST2 = BAI_aggregated2 * 1e-2 ) # BAI from mm2 to cm2
      
      if(is.null(standYearCASTANEATable)){
        standYearCASTANEATable = table_line
      }else{
        standYearCASTANEATable = rbind(standYearCASTANEATable, table_line)
      }
      
    }
    setTxtProgressBar(pb, i) # change progress bar
  }
  rm(yearlyResultsSimu, nameSimu, nameSimuPreli, yearlyResultsSimuPreli, yearlyResultsSimuPreliYear, yearlyResultsSimuPreliYearSpecies, 
     yearlyResultsSimuYear, yearlyResultsSimuYearSpecies, table_line, GPP_aggregated, BAI_aggregated, BAI_aggregated2, 
     nTreeSpecies_est, nTreeSpecies_Plot, ratioSpecies)
  close(pb)
  
  # NULL BAI
  # standYearCASTANEATable = standYearCASTANEATable[ standYearCASTANEATable$code_site != "bg_haut_ph_1a", ]
  
  return(standYearCASTANEATable)
  
}








# from a treeYearTable (that has information for every tree and every year), construct a table that has an average information for every tree for a set of periods
# Periodes : 1996-2002, 2003-2006, 2007-2013, 1996-2013
makeTreePeriodTable = function(treeYearTable, periodList){
  
  treeGlobalIdList = unique(treeYearTable$treeGlobalId)
  
  nTreePeriodTot = length(treeGlobalIdList) * length(periodList)
  
  treePeriodTable = tibble( treeGlobalId = "",
                            code_site = "",
                            site = "",
                            triplet = "",
                            composition = "",
                            period = "",
                            
                            GPPy_m2_sim = 0,
                            GPPy_abs_sim = 0,
                            NPPy_abs_sim = 0,
                            Rautoy_abs_sim = 0,
                            
                            BAIy_sim = 0,
                            BAIy_mes = 0,
                            WVIy_sim = 0,
                            WVIy_mes = 0,
                            species = "",
                            dbh = 0,
                            hauteur = 0,
                            hcb = 0,
                            span = 0,
                            .rows = nTreePeriodTot)
  
  index_line = 1
  nIter = length(treeGlobalIdList)
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  for(i in 1:nIter){ # i = 1
    aTreeGlobalId  = treeGlobalIdList[i]
    
    treeYearTable_tree = subset(treeYearTable, treeGlobalId == aTreeGlobalId)
    for(j in 1:length(periodList)){
      period = periodList[[j]]
      nYearsPeriod = period[2] - period[1] + 1
      treeYearTable_tree_periodSubset = subset(treeYearTable_tree, year >= period[1] & year <= period[2])
      
      treePeriodTable[index_line, 1] = aTreeGlobalId
      treePeriodTable[index_line, 2] = unique(treeYearTable_tree$code_site)
      treePeriodTable[index_line, 3] = unique(treeYearTable_tree$site)
      treePeriodTable[index_line, 4] = unique(treeYearTable_tree$triplet)
      treePeriodTable[index_line, 5] = unique(treeYearTable_tree$composition)
      treePeriodTable[index_line, 6] = paste0(period, collapse = "_")
      
      treePeriodTable[index_line, 7] = sum(treeYearTable_tree_periodSubset$GPPm2_sim) / nYearsPeriod
      treePeriodTable[index_line, 8] = sum(treeYearTable_tree_periodSubset$GPPabs_sim) / nYearsPeriod
      treePeriodTable[index_line, 9] = sum(treeYearTable_tree_periodSubset$NPPabs_sim) / nYearsPeriod
      treePeriodTable[index_line, 10] = sum(treeYearTable_tree_periodSubset$RautoAbs_sim) / nYearsPeriod
      
      treePeriodTable[index_line, 11] = sum(treeYearTable_tree_periodSubset$BAI_sim) / nYearsPeriod
      treePeriodTable[index_line, 12] = sum(treeYearTable_tree_periodSubset$BAI_mes) / nYearsPeriod
      treePeriodTable[index_line, 13] = sum(treeYearTable_tree_periodSubset$WVI_sim) / nYearsPeriod
      treePeriodTable[index_line, 14] = sum(treeYearTable_tree_periodSubset$WVI_mes) / nYearsPeriod
      treePeriodTable[index_line, 15] = unique(treeYearTable_tree$species)
      treePeriodTable[index_line, 16] = unique(treeYearTable_tree$dbh)
      treePeriodTable[index_line, 17] = unique(treeYearTable_tree$hauteur)
      treePeriodTable[index_line, 18] = unique(treeYearTable_tree$hcb)
      treePeriodTable[index_line, 19] = unique(treeYearTable_tree$span)
      
      index_line = index_line + 1
      
    } #end loop on period
    
    
    setTxtProgressBar(pb, i) # change progress bar
    
  }
  
  close(pb)
  rm(treeYearTable_tree)
  return(treePeriodTable)
}




# from a standYearTable (that has information for every stand and every year), construct a table that has an average information for every stand for a set of periods
makeStandPeriodTable = function(standYearTable, periodList){
  standPeriodTable = NULL
  
  code_site_list = unique(standYearTable$code_site)
  nIter = length(code_site_list)
  
  pb <- txtProgressBar(min = 0, max = nIter, style = 3, width = 50, char = "=") 
  for(i in 1:nIter){ # i = 1
    a_code_site  = code_site_list[i]
    
    standYearTable_stand = subset(standYearTable, code_site == a_code_site)
    
    for(j in 1:length(periodList)){
      period = periodList[[j]]
      nYearsPeriod = period[2] - period[1] + 1
      standYearTable_stand_subset = subset(standYearTable_stand, year >= period[1] & year <= period[2])
      
      table_line = tibble( code_site = a_code_site,
                           site = unique(standYearTable_stand$site),
                           triplet = unique(standYearTable_stand$triplet),
                           composition = unique(standYearTable_stand$composition),
                           period = paste0(period, collapse = "_"),
                           
                           vegAbso_MJm2Y = sum(standYearTable_stand_subset$vegAbso_MJm2) / nYearsPeriod,
                           vegAbsorbance = sum(standYearTable_stand_subset$vegAbsorbance) / nYearsPeriod,
                           vegAbsoPAR_MJm2Y = sum(standYearTable_stand_subset$vegAbsoPAR_MJm2) / nYearsPeriod,
                           vegAbsorbancePAR = sum(standYearTable_stand_subset$vegAbsorbancePAR) / nYearsPeriod,
                           
                           GPPy_abs_sim = sum(standYearTable_stand_subset$GPPabs_sim) / nYearsPeriod,
                           NPPy_abs_sim = sum(standYearTable_stand_subset$NPPabs_sim) / nYearsPeriod,
                           Rautoy_abs_sim = sum(standYearTable_stand_subset$RautoAbs_sim) / nYearsPeriod,
                           # GPPy_abs_simCAST = sum(standYearTable_stand_subset$GPPabs_simCAST) / nYearsPeriod,
                           BAIy_sim = sum(standYearTable_stand_subset$BAI_sim) / nYearsPeriod,
                           BAIy_mes = sum(standYearTable_stand_subset$BAI_mes) / nYearsPeriod,
                           WVIy_sim = sum(standYearTable_stand_subset$WVI_sim) / nYearsPeriod,
                           WVIy_mes = sum(standYearTable_stand_subset$WVI_mes) / nYearsPeriod,
                           REWmin = sum(standYearTable_stand_subset$REWmin) / nYearsPeriod,
                           RU_level_min = sum(standYearTable_stand_subset$RU_level_min) / nYearsPeriod,
                           transpiration = sum(standYearTable_stand_subset$transpiration) / nYearsPeriod,
                           RU_shortage_max = sum(standYearTable_stand_subset$RU_shortage_max) / nYearsPeriod,
                           # BAIy_simCAST = sum(standYearTable_stand_subset$BAI_simCAST) / nYearsPeriod,
                           # BAIy_simCAST2 = sum(standYearTable_stand_subset$BAI_simCAST2) / nYearsPeriod
                           )
      
      if(is.null(standPeriodTable)){
        standPeriodTable = table_line
      }else{
        standPeriodTable = rbind(standPeriodTable, table_line)
      }
    }
    
    setTxtProgressBar(pb, i) # change progress bar
    
  }
  
  close(pb)
  rm(standYearTable_stand)
  return(standPeriodTable)
}




# from a stand-period table, give a set of coefficients of correlation between two variable of a table
getComparisonCoefficient = function(table, nameVar1, nameVar2, coefficientList = c("correlation", "r2", "1-r2", "RMSE", "MAPE")){

  res = NULL
  if("correlation" %in% coefficientList){
    res = tibble(test = "correlation", 
                 PDGvsMES = cor(table[[nameVar1]], table[[nameVar2]]))
  }
  
  if("r2" %in% coefficientList){
    res = rbind(res, tibble(test = "r2", 
                            PDGvsMES = cor(table[[nameVar1]], table[[nameVar2]])**2))
  }
  
  if("1-r2" %in% coefficientList){
    res = rbind(res, tibble(test = "1-r2", 
                            PDGvsMES = 1 - cor(table[[nameVar1]], table[[nameVar2]])**2))
  }
  
  if("RMSE" %in% coefficientList){
    res = rbind(res, tibble(test = "RMSE", 
                            PDGvsMES = rmse(table[[nameVar1]], table[[nameVar2]])))
  }
  
  if("MAPE" %in% coefficientList){
    res = rbind(res, tibble(test = "MAPE", 
                            PDGvsMES = mape(table[[nameVar1]], table[[nameVar2]])))
  }
  
  return(res)
}


# from a stand-period table, give a set of coefficients of correlation between two variable of a table
# coefficients are computed for subset of the table based on compositions and sites
getComparisonCoefficientPerSiteAndComposition = function(standPeriodTable, nameVar1, nameVar2,  coefficientList = c("correlation", "r2", "1-r2", "RMSE", "MAPE")){
  res = list()
  
  standPeriodTable_m = subset(standPeriodTable, composition == "m")
  standPeriodTable_ph = subset(standPeriodTable, composition == "ph")
  standPeriodTable_sp = subset(standPeriodTable, composition == "sp")
  
  standPeriodTable_mono = subset(standPeriodTable, composition %in% c("ph", "sp"))
  
  standPeriodTable_bg = subset(standPeriodTable, site == "bg")
  standPeriodTable_vl = subset(standPeriodTable, site == "vl")
  standPeriodTable_vtx = subset(standPeriodTable, site == "vtx")
  
  
  standPeriodTable_bg_m = subset(standPeriodTable, site == "bg" & composition == "m")
  standPeriodTable_vl_m = subset(standPeriodTable, site == "vl" & composition == "m")
  standPeriodTable_vtx_m = subset(standPeriodTable, site == "vtx" & composition == "m")
  
  standPeriodTable_bg_ph = subset(standPeriodTable, site == "bg" & composition == "ph")
  standPeriodTable_vl_ph = subset(standPeriodTable, site == "vl" & composition == "ph")
  standPeriodTable_vtx_ph = subset(standPeriodTable, site == "vtx" & composition == "ph")
  
  standPeriodTable_bg_sp = subset(standPeriodTable, site == "bg" & composition == "sp")
  standPeriodTable_vl_sp = subset(standPeriodTable, site == "vl" & composition == "sp")
  standPeriodTable_vtx_sp = subset(standPeriodTable, site == "vtx" & composition == "sp")
  
  
  res$global = getComparisonCoefficient(standPeriodTable, nameVar1, nameVar2, coefficientList)
  res$comp_m = getComparisonCoefficient(standPeriodTable_m, nameVar1, nameVar2, coefficientList)
  res$comp_ph = getComparisonCoefficient(standPeriodTable_ph, nameVar1, nameVar2, coefficientList)
  res$comp_sp = getComparisonCoefficient(standPeriodTable_sp, nameVar1, nameVar2, coefficientList)
  
  res$comp_mono = getComparisonCoefficient(standPeriodTable_mono, nameVar1, nameVar2, coefficientList)
  
  res$site_bg = getComparisonCoefficient(standPeriodTable_bg, nameVar1, nameVar2, coefficientList)
  res$site_vl = getComparisonCoefficient(standPeriodTable_vl, nameVar1, nameVar2, coefficientList)
  res$site_vtx = getComparisonCoefficient(standPeriodTable_vtx, nameVar1, nameVar2, coefficientList)
  
  res$bg_m = getComparisonCoefficient(standPeriodTable_bg_m, nameVar1, nameVar2, coefficientList)
  res$vl_m = getComparisonCoefficient(standPeriodTable_vl_m, nameVar1, nameVar2, coefficientList)
  res$vtx_m = getComparisonCoefficient(standPeriodTable_vtx_m, nameVar1, nameVar2, coefficientList)
  
  res$bg_ph = getComparisonCoefficient(standPeriodTable_bg_ph, nameVar1, nameVar2, coefficientList)
  res$vl_ph = getComparisonCoefficient(standPeriodTable_vl_ph, nameVar1, nameVar2, coefficientList)
  res$vtx_ph = getComparisonCoefficient(standPeriodTable_vtx_ph, nameVar1, nameVar2, coefficientList)
  
  res$bg_sp = getComparisonCoefficient(standPeriodTable_bg_sp, nameVar1, nameVar2, coefficientList)
  res$vl_sp = getComparisonCoefficient(standPeriodTable_vl_sp, nameVar1, nameVar2, coefficientList)
  res$vtx_sp = getComparisonCoefficient(standPeriodTable_vtx_sp, nameVar1, nameVar2, coefficientList)
  
  return(res)
}



# add WVI (wood volume increment, in m3) in yearlyResults
# based on simulated BAI (Basal Area Increment) per tree and year and hauteur of GMAP data per tree
# set averagePerSpecies to TRUE if it should use average height per species per stand
# 6.11.2023
addWVIfromBAIYearlyResultsAndGMAPheight = function(simuList, dGMAP_horsprod, averagePerSpecies = FALSE){
  
  # pour chaque site
  for(aCodeSite in names(simuList)){
    aCodeSite_clean = getCodeSiteFromSimulationName(aCodeSite)
    simuList[[aCodeSite]]$yearlyResults$WVI = 0
    
    if(averagePerSpecies){
      hauteurMoyenneSp = list()
      
      hauteurMoyenneSp[["sapin"]] = getAverageHeightBasedOnIndividualTreeVolumes(dGMAP_horsprod, aCodeSite_clean, speciesName = "sapin")
      hauteurMoyenneSp[["hetre"]] = getAverageHeightBasedOnIndividualTreeVolumes(dGMAP_horsprod, aCodeSite_clean, speciesName = "hetre")
    }
    
    aSimuYearlyResults = simuList[[aCodeSite]]$yearlyResults
    
    # pour chaque arbre
    for(i in 1:length(unique(aSimuYearlyResults$idFmCell))){
      
      aIdFmCell = unique(aSimuYearlyResults$idFmCell)[i]
      
      if(!averagePerSpecies){
        aTreeGlobalId = paste0(aCodeSite_clean, "_", aIdFmCell)
        
        if(length(aIdFmCell) <= 0){
          browser()
        }
        if(aIdFmCell > 1000000){
          simuList[[aCodeSite]]$yearlyResults$WVI[simuList[[aCodeSite]]$yearlyResults$idFmCell == aIdFmCell] = -1
          next
        }
        
        dGMAP_tree = subset(dGMAP_horsprod, treeGlobalId == aTreeGlobalId)
        hauteurTree = dGMAP_tree$htot
      }else{
        if(unique(subset(aSimuYearlyResults, idFmCell == aIdFmCell)$speciesName) == "Abies alba"){
          hauteurTree = hauteurMoyenneSp[["sapin"]]
        }
        if(unique(subset(aSimuYearlyResults, idFmCell == aIdFmCell)$speciesName) == "Fagus sylvatica"){
          hauteurTree = hauteurMoyenneSp[["hetre"]]
        }
      }
      
      # WVI Wood Volume Index in m3
      # BAI in mm2 to m2
      # height in m 
      simuList[[aCodeSite]]$yearlyResults$WVI[simuList[[aCodeSite]]$yearlyResults$idFmCell == aIdFmCell] =
        simuList[[aCodeSite]]$yearlyResults$BAI[simuList[[aCodeSite]]$yearlyResults$idFmCell == aIdFmCell] * (0.001)**2 * hauteurTree 
        
    }
  }
  return(simuList)
}




# add WVI (wood volume increment, in m3) in yearlyResults for stand scale simulation
# based on simulated mean BAI (Basal Area Increment) per tree per stand and year and hauteur of GMAP data per tree
# 23.11.2023
addWVIfromBAIYearlyResultsAndGMAPheight_standScale = function(simuList, dGMAP_horsprod){
  
  # pour chaque site
  for(aCodeSite in names(simuList)){
    aCodeSite_clean = getCodeSiteFromSimulationName(aCodeSite)
    simuList[[aCodeSite]]$yearlyResults$WVI = 0
    
    aSimuYearlyResults = simuList[[aCodeSite]]$yearlyResults
    
    if(length(unique(aSimuYearlyResults$species)) > 1){
      stop("More thant one speciesi in simuList[[aCodeSite]]")
    }
    
    
    speciesName = ifelse(unique(aSimuYearlyResults$species) == 1, "sapin", 
                         ifelse(unique(aSimuYearlyResults$species) == 3, "hetre", stop("species is not hetre nor sapin.")))
    
    hauteurTree = getAverageHeightBasedOnIndividualTreeVolumes(dGMAP_horsprod, aCodeSite_clean, speciesName = speciesName)
    
    # WVI Wood Volume Index in m3
    # BAI from mm2 to m2
    # height in m
    simuList[[aCodeSite]]$yearlyResults$WVI =
      simuList[[aCodeSite]]$yearlyResults$BAI * 0.001**2 * hauteurTree 
     
  }
  return(simuList)
}




# Compute the average height of a stand (in cm) based on its tree volume sum and basal area
# Need GMAP data
getAverageHeightBasedOnIndividualTreeVolumes = function(dGMAP_horsprod, a_code_site, speciesName = NULL){
  
  if(!a_code_site %in% dGMAP_horsprod$code_site){
    stop("a_code_site is not found in dGMAP_horsprod.")
  }
  
  dGMAP_horsprod_code_site = subset(dGMAP_horsprod, code_site == a_code_site)
  
  if(is.null(speciesName)){
    dGMAP_horsprod_code_site_species = dGMAP_horsprod_code_site
  }else{
    if(speciesName == "sapin"){
      dGMAP_horsprod_code_site_species = subset(dGMAP_horsprod_code_site, essence == "sapin")
    }else if(speciesName == "hetre"){
      dGMAP_horsprod_code_site_species = subset(dGMAP_horsprod_code_site, essence == "hetre")
    }else{
      stop("speciesName must be hetre or sapin.")
    }
  }
  
  volumeTotal = 0
  basalAreaTotal = 0
  
  # pour chaque arbre, ajouter le volume et la surface terriere
  for(i in 1:dim(dGMAP_horsprod_code_site_species)[1]){
    dbh = dGMAP_horsprod_code_site_species[i, ]$circonference / pi # in cm
    hauteur = dGMAP_horsprod_code_site_species[i, ]$htot * 100 # in cm
    basalAreaTree = pi * (dbh/2) **2  # in cm2 
    volumeTree = basalAreaTree * hauteur # in cm3
    
    volumeTotal = volumeTotal + volumeTree
    basalAreaTotal = basalAreaTotal + basalAreaTree
  }
  
  averagedHeight = volumeTotal / basalAreaTotal * 0.01 # in m
  
  return(averagedHeight)
}



# ..............................................................................
# GMAP metainfo ----------------------------------------------------------------



# extract GMAP code_site from the name of the simulation
getCodeSiteFromSimulationName = function(nameSimu){
  
  if(length(nameSimu) > 1){
    return(vapply(nameSimu, FUN = getCodeSiteFromSimulationName, FUN.VALUE = "", USE.NAMES = F))
  }
  
  nameSimuSplit = strsplit(nameSimu, "_")[[1]]
  firstIndex = which(vapply(nameSimuSplit, FUN = function(x) x %in% c("vl", "bg", "vtx"), FUN.VALUE = T, USE.NAMES = F))
  hasFirstCharacterAsNumeric = sapply(nameSimuSplit, FUN = function(x) !is.na(suppressWarnings(as.numeric(strsplit(x, split = "")[[1]][1]))), USE.NAMES = F)
  lastIndex = min(which(hasFirstCharacterAsNumeric))
  code_site = paste0(nameSimuSplit[firstIndex:lastIndex], collapse = "_")
  
  return(code_site)
}

# extract GMAP site from code_site
getSiteFromCodeSite = function(code_site){
  for(site in c("bg", "vl", "vtx")){
    sitePattern = paste0(site, "_")
    if(grepl(code_site, pattern = sitePattern)){
      return(site)
    }
  }
}

# extract GMAP composition from code_site
getCompositionFromCodeSite = function(code_site){
  for(composition in c("m", "ph", "sp")){
    compositionPattern = paste0("_", composition, "_")
    if(grepl(code_site, pattern = compositionPattern)){
      return(composition)
    }
  }
}

# extract GMAP subsite from code_site
getSubsiteFromCodeSite = function(code_site){
  for(subsite in c("bas", "inter", "haut")){
    subsitePattern = paste0("_", subsite, "_")
    if(grepl(code_site, pattern = subsitePattern)){
      return(subsite)
    }
  }
}

# extract GMAP subsite number from code_site
getNumberFromCodeSite = function(code_site){
  splitCodeSite = strsplit(code_site, split = "_")[[1]]
  nElements = length(splitCodeSite)
  return(splitCodeSite[nElements])
}


# extract GMAP triplet from code_site
getTripletFromCodeSite = function(code_site){
  site = getSiteFromCodeSite(code_site)
  subsite = getSubsiteFromCodeSite(code_site)
  number = getNumberFromCodeSite(code_site)
  return(paste(site, subsite, number, sep = "_"))
}




# ONF inventories meta-info ----

# get RU in integer based on codesite
getRUfromCode_site = function(a_code_site){
  
  if(length(a_code_site) > 1){
    return( vapply(a_code_site, FUN = getRUfromCode_site, FUN.VALUE = 0, USE.NAMES = F) )
  }
  
  splitlist = strsplit(a_code_site, split = "_")[[1]]
  RUstr = splitlist[grepl(splitlist, pattern = "RU")]
  
  # only the number
  RUstr = as.integer(substr(RUstr, 3 , nchar(RUstr)))
  
  return(RUstr)
}



# get tree class size in string based on codesite
getClassSizefromCode_site = function(a_code_site){
  
  if(length(a_code_site) > 1){
    return( vapply(a_code_site, FUN = getClassSizefromCode_site, FUN.VALUE = "", USE.NAMES = F) )
  }
  
  splitlist = strsplit(a_code_site, split = "_")[[1]]
  treeSizeStr = splitlist[grepl(splitlist, pattern = "PB|BM|GB")]
  # treeSizeStr = substr(treeSizeStr, 2, nchar(treeSizeStr))
  
  
  return(treeSizeStr)
}

# get composition in string based on codesite
getCompositionfromCode_site = function(a_code_site){
  
  if(length(a_code_site) > 1){
    return( vapply(a_code_site, FUN = getCompositionfromCode_site, FUN.VALUE = "", USE.NAMES = F) )
  }
  
  splitlist = strsplit(a_code_site, split = "_")[[1]]
  compositionStr = splitlist[grepl(splitlist, pattern = "het|sap|HET|SAP")]
  compositionStr = substr(compositionStr, 2, nchar(compositionStr))
  
  
  return(compositionStr)
}

# get composition index based on codesite
getCompositionIndexfromCode_site = function(a_code_site){
  
  if(length(a_code_site) > 1){
    return( vapply(a_code_site, FUN = getCompositionIndexfromCode_site, FUN.VALUE = 0, USE.NAMES = F) )
  }
  
  splitlist = strsplit(a_code_site, split = "_")[[1]]
  compositionIndexStr = splitlist[grepl(splitlist, pattern = "het|sap|HET|SAP")]
  compositionIndexStr = substr(compositionIndexStr, 1, 1)
  compositionIndexStr = as.integer(compositionIndexStr)
  
  return(compositionIndexStr)
}



# get ntree in integer based on codesite
getNTreefromCode_site = function(a_code_site){
  
  if(length(a_code_site) > 1){
    return( vapply(a_code_site, FUN = getNTreefromCode_site, FUN.VALUE = 0, USE.NAMES = F) )
  }
  
  splitlist = strsplit(a_code_site, split = "__")[[1]]
  nTreeStr = splitlist[grepl(splitlist, pattern = "Narbres")]
  nTreeStr = substr(nTreeStr, 9, nchar(nTreeStr))
  nTreeStr = as.integer(nTreeStr)
  return(nTreeStr)
}

# get nha in integer based on codesite
getNhafromCode_site = function(a_code_site){
  
  if(length(a_code_site) > 1){
    return( vapply(a_code_site, FUN = getNhafromCode_site, FUN.VALUE = 0, USE.NAMES = F) )
  }
  
  splitlist = strsplit(a_code_site, split = "__")[[1]]
  nhaStr = splitlist[grepl(splitlist, pattern = "Nha")]
  nhaStr = substr(nhaStr, 5, nchar(nhaStr))
  nhaStr = as.integer(nhaStr)
  return(nhaStr)
}