# This work © 2024 by Camille Rouet is licensed under [CC BY-NC 4.0](http://creativecommons.org/licenses/by-nc/4.0/).

# Read a the PDG inventory option file
# It is returned as a list of table
readInventoryOptionFile = function(filePath, delimTable = ";"){
  
  # a. READ ALL LINES
  lines = read_lines(file = filePath)
  
  
  # b. EXTRACT EMPTY AND COMMENT LINES
  lines_empty = substr(lines, 1, 1) == ";"
  lines_empty_index = which(lines_empty)
  
  comment_lines = substr(lines, 1, 1) == "#"
  comment_lines_index = which(comment_lines)
  
  
  
  # c.1. PREPARE TABLE EXTRACTION
  # extract table entry lines
  lines_table_entries = substr(lines,1,3) == ">T>"
  
  # make a table with table entry name and entry line index
  table_entries_names = gsub(lines[lines_table_entries], pattern = ";| |>T>", replacement = "")
  table_entries_index_vec = which(lines_table_entries)
  table_entries_index = as_tibble(t(data.frame( table_entries_index_vec ))) ; colnames(table_entries_index) = table_entries_names
  
  # Find the column names line
  tables_entries_column_index_vec = table_entries_index_vec+1
  
  # pass comment lines
  while(sum(tables_entries_column_index_vec %in% comment_lines_index) >= 1){
    tables_entries_column_index_vec[tables_entries_column_index_vec %in% comment_lines_index] = tables_entries_column_index_vec[tables_entries_column_index_vec %in% comment_lines_index]+ 1
  }
  table_entries_column_lines = lines[tables_entries_column_index_vec]
  
  count_column = function(entry_column_line){
    entry_columns_split = strsplit(entry_column_line, split = ";")[[1]]
    entry_columns_split_clean = entry_columns_split[entry_columns_split != ""]
    return( length(entry_columns_split_clean))
  }
  
  table_entries_nbColumns_vec = vapply(table_entries_column_lines, FUN = count_column, FUN.VALUE = 0, USE.NAMES = F)
  table_entries_nbColumns = as_tibble(t(data.frame( table_entries_nbColumns_vec ))) ; colnames(table_entries_nbColumns) = table_entries_names
  
  
  
  
  
  # import a table that is an entry of the inventory file
  readTableInInventoryOptionFile = function(entryName, 
                                            filePath_ = filePath, delimTable_ = delimTable, 
                                            lines_ = lines,
                                            table_entries_index_ = table_entries_index, table_entries_nbColumns_ = table_entries_nbColumns, lines_empty_index_ = lines_empty_index){
    
    if(!entryName %in% colnames(table_entries_index_)){
      stop(paste0("entryName ", entryName, " is not in table_entries_index"))
    }
    
    if( sum(lines_empty_index_ > table_entries_index_[[entryName]] ) > 0 ){
      firstEmptyLine = min(lines_empty_index_[lines_empty_index_ > table_entries_index_[[entryName]]])
      tableNbLines = firstEmptyLine - table_entries_index_[[entryName]] - 2
    }else{
      tableNbLines = Inf
    }
    
    # Count the number of comment line to ignore it
    tableLines = lines_[ (table_entries_index_[[entryName]]+1) : (table_entries_index_[[entryName]] + 1 + tableNbLines)]
    tableNbLines = tableNbLines - sum(substr(tableLines, 1,1) == "#")
    
    entryTable = read_delim(filePath_, delim = delimTable_, skip = table_entries_index_[[entryName]], 
                            n_max = tableNbLines, 
                            col_names = T, show_col_types = F,
                            comment = "#")
    entryTable = entryTable[, 1:table_entries_nbColumns_[[entryName]] ]
    return(entryTable)
  }
  
  
  
  
  # c.2. PREPARE LIST EXTRACTION
  # extract list entry lines
  lines_list_entries = substr(lines,1,3) == ">L>"
  
  # make a table with table entry name and entry line index
  list_entries_names = gsub(lines[lines_list_entries], pattern = ";| |>L>", replacement = "")
  list_entries_index_vec = which(lines_list_entries)
  list_entries_index = as_tibble(t(data.frame( list_entries_index_vec ))) ; colnames(list_entries_index) = list_entries_names
  
  # import a table that is an entry of the inventory file
  readListInInventoryOptionFile = function(entryName, 
                                           filePath_ = filePath, delimTable_ = delimTable, 
                                           list_entries_index_ = list_entries_index, lines_empty_index_ = lines_empty_index){
    
    if(!entryName %in% colnames(list_entries_index_)){
      stop(paste0("entryName ", entryName, " is not in list_entries_index_"))
    }
    
    if( sum(lines_empty_index_ > list_entries_index_[[entryName]] ) > 0 ){
      firstEmptyLine = min(lines_empty_index_[lines_empty_index_ > list_entries_index_[[entryName]]])
      listNbLines = firstEmptyLine - list_entries_index_[[entryName]] - 1
    }else{
      listNbLines = Inf
    }
    
    entryList = read_delim(filePath_, delim = delimTable_, skip = list_entries_index_[[entryName]], 
                           n_max = listNbLines,
                           col_names = F,
                           show_col_types = F)
    variableNames = entryList$X1
    entryList = as_tibble(t(entryList$X2))
    colnames(entryList) = variableNames
    
    for(col in colnames(entryList)){
      numericalConvertValue = suppressWarnings(as.numeric(entryList[[col]][1]))
      booleanConvertValue = suppressWarnings(as.logical(entryList[[col]][1]))
      
      if( !is.na(numericalConvertValue) ){
        entryList[[col]] = numericalConvertValue
      }else if( !is.na(booleanConvertValue) ){
        entryList[[col]] = booleanConvertValue
      }
    }
    
    return(entryList)
  }
  
  
  
  
  # c.3. PREPARE VECTOR EXTRACTION
  # extract list entry lines
  lines_vec_entries = substr(lines,1,3) == ">V>"
  
  # make a table with table entry name and entry line index
  vec_entries_names = gsub(lines[lines_vec_entries], pattern = ";| |>V>", replacement = "")
  vec_entries_index_vec = which(lines_vec_entries)
  vec_entries_index = as_tibble(t(data.frame( vec_entries_index_vec ))) ; colnames(vec_entries_index) = vec_entries_names
  
  
  # import a vector that is an entry of the inventory file
  readVectorInInventoryOptionFile = function(entryName, 
                                             filePath_ = filePath, delimTable_ = delimTable, 
                                             lines_ = lines,
                                             vec_entries_index_ = vec_entries_index, lines_empty_index_ = lines_empty_index){
    
    if(!entryName %in% colnames(vec_entries_index_)){
      stop(paste0("entryName ", entryName, " is not in vec_entries_index"))
    }
    
    entryVector = read_delim(filePath_, delim = delimTable_, skip = vec_entries_index_[[entryName]], 
                             n_max = 1, 
                             col_names = F, show_col_types = F,
                             comment = "#")
    entryVector = entryVector[!is.na(entryVector)]
    
    return(entryVector)
  }
  
  
  
  
  # d. READ TABLES AND LINES
  demographic_parameters = readTableInInventoryOptionFile("demographic_parameters")
  genetics_parameters = readTableInInventoryOptionFile("genetics_parameters")
  PDG_simulation_parameters = readListInInventoryOptionFile("PDG_simulation_parameters")
  CASTANEA_simulation_parameters = readListInInventoryOptionFile("CASTANEA_simulation_parameters")
  distance_class_bounds_for_SGS = readVectorInInventoryOptionFile("distance_class_bounds_for_SGS")
  
  
  return( list(demographic_parameters = demographic_parameters, 
               PDG_simulation_parameters = PDG_simulation_parameters, 
               CASTANEA_simulation_parameters = CASTANEA_simulation_parameters,
               genetics_parameters = genetics_parameters,
               distance_class_bounds_for_SGS = distance_class_bounds_for_SGS))
}




# 
readCASTANEASpeciesFile = function(pathtofile){
  lines = read_delim(pathtofile, delim = ";")
  table = suppressWarnings(read_delim(pathtofile, delim = "\t"))
  
  # global parameters
  firstLineParam = which( apply(lines, FUN = function(x) grepl(pattern = "MRN", x), MARGIN = 1) )
  lastLineParam = which( apply(lines, FUN = function(x) grepl(pattern = "Asoil", x), MARGIN = 1) )
  linesParam_tmp = apply(lines[firstLineParam:lastLineParam, ], FUN = function(x) strsplit(gsub(x = x, pattern = '\t', replacement = ""), split = "=")[[1]], MARGIN = 1)
  linesParam = data.frame( t ( as.numeric(linesParam_tmp[2, ])))
  colnames(linesParam) = linesParam_tmp[1, ]
  
  # species table
  beforeFirstLineTable = which( apply(lines, FUN = function(x) grepl(pattern = "speciesName", x), MARGIN = 1) )
  lastLineTable = dim(table)[1]
  table_table = table[(beforeFirstLineTable+1):lastLineTable, ]
  colnames(table_table) = table[beforeFirstLineTable, ]
  
  for(col in colnames(table_table)){
    echantillon_col = suppressWarnings(as.numeric(table_table[[col]][1]))
    if( !is.na( echantillon_col ) ){
      table_table[[col]] = as.numeric(table_table[[col]])
    }
  }
  
  return(list(parameters = linesParam, speciesTable = table_table))
}

# Compute RU (water useful reserve) and add it into the table
computeAndAddRU = function(soil_table, signif = NULL){
  val = (soil_table$wfc - soil_table$wilt) * soil_table$soilHeight * (1 - soil_table$stoneContent) * soil_table$bulk
  
  if(!is.null(signif)){
    val = signif(val, signif)
  }
  
  soil_table$RU = val
  return(soil_table)
}


# compute cell id, col index and line index
computeCellInfos = function(cell_line, nCellSide){
  cells_new = tibble(cell_line, .rows = nCellSide**2)
  IdCell=1:nCellSide**2
  
  clign = rep(1:nCellSide-1, each = nCellSide)
  ccol = rep(1:nCellSide-1, nCellSide)
  cells_new$clign = clign
  cells_new$ccol = ccol
  cells_new$IdCell = IdCell
  return(cells_new)
}


# find cell id and relative positon in cell of a x,y coordinate
# [id, rx, ry]
find_cell_and_rpositon = function(x, y, cells, cellWidth){
  
  # If multi
  if(length(x) > 1){
    if(length(x) != length(y)){
      stop("x and y has different lengths")
    }
    
    
    resTable = tibble("cellId" = 0, "res_x" = 0, "res_y" = 0, .rows = length(x))
    
    for(i in 1:length(x)){
      res = find_cell_and_rpositon(x[i], y[i], cells, cellWidth)
      resTable[i, 1] = res[1]
      resTable[i, 2] = res[2]
      resTable[i, 3] = res[3]
    }
    return(resTable)
  }
  
  nCells = length(unique(cells$IdCell))
  nCells_side = sqrt(nCells)
  sideWidth = nCells_side * cellWidth
  
  if(max(x) > sideWidth){
    stop("x position greater than width of plot has been found")
  }
  if(max(y) > sideWidth){
    stop("y position greater than width of plot has been found")
  }
  
  coll = floor(x / cellWidth)
  lign = floor(y / cellWidth)
  res_x = x - coll * cellWidth
  res_y = cellWidth - (y - lign * cellWidth)
  # res_y = y - lign * cellWidth # OLD (changed by cr-15-12-2023)
  
  
  cellId = cells[cells$ccol == coll & cells$clign == lign, ]$IdCell
  if(length(cellId) == 0){
    stop("Regularization: cell was not found for a tree.")
  }
  return(c(cellId, res_x, res_y))
}



# Compute altitude at a given point (x, y) according to slope and aspect
# referenceAltitude is altitude at (x, y) = (0, 0)
compAltitude = function(x, y, slope_deg, aspect_deg, referenceAltitude){
  degToRad = 2 * pi  / 360
  return(referenceAltitude + (-x * sin( aspect_deg * degToRad) - y * cos(aspect_deg * degToRad) ) * tan(slope_deg * degToRad))
}



# Compute a set of regular spaced altitude classes, knowing extreme X and Y position and slope, ascptc, reference_altitude..
getAltittudeClases = function(xMin, xMax, yMin, yMax, slope_deg, aspect_deg, referenceAltitude, nbAltitudeClasses = 5){
  
  # four extreme point altitudes
  altitudes = c(compAltitude(xMin, yMin, slope_deg, aspect_deg, referenceAltitude),
                compAltitude(xMin, yMax, slope_deg, aspect_deg, referenceAltitude),
                compAltitude(xMax, yMin, slope_deg, aspect_deg, referenceAltitude),
                compAltitude(xMax, yMax, slope_deg, aspect_deg, referenceAltitude))
  
  minimalAltitude = 0.99 * min(altitudes) + 0.01 * max(altitudes)
  maximumAltitude = 0.99 * max(altitudes) + 0.01 * min(altitudes)
  
  altitudeClasses = seq(minimalAltitude, maximumAltitude, length.out = nbAltitudeClasses)
  
  return(altitudeClasses)
}





# Method of allometry Forrester, aLeafAreaDBH and bLeafAreaDBH are from CAsTANEa species file
getLeafAreafromDBH = function(dbh, aLeafAreaDBH, bLeafAreaDBH){
  logLeafArea = aLeafAreaDBH + bLeafAreaDBH * log(dbh)
  leafAreaTree = exp(logLeafArea) # in m2 per tree
  return(leafAreaTree)
}



# Compute the LAI value of a species
# trees is the table with all species trees
# LAI_LIDAR is the measured LAI of the allspecies plot
# standArea is the area of the allspecies plot
computeLAI_specific = function(trees, idSpecies, LAI_LIDAR, standArea, standArea_sp, CASTANEASpeciesFile, demographic_parameters){
  
  CASTANEASpeciesTable = CASTANEASpeciesFile$speciesTable
  
  speciesFullNames = demographic_parameters$speciesFullNames
  idSp = demographic_parameters$idSp
  
  
  # Compute leaf area of trees and species trees based on forrester allometry
  sumLeaf_sp_th = 0
  sumLeaf_th = 0
  for(itree in 1:dim(trees)[1]){
    tree = trees[itree, ]
    
    idSpTree = tree$sp
    aLeafAreaDBH = CASTANEASpeciesTable[CASTANEASpeciesTable$speciesName == speciesFullNames[which(idSp == idSpTree)], ]$aLeafAreaDBH
    bLeafAreaDBH = CASTANEASpeciesTable[CASTANEASpeciesTable$speciesName == speciesFullNames[which(idSp == idSpTree)], ]$bLeafAreaDBH
    
    treeLeafArea = getLeafAreafromDBH(tree$dbh, aLeafAreaDBH, bLeafAreaDBH)
    sumLeaf_th = sumLeaf_th + treeLeafArea
    
    if(tree$sp == idSpecies){
      sumLeaf_sp_th = sumLeaf_sp_th + treeLeafArea
    }
  }
  
  # Compute leaf area as measured by LIDAR
  sumLeaf_real = LAI_LIDAR * standArea
  
  # Evalueate the variation between measured leaf area and allometry-based leaf area
  ratioLeafAreaCorrection = sumLeaf_real / sumLeaf_th
  
  # Apply the variation to allometry-based specific leaf area
  sumLeaf_sp_real = sumLeaf_sp_th * ratioLeafAreaCorrection
  LAI_sp = sumLeaf_sp_real / standArea_sp
  # note cahier : LAI_sp is the same in the original demography-fine plot and the new regularized monospecific plot
  
  return(LAI_sp)
}


# Compute missing value of variableY based on existing values of variableX
# This suppose that every variableX are valid measure and that some values of variableY do exist
# Methods are lm (linear model with slope & intercept), proportionnal_lm (linear model with slope and null intercept) and proportionnal_mean (using a mean proportionnality coefficients)
computeMissingVariableProportionnaly = function(tree_table, variableX, variableY, variableSpecies, variableYisMissing, method = "lm"){
  
  treesToModify = tree_table[variableYisMissing, ]
  trees_valid = tree_table[!variableYisMissing, ]
  
  # in order of apparition of species
  species_list = unique(treesToModify[[variableSpecies]])
  
  # also in order of apparition of species
  intercept_list = c()
  slope_list = c()
  
  # A model for each species is calibrated
  for(species in species_list){
    trees_valid_sp = trees_valid[trees_valid[[variableSpecies]] == species, ]
    
    # if the tree of the species is alone, we use information from other species trees
    if(dim(trees_valid_sp)[1] <= 0){
      warning("The trees to interpolate does not have valid trees from their species. Hence, all trees are used.")
      tree_table$species_bis = 1
      return(computeMissingVariableProportionnaly(tree_table, variableX, variableY, "species_bis", variableYisMissing, method))
    }
    
    
    if(method == "lm"){
      
      # If the number of trees of this species is insuficient, we use a proportionnal model
      if(dim(trees_valid_sp)[1] <= 2){
        warning("The number of trees is not sufficient for linear model. A proportionnal model is used.")
        return(computeMissingVariableProportionnaly(tree_table, variableX, variableY, variableSpecies, variableYisMissing, method = "proportionnal_mean"))
      }
      
      model = lm(data = trees_valid_sp, as.formula(paste0(variableY, " ~ ", variableX) ))
      intercept = model$coefficients[[1]]
      slope = model$coefficients[[2]]
      
      if(slope < 0 | intercept < 0){
        warning(paste0("Results from linear model are too extreme. Hence a proportionnal model is used. slope: ", round(slope,3), " intercept: ", round(intercept, 3)))
        return(computeMissingVariableProportionnaly(tree_table, variableX, variableY, variableSpecies, variableYisMissing, method = "proportionnal_mean"))
      }
      
      intercept_list = c(intercept_list, intercept)
      slope_list = c(slope_list, slope)
    }else if(method == "proportionnal_lm"){ # use methods of least squared error. compared to proportionnal_mean distant values "pushes more"
      model = lm(data = trees_valid_sp, as.formula(paste0(variableY, " ~ 0 + ", variableX) ))
      ratio = model$coefficients
      slope_list = c(slope_list, ratio)
      intercept_list = c(intercept_list, 0)
    }else if(method == "proportionnal_mean"){ # use mean proportion. compared to proportionnal_lm all values proportions have "the same weight"
      ratio = mean(trees_valid_sp[[variableY]] / trees_valid_sp[[variableX]])
      slope_list = c(slope_list, ratio)
      intercept_list = c(intercept_list, 0)
    }else{
      stop("Valid methods are lm, proportionnal_lm and proportionnal_mean")
    }
  }
  
  # find the coefficients for each tree
  missingTreesSpecies = tree_table[variableYisMissing, ][[variableSpecies]]
  missingTreesInSpeciesList = match(missingTreesSpecies, species_list)
  
  interceptPerMissingTrees = intercept_list[ missingTreesInSpeciesList ]
  slopePerMissingTrees = slope_list[ missingTreesInSpeciesList ]
  
  # compute missing values with linear model for species
  tree_table[variableYisMissing, ][[variableY]] = interceptPerMissingTrees + tree_table[variableYisMissing, ][[variableX]] * slopePerMissingTrees
  
  return(tree_table)
}










# From a table of irregular, plurispecific trees, compute a table with regular trees of a given species
# idSpecies is the species identifier found in treesTable$sp
computeRegDemoMonoSpInventory = function(treesTable, cellsTable, idSpecies, allele1_1, allele2_1, allele1_2, allele2_2, pdgInventoryOptions, CASTANEASpeciesFile, thisScriptOptions, LAI_LIDAR){
  
  PDG_simulation_parameters = pdgInventoryOptions$PDG_simulation_parameters
  demographic_parameters = pdgInventoryOptions$demographic_parameters
  plotLocationParameters = pdgInventoryOptions$plotLocationParameters
  
  idSp = demographic_parameters$idSp
  
  slope_deg = plotLocationParameters$slope_deg
  aspect_deg = plotLocationParameters$aspect_deg
  referenceAltitude = plotLocationParameters$referenceAltitude
  
  nTree_regularization = thisScriptOptions$nTree_regularization
  nbAltitudeClasses = thisScriptOptions$nbAltitudeClasses
  useSpeciesLAIduringRegularization = thisScriptOptions$useSpeciesLAIduringRegularization
  targetCellWidth = PDG_simulation_parameters$targetCellWidth
  
  
  
  res = list()
  
  if(idSpecies == idSp[1]){
    allele1_1 = tibble( data.frame(t(allele1_1[1, ])), .rows = nTree_regularization)
    allele2_1 = tibble( data.frame(t(allele2_1[1, ])), .rows = nTree_regularization)
  }else{
    allele1_2 = tibble( data.frame(t(allele1_2[1, ])), .rows = nTree_regularization)
    allele2_2 = tibble( data.frame(t(allele2_2[1, ])), .rows = nTree_regularization)
  }
  
  # dichotomy of tree table
  trees_sp = subset(treesTable, sp == idSpecies)
  
  
  # count trees
  nTrees = dim(treesTable)[1]
  nTrees_sp = dim(trees_sp)[1]
  
  # count stem_surface
  stem_surface = sum(pi * (treesTable$dbh / 2)**2) # cm2
  stem_surface_sp = sum(pi * (trees_sp$dbh / 2)**2) # cm2
  
  # compute mean height to conserve mean tree volume with regular trees
  volume_per_tree_sp = trees_sp$height * 100 * pi * (trees_sp$dbh / 2)**2 # in cm3
  meanHeight_sp = sum(volume_per_tree_sp) / stem_surface_sp * 0.01 # in m
  
  # mixed plot attributes
  standArea = dim(cellsTable)[1] * targetCellWidth ** 2 # m2
  nha = dim(treesTable)[1] / (standArea/10000) # tree/ha
  gha = stem_surface / standArea # cm2 / m2
  
  proportion_sp = stem_surface_sp / stem_surface
  
  # species plots attributes
  standArea_sp = standArea * proportion_sp
  
  # Gha (with conservation of Gha between mixed inventories and dichotomized inventories)
  mean_basal_area_sp = stem_surface_sp / nTrees_sp
  
  gha_sp = nTrees_sp * mean_basal_area_sp / standArea_sp
  
  error_sp = gha_sp - gha
  
  if(abs(error_sp) > 1e-8){
    stop("Gha was not conserved for species during DBH regularization")
  }
  
  # Regularisation with quadratic mean DBH (with conservation basal area)
  quadraticMeanDbh_sp = sqrt(mean(trees_sp$dbh**2))
  mean_basal_area_bis_sp = (pi * (quadraticMeanDbh_sp/2)**2)
  
  error_bis_sp = mean_basal_area_bis_sp - mean_basal_area_sp
  if(abs(error_bis_sp) > 1e-8){
    stop("Gha was not conserved for species during DBH regularization")
  }
  
  # Regularization of trees position
  # On conserve le Nha_sp, le Gha et le dbh_sp moyen (quadratiquement)
  # On place X arbres en carré et on calcule ensuite la surface correspondante à X arbres
  
  nTree_byrow = sqrt(nTree_regularization)
  basal_area_new_sp = nTree_regularization * mean_basal_area_bis_sp
  
  # new stand area
  standArea_new_sp = basal_area_new_sp / gha_sp
  standWidth_new_sp = sqrt(standArea_new_sp)
  
  # Les cellules sont organisées en carré et leur taille est modifiée pour matcher exactement avec la surface de chaque espèce
  # Gha est ainsi conservé exactement
  # La largeur des cellules doit être proche de targetCellWidth
  nCellSide_sp = round(standWidth_new_sp / targetCellWidth)
  cellWidth_sp = standWidth_new_sp / nCellSide_sp
  
  # a. placer les cellules
  
  cells_sp = computeCellInfos(cellsTable[1,], nCellSide_sp)
  
  
  # b. placer les arbres
  
  # creer un nouvel arbre à partir des données des autres
  tmp_trees = rbind(treesTable, NA)
  tmp_trees[dim(tmp_trees)[1], ]$sp = idSpecies
  tmp_trees[dim(tmp_trees)[1], ]$dbh = quadraticMeanDbh_sp
  
  # cr-08.11.2032 set a computed average height to conserve overall volumes
  tmp_trees[dim(tmp_trees)[1], ]$height = meanHeight_sp
  
  # tmp_trees = computeMissingVariableProportionnaly(tmp_trees, "dbh", "height", "sp", 
  #                                                  variableYisMissing = c(rep(F, dim(treesTable)[1]), T),
  #                                                  method = "lm")
  tmp_trees = computeMissingVariableProportionnaly(tmp_trees, "height", "hcb", "sp", 
                                                   variableYisMissing = c(rep(F, dim(treesTable)[1]), T),
                                                   method = "lm")
  tmp_trees = computeMissingVariableProportionnaly(tmp_trees, "dbh", "age", "sp", 
                                                   variableYisMissing = c(rep(F, dim(treesTable)[1]), T),
                                                   method = "lm")
  
  tmp_trees$age = round(tmp_trees$age)
  
  # créer un tableau et le remplir
  trees_new_sp = tibble(tmp_trees[dim(tmp_trees)[1], ], .rows = nTree_regularization)
  
  trees_new_sp$idTree = 1:nTree_regularization
  trees_new_sp$pop = 1
  trees_new_sp$LAIinit = 0
  
  # placer les arbres
  tree_spacing_sp = standWidth_new_sp / nTree_byrow
  trees_new_sp$x_abs = tree_spacing_sp / 2 + rep(0:(nTree_byrow-1) * tree_spacing_sp, times = nTree_byrow)
  trees_new_sp$y_abs = tree_spacing_sp / 2 + rep(0:(nTree_byrow-1) * tree_spacing_sp, each = nTree_byrow)
  
  # ggplot(trees_new_sp1, aes(x = x_abs, y = y_abs, label = idTree)) + geom_point() + geom_label() + 
  #   coord_cartesian(xlim = c(0, standWidth_new_sp1), ylim = c(0, standWidth_new_sp1))
  
  
  # find cell and position for every trees
  for(i in 1:nTree_regularization){
    tree = trees_new_sp[i, ]
    cell_and_position = find_cell_and_rpositon(tree$x_abs, tree$y_abs, cells_sp, cellWidth_sp)
    trees_new_sp[i, ]$cell_IDs = cell_and_position[1]
    trees_new_sp[i, ]$xnorm = cell_and_position[2]
    trees_new_sp[i, ]$ynorm = cell_and_position[3]
  }
  
  
  nbColumns_sp = nCellSide_sp
  nbLines_sp = nCellSide_sp
  xMin_sp = 0; xMax_sp = xMin_sp + cellWidth_sp * nbColumns_sp ; yMin_sp = 0; yMax_sp = yMin_sp + cellWidth_sp * nbLines_sp
  
  altitudeClasses_sp = getAltittudeClases(xMin_sp, xMax_sp, yMin_sp, yMax_sp, slope_deg, aspect_deg, referenceAltitude, nbAltitudeClasses)
  
  
  
  
  
  # redefine species LAI
  standArea = dim(cellsTable)[1] * targetCellWidth **2
  LAI_sp = computeLAI_specific(trees, idSpecies, LAI_LIDAR, standArea, standArea_sp, CASTANEASpeciesFile, demographic_parameters)
  
  if(useSpeciesLAIduringRegularization){
    LAI_reg = LAI_sp
  }else{
    LAI_reg = LAI_LIDAR
  }
  
  
  res$cells_sp = cells_sp
  res$trees_sp = trees_new_sp
  res$cellWidth_sp = cellWidth_sp
  res$nbLines_sp = nbLines_sp
  res$nbColumns_sp = nbColumns_sp
  res$altitudeClasses_sp = altitudeClasses_sp
  res$allele1_1 = allele1_1
  res$allele2_1 = allele2_1
  res$allele1_2 = allele1_2
  res$allele2_2 = allele2_2
  res$LAI_reg = LAI_reg
  res$proportion_sp = proportion_sp
  res$quadraticMeanDbh_sp = quadraticMeanDbh_sp
  res$gha = gha # gha for mixed plot and gha for demographically regular monospecific plot
  
  return(res)
  
} # end computeRegDemoMonoSpInventory


# compute tree and cell tables for a demographically regular and plurispecific inventory
# from a demographuccaly irregular and plurispecicic inventory
computeRegDemoPluriSpInventory = function(trees_table, cells_table, allele1_1, allele2_1, allele1_2, allele2_2, pdgInventoryOptions, CASTANEASpeciesFile, thisScriptOptions, LAI_LIDAR){
  
  PDG_simulation_parameters = pdgInventoryOptions$PDG_simulation_parameters
  demographic_parameters = pdgInventoryOptions$demographic_parameters
  plotLocationParameters = pdgInventoryOptions$plotLocationParameters
  
  idSp = demographic_parameters$idSp
  
  
  slope_deg = plotLocationParameters$slope_deg
  aspect_deg = plotLocationParameters$aspect_deg
  referenceAltitude = plotLocationParameters$referenceAltitude
  
  nTree_regularization = thisScriptOptions$nTree_regularization
  nbAltitudeClasses = thisScriptOptions$nbAltitudeClasses
  useSpeciesLAIduringRegularization = thisScriptOptions$useSpeciesLAIduringRegularization
  targetCellWidth = PDG_simulation_parameters$targetCellWidth
  
  res = list()
  
  # get regularized inventories for each species
  res_species = list()
  
  # new trees table
  oneTreePerSpecies = NULL
  
  nTree_sp = c()
  propTree_sp = c()
  meanDbh_sp = c()
  meanTreeBA_sp = c()
  propBA_sp = c()
  basalarea_sp = c()
  surface_sp = c()
  idSp_bis = c()
  gha = NULL
  
  # first pass on species
  #   we compute regular monospecific inventories to have value of proportion and quadratic mean dbh
  for(idSpecies in unique(trees_table$sp)){
    code_asp = paste0("sp_", idSpecies)
    res_species[[code_asp]] = computeRegDemoMonoSpInventory(treesTable = trees_table, cellsTable = cells_table, idSpecies = idSpecies, 
                                                            allele1_1, allele2_1, allele1_2, allele2_2, 
                                                            pdgInventoryOptions, CASTANEASpeciesFile, thisScriptOptions, LAI_LIDAR)
    
    # get one examplary of each tree
    oneTree_asp = res_species[[code_asp]]$trees_sp[1, ]
    oneTreePerSpecies = rbind(oneTreePerSpecies, oneTree_asp)
    
    # get characteristic of stands
    idSp_bis = c(idSp_bis, oneTree_asp$sp)
    meanDbh_sp = c(meanDbh_sp, res_species[[code_asp]]$quadraticMeanDbh_sp)
    propBA_sp = c(propBA_sp, res_species[[code_asp]]$proportion_sp)
    
    gha = res_species[[code_asp]]$gha
  }
  
  meanTreeBA_sp = pi * (meanDbh_sp/2)**2
  
  # combien d'arbre pour la premiere espèce ?
  # propTreeFirstSp_th, with 'th' for theorical (because we will have to round the real number of trees)
  if( length(unique(trees_table$sp)) == 1 ){
    propTreeFirstSp_th = 1
  }else if(length(unique(trees_table$sp)) == 2){
    prop1 = propBA_sp[1]
    prop2 = propBA_sp[2]
    treeBA1 = meanTreeBA_sp[1]
    treeBA2 = meanTreeBA_sp[2]
    
    # this formulas gives the ratio of trees of the species one depending of the Basal area proportion and meanTreeBasalArea of each species
    propTreeFirstSp_th = 1 / ( 1 + prop2 / prop1 * treeBA1 / treeBA2)
  }else{
    stop("If number of species if greater than 2, use the general formula (not implemented, see comments)")
    
    # This is the general case formula with N species (demonstration in Cahier n°8 of Camille Rouet 27/09/2023) :
    # propTreeFirstSp_th = 1 / (   1
    #                         + prop2 / prop1 * treeBA1 / treeBA2 
    #                         + prop3 / prop1 * treeBA1 / treeBA3
    #                         + ...
    #                         + propN / prop1 * treeBA1 / treeBAN)
  }
  
  nTreeFirstSp = round(propTreeFirstSp_th * nTree_regularization)
  nTree_sp = c(nTreeFirstSp, nTree_regularization - nTreeFirstSp)
  propTree_sp = nTree_sp / nTree_regularization
  basalarea_sp = nTree_sp * meanTreeBA_sp
  surface_sp = basalarea_sp / gha
  
  # Placement des arbres un par un pour maximiser la dispersion
  
  
  # # Strategie 1 : Il faut : nTree_sp le nombre d'arbre de chaque espèce (a decrementer a chaque ajout d'arbre)
  # # period_sp Le temps de roulement avant d'ajouter un nouvel arbre pour chaque espece (a decrementer à chaque ajout d'arbre, a remettre à zero quand on place l'arbre d'une espèce)
  # 
  # # 1. on fait baisser le temps de roulement de tous jusqu'à ce que le plus petit soit à 1
  # # 2. a chaque pas, on fait baisser tous les temps de roulement
  # #   on ajoute un arbre par espèce qui est à 0
  # #   les espèces arrivées à zéro sont remis à leur chiffre
  # #   on retire un arbre aux arbres decompté de l'espèce qu'on a ajouté
  # #   si une espèce arrivé à zéro pour son nombre d'arbre ?
  # 
  # period_sp = round(1/propTree_sp)
  # nTree_sp_decremental = nTree_sp
  # period_sp_decremental = period_sp # the minimum value will be at one
  # 
  # trees_new = NULL
  # 
  # while(sum(nTree_sp_decremental) != 0){
  #   period_sp_decremental = period_sp_decremental - min(period_sp_decremental) # the minimum value will be at 0
  #   
  #   for(i_sp in 1:length(period_sp_decremental)){
  #     
  #     if(period_sp_decremental[i_sp] == 0 & nTree_sp_decremental[i_sp] > 0){
  #       one_tree_sp = oneTreePerSpecies[oneTreePerSpecies$sp == idSp_bis[i_sp], ]
  #       trees_new = rbind(trees_new, one_tree_sp)
  #       nTree_sp_decremental[i_sp] = nTree_sp_decremental[i_sp] - 1
  #       period_sp_decremental[i_sp] = period_sp[i_sp]
  #     }
  #   }
  # }
  
  
  # # check the number of trees
  # if(dim(trees_new)[1] != nTree_regularization){
  #   stop("Placing trees gave the wrong number of trees !")
  # }
  # 
  # # check the number of trees per species
  # for(asp in unique(trees_new$sp)){
  #   nTree_asp_2 = dim(trees_new[trees_new$sp == asp, ])[1]
  #   if(nTree_asp_2 != nTree_sp[idSp_bis == asp]){
  #     stop("Placing trees gave the wrong number of trees for a species !")
  #   }
  # }
  # 
  # trees_new$idTree = 1:nTree_regularization
  # 
  # 
  # 
  # 
  
  
  
  
  # Stratégie 2 : je place tous les arbres espèce par espèce un-à-un, PUIS, je place les arbres restants au hasard
  trees_new = NULL
  nTree_sp_decremental = nTree_sp
  i_sp = 1
  # place all trees in row until one species is at 0
  while(min(nTree_sp_decremental) != 0){ # if no species reach 0
    
    if(nTree_sp_decremental[i_sp] > 0){
      one_tree_asp = oneTreePerSpecies[oneTreePerSpecies$sp == idSp_bis[i_sp], ]
      trees_new = rbind(trees_new, one_tree_asp)
      nTree_sp_decremental[i_sp] = nTree_sp_decremental[i_sp] - 1
    }
    
    i_sp = i_sp%%length(idSp_bis) + 1 # increment through species indefinitly
  }
  
  
  # fill randomly with left trees
  for(i_sp in which(nTree_sp_decremental != 0)){
    one_tree_asp = oneTreePerSpecies[oneTreePerSpecies$sp == idSp_bis[i_sp], ]
    
    for(rep in 1:nTree_sp_decremental[i_sp]){
      
      nbPlacedTrees =  dim(trees_new)[1]
      if(is.null(trees_new)){
        nbPlacedTrees = 0
      }
      max_newtree_position = nbPlacedTrees + 1
      
      index_reputting = sample(1:max_newtree_position, size = 1)
      
      trees_new_split1 = NULL
      if(index_reputting != 1){
        trees_new_split1 = trees_new[1:(index_reputting-1), ]
      }
      
      trees_new_split2 = NULL
      if(index_reputting != max_newtree_position){
        trees_new_split2 = trees_new[index_reputting:dim(trees_new)[1], ]
      }
      
      trees_new = rbind(trees_new_split1, one_tree_asp, trees_new_split2)
    }
  }
  
  trees_new$idTree = 1:nTree_regularization
  
  
  
  # check the number of trees
  if(dim(trees_new)[1] != nTree_regularization){
    stop("Placing trees gave the wrong number of trees !")
  }
  
  
  
  # check the number of trees per species
  for(asp in unique(trees_new$sp)){
    nTree_asp_2 = dim(trees_new[trees_new$sp == asp, ])[1]
    if(nTree_asp_2 != nTree_sp[which(idSp_bis == asp)]){
      stop("Placing trees gave the wrong number of trees for a species !")
    }
  }
  
  
  
  # Ajouter un peu de décalage au hasard pour éviter d'avoir trop de colonne et ligne d'arbre de la même espèce
  # tout en gardant le mixage
  nSwap = 6
  swapLength = 1
  for(i in 1:nSwap){
    # maxSwapLength = round(nTree_regularization / 4)
    indices = 1:nTree_regularization
    
    index_retiring = sample(1:(nTree_regularization - swapLength + 1), size = 1)
    
    retireTreeIndices = index_retiring:(index_retiring+swapLength-1)
    retireTreeBool = indices %in% retireTreeIndices
    
    
    trees_kept = trees_new[retireTreeBool, ]
    trees_new = trees_new[!retireTreeBool, ]
    
    index_reputting = sample(1:(nTree_regularization - swapLength + 1), size = 1)
    
    trees_new_split1 = NULL
    if(index_reputting != 1){
      trees_new_split1 = trees_new[1:(index_reputting-1), ]
    }
    
    trees_new_split2 = NULL
    if(index_reputting != nTree_regularization){
      trees_new_split2 = trees_new[index_reputting:dim(trees_new)[1], ]
    }
    trees_new = rbind(trees_new_split1, trees_kept, trees_new_split2)
    
  }
  
  # recheck the number of trees
  if(dim(trees_new)[1] != nTree_regularization){
    stop("Swaping trees gave the wrong number of trees !")
  }
  
  # recheck the number of trees per species
  for(asp in unique(trees_new$sp)){
    nTree_asp_2 = dim(trees_new[trees_new$sp == asp, ])[1]
    if(nTree_asp_2 != nTree_sp[which(idSp_bis == asp)]){
      stop("Swaping trees gave the wrong number of trees for a species !")
    }
  }
  
  
  # recreer les tableaux alleliques à la bonne taille
  for(asp in unique(trees_table$sp)){
    nTree_asp = nTree_sp[which(idSp_bis == asp)]
    if(asp == idSp[1]){
      allele1_1 = tibble( data.frame(t(allele1_1[1, ])), .rows = nTree_asp)
      allele2_1 = tibble( data.frame(t(allele2_1[1, ])), .rows = nTree_asp)
    }else{
      allele1_2 = tibble( data.frame(t(allele1_2[1, ])), .rows = nTree_asp)
      allele2_2 = tibble( data.frame(t(allele2_2[1, ])), .rows = nTree_asp)
    }
  }
  
  
  
  
  # new plot properties
  plotArea_new = sum(surface_sp)
  
  plotWidth_new = sqrt(plotArea_new)
  
  # new cells table
  nCellSide_new = round(plotWidth_new / targetCellWidth)
  cellWidth_new = plotWidth_new / nCellSide_new
  cells_new = computeCellInfos(cells_table[1,], nCellSide_new)
  
  # Vrai placement des arbres
  nTree_byrow = sqrt(nTree_regularization)
  tree_spacing_sp = plotWidth_new / nTree_byrow
  trees_new$x_abs = tree_spacing_sp / 2 + rep(0:(nTree_byrow-1) * tree_spacing_sp, times = nTree_byrow)
  trees_new$y_abs = tree_spacing_sp / 2 + rep(0:(nTree_byrow-1) * tree_spacing_sp, each = nTree_byrow)
  
  
  
  
  # find cell and position for every trees
  for(i in 1:nTree_regularization){
    tree = trees_new[i, ]
    cell_and_position = find_cell_and_rpositon(tree$x_abs, tree$y_abs, cells_new, cellWidth_new)
    trees_new[i, ]$cell_IDs = cell_and_position[1]
    trees_new[i, ]$xnorm = cell_and_position[2]
    trees_new[i, ]$ynorm = cell_and_position[3]
  }
  
  # recompute altitude classes
  nbColumns_new = nCellSide_new
  nbLines_new = nCellSide_new
  xMin_new = 0; xMax_new = xMin_new + cellWidth_new * nbColumns_new ; yMin_new = 0; yMax_new = yMin_new + cellWidth_new * nbLines_new
  
  altitudeClasses_new = getAltittudeClases(xMin_new, xMax_new, yMin_new, yMax_new, slope_deg, aspect_deg, referenceAltitude, nbAltitudeClasses)
  
  res$cells_new = cells_new
  res$trees_new = trees_new
  res$cellWidth_new = cellWidth_new
  res$nbLines_new = nbLines_new
  res$nbColumns_new = nbColumns_new
  res$altitudeClasses_new = altitudeClasses_new
  res$allele1_1 = allele1_1
  res$allele2_1 = allele2_1
  res$allele1_2 = allele1_2
  res$allele2_2 = allele2_2
  
  
  return(res)
  
} # end computeRegDemoPluriSpInventory



# From a table of irregular plurispecific trees, gives a table of irrgular monospecifi trees of species idSpecies
# not-targeted trees are converted in the target species (their idTree is marked to be reckonizable)
# dbh, and tree position are conserved 
# New trees height, hcb and age are computed with linear interpolation based on dbh
computeIrregDemoMonoSpInventory = function(treesTable, cellsTable, idSpecies, allele1_1, allele2_1, allele1_2, allele2_2, pdgInventoryOptions, CASTANEASpeciesFile, LAI_LIDAR){
  
  
  PDG_simulation_parameters = pdgInventoryOptions$PDG_simulation_parameters
  demographic_parameters = pdgInventoryOptions$demographic_parameters
  
  cellWidth = PDG_simulation_parameters$targetCellWidth
  
  idSp = demographic_parameters$idSp
  
  res = list()
  
  # isolate species trees
  trees_sp = subset(treesTable, sp == idSpecies)
  trees_another_sp = subset(treesTable, sp != idSpecies)
  
  # count trees
  nTrees = dim(treesTable)[1]
  nTrees_sp = dim(trees_sp)[1]
  
  # count stem_surface
  stem_surface = sum(pi * (treesTable$dbh / 2)**2) # cm2
  stem_surface_sp = sum(pi * (trees_sp$dbh / 2)**2) # cm2
  stem_surface_another_sp = stem_surface - stem_surface_sp
  
  # mixed plot attributes
  standArea = dim(cellsTable)[1] * cellWidth ** 2 # m2
  nha = dim(treesTable)[1] / (standArea/10000) # tree/ha
  gha = stem_surface / standArea # cm2 / m2
  
  proportion_sp = stem_surface_sp / stem_surface
  
  # species plots attributes
  standArea_sp = standArea * proportion_sp
  
  # dbh and mean basal area
  mean_basal_area_sp = stem_surface_sp / nTrees_sp
  quadraticMeanDbh_sp = sqrt(mean(trees_sp$dbh**2))
  
  LAI_sp = computeLAI_specific(trees, idSpecies, LAI_LIDAR, standArea, standArea_sp, CASTANEASpeciesFile, demographic_parameters)
  
  
  
  method = "1"
  
  
  if(method == "1"){
    # Choix 1 : Transform species into another to get monospecific
    treesTable$targetSpecies = treesTable$sp == idSpecies
    treesTable$sp = idSpecies
    treesTable[!treesTable$targetSpecies, ]$idTree = treesTable[!treesTable$targetSpecies, ]$idTree + 1000000
    
    treesTable = computeMissingVariableProportionnaly(treesTable, "dbh", "height", "sp", !treesTable$targetSpecies, method = "proportionnal_mean")
    treesTable = computeMissingVariableProportionnaly(treesTable, "dbh", "height2013", "sp", !treesTable$targetSpecies, method = "proportionnal_mean")
    treesTable = computeMissingVariableProportionnaly(treesTable, "dbh", "hcb", "sp", !treesTable$targetSpecies, method = "proportionnal_mean")
    treesTable = computeMissingVariableProportionnaly(treesTable, "dbh", "hcb2013", "sp", !treesTable$targetSpecies, method = "proportionnal_mean")
    treesTable = computeMissingVariableProportionnaly(treesTable, "dbh", "age", "sp", !treesTable$targetSpecies, method = "proportionnal_mean")
    treesTable = computeMissingVariableProportionnaly(treesTable, "dbh", "age2013", "sp", !treesTable$targetSpecies, method = "proportionnal_mean")
    treesTable$age = round(treesTable$age)
  }else{ # in dev..
    # Choix 2 : on ajoute des arbres de dbh moyen de l'espèce de sorte à combler la surface terrière manquante.
    # On place les individus aléatoires selon plusieurs jets, de sorte à maximiser l'écart entre les arbres
    nNewTrees = round(stem_surface_another_sp / mean_basal_area_sp)
    basal_area_newTrees = stem_surface_another_sp / nNewTrees
    dbh_newTrees = 2 * sqrt(basal_area_newTrees / pi)
    
    nResampling = 1000
    
    sampleTable = tibble(i_sample = 0, x_abs = 0, y_abs = 0, .rows = 0)
    sampleSumDistance = tibble(i_sample = 0, sumDistance = 0, .rows = 0)
    
    for(i_sample in 1:nResampling){
      for(i_tree in 1:nNewTrees){
        
      }
      # compute sumDistance
    }
    
    
    # vérifier la conservation du LAI de l'espèce
    # vérifier le gha
  }
  
  nTrees_new = dim(treesTable)[1]
  if(idSpecies == idSp[1]){
    allele1_1 = tibble( data.frame(t(allele1_1[1, ])), .rows = nTrees_new)
    allele2_1 = tibble( data.frame(t(allele2_1[1, ])), .rows = nTrees_new)
  }else{
    allele1_2 = tibble( data.frame(t(allele1_2[1, ])), .rows = nTrees_new)
    allele2_2 = tibble( data.frame(t(allele2_2[1, ])), .rows = nTrees_new)
  }
  
  res$trees_sp = treesTable
  res$allele1_1 = allele1_1
  res$allele2_1 = allele2_1
  res$allele1_2 = allele1_2
  res$allele2_2 = allele2_2
  res$LAI_sp = LAI_sp
  return(res)
} # end computeIrregDemoMonoSpInventory




printInventoryGraph = function(treesTable, folderPath, fileName, forestPlotWidth, cellWidth, plotCrown = FALSE, useRDI = FALSE, ratioOnRadius = 1){
  
  if(plotCrown){
    # add crown projection
    if(mean(unique(treesTable$speciesName) %in% c("hetre",  "sapin")) != 1){
      stop("Species in treesTable are not hetre or sapin")
    }
    
    # from CastaneaSpecies_2022_05.txt, crownRadius1 and crownRadius2
    treesTable$crownRadius = ifelse(treesTable$speciesName == "hetre", 0.1082, 0.08151) * treesTable$dbh + ifelse(treesTable$speciesName == "hetre", 1.04, 0.69535) 
    if(useRDI){
      standArea_ha = forestPlotWidth**2 / 10000
      nha = dim(treesTable)[1] / standArea_ha
      meanDbh = mean(treesTable$dbh)
      RDIbeech = max(1, nha / exp(12.95  - 1.941 * log(meanDbh)))
      RDIfir = max(1, nha/ exp(12.621 - 1.779 * log(meanDbh)))
      
      # convert to crown area
      treesTable$crownArea = pi*treesTable$crownRadius**2
      
      # RDI effect
      treesTable$crownArea = ifelse(treesTable$speciesName == "hetre", treesTable$crownArea / RDIbeech, treesTable$crownArea / RDIfir)
      treesTable$crownArea = vapply(treesTable$crownArea, FUN = function(x) max(5, x), FUN.VALUE = 0)
      
      # convert to crown radius
      treesTable$crownRadius = sqrt(treesTable$crownArea  / pi)
    }
    
    npc = 40 # number of point circle
    data_polygon = tibble(x = 0, y  = 0, id = 0, targetTrees = T, done = F, .rows = dim(treesTable)[1] * npc )
    for(i in 1:dim(treesTable[1])){
      tree = treesTable[i, ]
      cr = tree$crownRadius
      x_polygon = tree$x_abs + cr * cos( seq(0, 2 * pi, length.out = npc) ) 
      y_polygon = tree$y_abs + cr * sin( seq(0, 2 * pi, length.out = npc) ) 
      id_polygon = rep(tree$idTree, each = npc)
      data_polygon[((i-1) * npc  + 1) :(i * npc), ] = tibble(x = x_polygon, y  = y_polygon, id = id_polygon, targetTrees = T, done = T)
    }
  }
  
  pointSizeMax = max(treesTable$dbh) / 10
  pointSizeMin = min(treesTable$dbh) / 10
  
  aplot = ggplot(treesTable, aes(x = x_abs, y = y_abs, color = speciesName, size = dbh, shape = targetTrees)) + 
    scale_size_continuous(range = c(pointSizeMin, pointSizeMax)) + geom_point() +
    coord_cartesian(xlim = c(0, forestPlotWidth), ylim = c(0, forestPlotWidth)) + 
    scale_color_manual(breaks = c("hetre", "sapin"), values=c("darkcyan", "darkgreen"))+ 
    scale_shape_manual(breaks = c(TRUE, FALSE), values=c(19, 1))
  
  if(plotCrown){
    aplot = aplot + geom_polygon(data = data_polygon, mapping =  aes(x = x, y = y, group = id), fill = rgb(0,0,0,0.05), color = rgb(0,0,0, 0.1), size = 0.1)
  }
  
  
  nLin = floor(forestPlotWidth / cellWidth)
  x = rep(seq(0, forestPlotWidth, by = cellWidth), times = nLin+1)
  y = rep(seq(0, forestPlotWidth, by = cellWidth), each = nLin+1)
  pointTable = tibble(x_abs = x, y_abs = y)
  
  # add PDGCell grid point
  aplot + geom_point(data = pointTable, shape = 16, color = rgb(0,0,0,0.1), size = 0.5)
  
  if(!dir.exists(folderPath)){
    dir.create(path = folderPath)
  }
  
  aplot = aplot + guides(shape="none") +  theme(legend.position="bottom")
  imageRatio = 0.9
  height = 960
  width = height * imageRatio
  ggsave(plot = aplot, path = folderPath, filename = paste0(fileName, ".jpg"), dpi = 0.5 * 300, width = width, height = height, units = "px")
}



# Write the inventory for given lists of cells and table
# allele1_1 is allele1 for sp 1, of dimension nbInvidiual(sp1) * (nbLocus_FCRITBB + nbLocus_g1max + nbMsat+nbSNP)
# allele2_1 is allele2 for sp 1, same dimension.
# allele1_2 and allele2_2 are the same for sp 2 
writeInventoryPDG = function(inventoryFolder, inventoryName, description,
                             cellsTable, treesTable, 
                             pdgInventoryOptions,
                             altitudeClasses_, allele1_1, allele2_1, allele1_2, allele2_2, thisScriptOptions, ratioOnRadiusPlot = 1){
  
  if(!"pdgPlotParameters" %in% names(pdgInventoryOptions)){
    stop("pdgPlotParameters not in pdgInventoryOptions")
  }
  
  if(!"genetics_variables" %in% names(pdgInventoryOptions)){
    stop("genetics_variables not in pdgInventoryOptions")
  }
  
  # Load parameters tables
  demographic_parameters = pdgInventoryOptions$demographic_parameters
  genetics_parameters = pdgInventoryOptions$genetics_parameters
  PDG_simulation_parameters = pdgInventoryOptions$PDG_simulation_parameters
  distanceClassBoundsForSGS = pdgInventoryOptions$distance_class_bounds_for_SGS
  pdgPlotParameters = pdgInventoryOptions$pdgPlotParameters
  plotLocationParameters = pdgInventoryOptions$plotLocationParameters
  
  genetics_variables = pdgInventoryOptions$genetics_variables
  
  
  # reformat
  demographic_parameters$pollenKernel = ifelse(sapply(demographic_parameters$pollenKernel, FUN = function(x) x == TRUE), "true", "false")
  demographic_parameters$selfing = ifelse(sapply(demographic_parameters$selfing, FUN = function(x) x == TRUE), "true", "false")
  
  
  # LOAD genetics parameters
  gen_param_u = genetics_parameters[1,]
  numberOfGeneticParameters = gen_param_u$numberOfGeneticParameters
  numberOfTraits = gen_param_u$numberOfTraits
  nbLocus_FCRITBB=	gen_param_u$nbLocus_FCRITBB
  nbAllelePerLoc_FCRITBB= gen_param_u$nbAllelePerLoc_FCRITBB 
  interStepEnvirVar_FCRITBB= gen_param_u$interStepEnvirVar_FCRITBB  
  herit_FCRITBB=	gen_param_u$herit_FCRITBB
  target_mean_FCRITBB= gen_param_u$target_mean_FCRITBB 
  target_CV_FCRITBB=	gen_param_u$target_CV_FCRITBB
  alleleEffectMultiplCoeffFCRITBB= gen_param_u$alleleEffectMultiplCoeffFCRITBB 
  nbLocus_g1max= gen_param_u$nbLocus_g1max
  nbAllelePerLoc_g1max= gen_param_u$nbAllelePerLoc_g1max 
  interStepEnvirVar_g1max= gen_param_u$interStepEnvirVar_g1max 
  herit_g1max=	gen_param_u$herit_g1max
  target_mean_g1max= gen_param_u$target_mean_g1max     
  target_CV_g1max=	gen_param_u$target_CV_g1max     
  alleleEffectMultiplCoeffg1max= gen_param_u$alleleEffectMultiplCoeffg1max 
  nbMsat= gen_param_u$nbMsat
  nbAllPerMsat= gen_param_u$nbAllPerMsat
  nbSNP= gen_param_u$nbSNP
  nbAllPerSNP= gen_param_u$nbAllPerSNP
  distanceClassBoundsForSGS = pdgInventoryOptions$distance_class_bounds_for_SGS
  targetSd_FCRITBB=target_CV_FCRITBB*target_mean_FCRITBB
  targetVariance_FCRITBB=targetSd_FCRITBB^2
  totEnvirVar_FCRITBB=	targetVariance_FCRITBB*(1-herit_FCRITBB)  
  targetSd_g1max=target_CV_g1max*target_mean_g1max 
  targetVariance_g1max=targetSd_g1max^2
  totEnvirVar_g1max=	targetVariance_g1max*(1-herit_g1max)   
  
  # LOAD script options
  outputInConsole = thisScriptOptions$outputInConsole
  initYear = thisScriptOptions$initYear
  
  
  
  if("x_abs" %in% names(treesTable)){
    plotWidth = sqrt(pdgPlotParameters$nlin * pdgPlotParameters$ncol * pdgPlotParameters$cellWidth**2)
    
    treesTable$speciesName = demographic_parameters$speciesFrenchNames[ match(treesTable$sp, demographic_parameters$idSp) ]
    treesTable$targetTrees = treesTable$idTree < 1000000
    
    printInventoryGraph(treesTable, paste0(inventoryFolder, "0_plots_image"), inventoryName, plotWidth, pdgPlotParameters$cellWidth, F, T, ratioOnRadius = ratioOnRadiusPlot)
    printInventoryGraph(treesTable, paste0(inventoryFolder, "0_plots_image/withCrowns"), inventoryName, plotWidth, pdgPlotParameters$cellWidth, T, T, ratioOnRadius = ratioOnRadiusPlot)
  }
  
  
  # make sure that even if error occurs, output will be redirected to console
  on.exit(while(sink.number() > 0){sink()}, add = TRUE, after = FALSE)
  
  # Open writing
  if(!outputInConsole)
    sink(paste0(inventoryFolder, inventoryName))
  
  
  cat(description)
  
  #HEADERS and simple parameters / EN TETES et parametres simples
  cat("\n\n\t#<General>");
  
  # 4 Write PARAMETERS
  if(plotLocationParameters$slope_deg > 0){
    simulationType = "SIMULATION_TYPE_ALTITUDE_GRADIENT"
  }else{
    simulationType = "SIMULATION_TYPE_FIXED_ALTITUDE"
  }
  
  PDG_simulation_parameters_print = PDG_simulation_parameters[, -which(colnames(PDG_simulation_parameters) == "targetCellWidth")]
  PDG_simulation_parameters_print$initYear = NULL
  
  cat("\n\n\t#SIMULATION parameters");
  cat("\nyear = ", initYear, sep = "")
  cat("\n")
  cat( paste0( paste0(colnames(PDG_simulation_parameters_print), " = ", c(PDG_simulation_parameters_print[1, ])), collapse = "\n" ) )
  cat("\n")
  cat("\nsimulationType = ", simulationType, sep = "")
  
  if (simulationType=="SIMULATION_TYPE_ALTITUDE_GRADIENT") 
  {
    result="\naltitudeClasses={"
    for (i in 1:length(altitudeClasses_)) {
      result=paste(result,altitudeClasses_[i], ";",sep="")
    }
    rm(i)
    result=paste(result,"}",sep="")
    cat(result)
    
  }
  
  cat("\n\n\t#PLOT parameters\n")
  stringLineText=paste("origin = (", pdgPlotParameters$xorigin,",",pdgPlotParameters$yorigin,")",sep=" ")
  cat(stringLineText)
  cat("\nnLin =", pdgPlotParameters$nlin, sep = " ")
  cat("\nnCol =", pdgPlotParameters$ncol, sep= " ")
  stringLineText=paste("\ncellWidth =", pdgPlotParameters$cellWidth, sep=" ")
  cat(stringLineText)
  cat(paste("\nreferenceAltitude =", plotLocationParameters$referenceAltitude, sep=" "))
  cat(paste("\nslope_deg =", plotLocationParameters$slope_deg, sep=" "))
  cat(paste("\naspect_deg =", plotLocationParameters$aspect_deg, sep=" "))
  cat("\nlatitude =", plotLocationParameters$latitude, sep=" ")
  cat("\nlongitude =", plotLocationParameters$longitude, sep=" ")
  cat("\nstandLAI =", round(pdgPlotParameters$LAI,3), sep = " ")
  
  cat("\n\n\t#GENETIC parameters")
  cat(paste("\nnumberOfGeneticParameters =", numberOfGeneticParameters, "\nnumberOfTraits =", numberOfTraits, sep=" "))
  cat(paste("\nnbLocus_FCRITBB =", nbLocus_FCRITBB, "\nnbLocus_g1max =", nbLocus_g1max, sep=" "))
  cat(paste("\nnbLocusMicrosat =", nbMsat, "\nnbLocusNeutralSNP =", nbSNP, sep=" "))
  
  #cat(paste("\ndistanceClassBoundsForSGS =", distanceClassBoundsForSGS, sep=" "))
  result="\ndistanceClassBoundsForSGS = {"
  for (i in 1:length(distanceClassBoundsForSGS)) {
    result=paste(result,distanceClassBoundsForSGS[i], ";",sep="")
  }
  rm(i)
  result=paste(result,"}",sep="")
  cat(result)
  
  cat("\n")
  cat("\n#[4 COMPUTED PARAMETERS]")
  
  if(herit_FCRITBB == 0){
    mean_FCRITBBval= genetics_variables$mean_FCRITBBval
    sd_FCRITBBval= genetics_variables$sd_FCRITBBval
  }
  
  if(is.nan(mean_FCRITBBval) | is.na(mean_FCRITBBval)){
    mean_FCRITBBval = -1
  }
  if(is.nan(sd_FCRITBBval) | is.na(sd_FCRITBBval)){
    sd_FCRITBBval = -1
  }
  
  cat(paste("\nmeanInitialGValue_FCRITBB =", mean_FCRITBBval, "\nsdInitialGValue_FCRITBB =", round(sd_FCRITBBval,2), sep=" "))
  
  if(herit_g1max == 0){
    mean_g1maxVal= genetics_variables$mean_g1maxVal
    sd_g1maxVal= genetics_variables$sd_g1maxVal
  }
  if(is.nan(mean_g1maxVal) | is.na(mean_g1maxVal)){
    mean_g1maxVal = -1
  }
  if(is.nan(sd_g1maxVal) | is.na(sd_g1maxVal)){
    sd_g1maxVal = -1
  }
  cat(paste("\nmeanInitialGValue_g1max =", mean_g1maxVal, "\nsdInitialGValue_g1max =", round(sd_g1maxVal,2), sep=" "))
  
  cat("\n")
  cat(paste("\nmeanTargetValue_FCRITBB =", target_mean_FCRITBB))
  cat(paste("\n#targetedCV_FCRITBB =", target_CV_FCRITBB, sep=" "))
  cat(paste("\nmeanTargetValue_g1max =", target_mean_g1max))
  cat(paste("\n#targetedCV_g1max =", target_CV_g1max, sep=" "))
  
  cat("\n")
  cat(paste("\nalleleEffectMultiplCoeffFCRITBB =", alleleEffectMultiplCoeffFCRITBB, sep=" "))
  cat(paste("\nalleleEffectMultiplCoeffg1max =", alleleEffectMultiplCoeffg1max, sep=" "))
  
  
  CASTANEA_simulation_parameters = pdgInventoryOptions$CASTANEA_simulation_parameters
  CASTANEA_simulation_parameters = CASTANEA_simulation_parameters[-which(colnames(CASTANEA_simulation_parameters) == 'Ca')]
  cat("\n\n#CASTANEA parameters")
  cat("\n")
  cat( paste0( paste0(colnames(CASTANEA_simulation_parameters), " = ", c(CASTANEA_simulation_parameters[1, ])), collapse = "\n" ) )
  
  
  
  
  # 4 Genetic Map / carte genetiqu
  
  cat("\n\n#Genetic Map \n#speciesName\tspeciesId\tcastaneaCode\tgeneticMap\tallelesNuclear\tallelesMCytoplasmic\tallelesPCytoplasmic");
  cat("\n", demographic_parameters$speciesFullNames[1],"\t", demographic_parameters$idSp[1],"\t", demographic_parameters$fmIdSp[1],"\t{}\t{", sep = "")
  stringGenetMap=character(0)
  
  for (i in 1:nbLocus_FCRITBB)
  {
    stringGenetMap=paste(stringGenetMap,"[",sep="")
    for (j in 1:(nbAllelePerLoc_FCRITBB-1))
    {
      stringGenetMap=paste(stringGenetMap,j,",",sep="")
    }
    stringGenetMap=paste(stringGenetMap,nbAllelePerLoc_FCRITBB,"];",sep="")
  }
  rm(i,j)
  cat(stringGenetMap)
  stringGenetMap=character(0)
  
  for (i in 1:nbLocus_g1max)
  {
    stringGenetMap=paste(stringGenetMap,"[",sep="")
    for (j in 1:(nbAllelePerLoc_g1max-1))
    {
      stringGenetMap=paste(stringGenetMap,j,",",sep="")
    }
    stringGenetMap=paste(stringGenetMap,nbAllelePerLoc_g1max,"];",sep="")
  }
  rm(i,j)
  cat(stringGenetMap)
  stringGenetMap=character(0)
  
  if(nbMsat>0){
    for (i in 1:nbMsat)
    {
      stringGenetMap=paste(stringGenetMap,"[",sep="")
      for (j in 1:(nbAllPerMsat-1))
      {
        stringGenetMap=paste(stringGenetMap,10+j,",",sep="")
      }
      stringGenetMap=paste(stringGenetMap,10+nbAllPerMsat,"];",sep="")
    }
    cat(stringGenetMap)
  }
  
  stringGenetMap=character(0)
  for (i in 1:nbSNP)
  {
    stringGenetMap=paste(stringGenetMap,"[",sep="")
    for (j in 1:(nbAllPerSNP-1))
    {
      stringGenetMap=paste(stringGenetMap,j,",",sep="")
    }
    stringGenetMap=paste(stringGenetMap,nbAllPerSNP,"];",sep="")
  }
  cat(stringGenetMap)
  cat("}\t{}\t{}")
  
  # WARNING : code was copypasted was each species, but genetic parameters could differ..
  cat("\n", demographic_parameters$speciesFullNames[2],"\t", demographic_parameters$idSp[2],"\t", demographic_parameters$fmIdSp[2],"\t{}\t{", sep = "")
  stringGenetMap=character(0)
  
  for (i in 1:nbLocus_FCRITBB)
  {
    stringGenetMap=paste(stringGenetMap,"[",sep="")
    for (j in 1:(nbAllelePerLoc_FCRITBB-1))
    {
      stringGenetMap=paste(stringGenetMap,j,",",sep="")
    }
    stringGenetMap=paste(stringGenetMap,nbAllelePerLoc_FCRITBB,"];",sep="")
  }
  cat(stringGenetMap)
  stringGenetMap=character(0)
  
  for (i in 1:nbLocus_g1max)
  {
    stringGenetMap=paste(stringGenetMap,"[",sep="")
    for (j in 1:(nbAllelePerLoc_g1max-1))
    {
      stringGenetMap=paste(stringGenetMap,j,",",sep="")
    }
    stringGenetMap=paste(stringGenetMap,nbAllelePerLoc_g1max,"];",sep="")
  }
  cat(stringGenetMap)
  stringGenetMap=character(0)
  
  if(nbMsat>0){
    for (i in 1:nbMsat)
    {
      stringGenetMap=paste(stringGenetMap,"[",sep="")
      for (j in 1:(nbAllPerMsat-1))
      {
        stringGenetMap=paste(stringGenetMap,10+j,",",sep="")
      }
      stringGenetMap=paste(stringGenetMap,10+nbAllPerMsat,"];",sep="")
    }
    cat(stringGenetMap)
  }
  
  stringGenetMap=character(0)
  for (i in 1:nbSNP)
  {
    stringGenetMap=paste(stringGenetMap,"[",sep="")
    for (j in 1:(nbAllPerSNP-1))
    {
      stringGenetMap=paste(stringGenetMap,j,",",sep="")
    }
    stringGenetMap=paste(stringGenetMap,nbAllPerSNP,"];",sep="")
  }
  cat(stringGenetMap)
  cat("}\t{}\t{}")
  
  
  
  
  # 4 Species demographics PARAMETERS
  #Species demographic parameters
  cat("\n\n# Species demographic parameters \n")
  demographic_parameters_print = cbind(demographic_parameters[, colnames(demographic_parameters) == "idSp"], demographic_parameters[which(colnames(demographic_parameters) == "deltaSeed") : which(colnames(demographic_parameters) == "maxNbOfGrowingSeedlings")])
  
  # dont print these columns
  demographic_parameters_print$maxSeedNumber = NULL
  demographic_parameters_print$areaOfOneTree = NULL
  
  cat("#", paste0(colnames(demographic_parameters_print), collapse = "\t"), sep = "")
  
  for(speciesIndex in 1:dim(demographic_parameters_print)[1]){
    cat("\n")
    cat(paste0(demographic_parameters_print[speciesIndex, ], collapse = "\t"))
  }
  
  
  
  # 4 allele effects / effets allelique
  # WARNING : allele effects was copied for each species, but genetic could differ between species
  cat("\n\n# allele effects for each trait")
  cat("\n#speciesId\tparameter\teffectCoefficient\tnuclear allele effect\tmcyto allele effect\tpcyto allele effect\theritability\tenvironmentalVariance\tinterEnvironmentalVariance	")
  
  stringLocusEffect = paste0("\n", demographic_parameters$idSp[1], "\tFCRITBB\t1\t{")
  for (i in 1:nbLocus_FCRITBB)
  {
    if(herit_FCRITBB>0) {
      all1=(target_mean_FCRITBB/nbLocus_FCRITBB-pdgInventoryOptions$effect_FCRITBBdraw[i])/2
      all1=round(all1*alleleEffectMultiplCoeffFCRITBB)
      all2=(target_mean_FCRITBB/nbLocus_FCRITBB+pdgInventoryOptions$effect_FCRITBBdraw[i])/2
      all2=round(all2*alleleEffectMultiplCoeffFCRITBB)
      #hom1=-effect_FCRITBB[i]
      #hom2=+effect_FCRITBB[i]
    } else {
      all1=0 
      all2=0
    }
    
    stringLocusEffect=paste(stringLocusEffect,"[",(i),",",all1,",",all2,"]",sep="")
    if (i<nbLocus_FCRITBB) stringLocusEffect=paste(stringLocusEffect,";",sep="")
  }
  stringLocusEffect=paste(stringLocusEffect,"}\t{}\t{}\t",herit_FCRITBB,"\t",sep="")
  if (totEnvirVar_FCRITBB==0) {
    stringLocusEffect=paste(stringLocusEffect,"0\t0")
  } else {
    stringLocusEffect=paste(stringLocusEffect,round(totEnvirVar_FCRITBB*alleleEffectMultiplCoeffFCRITBB),"\t0")
  }
  cat (stringLocusEffect)                        
  
  stringLocusEffect = paste0("\n", demographic_parameters$idSp[1], "\tg1max\t1\t{")
  for (i in 1:nbLocus_g1max)
  {
    if(herit_g1max>0) {
      all1=(target_mean_g1max/nbLocus_g1max-pdgInventoryOptions$effect_g1maxdraw[i])/2
      all1=round(all1*alleleEffectMultiplCoeffg1max)
      all2=(target_mean_g1max/nbLocus_g1max+pdgInventoryOptions$effect_g1maxdraw[i])/2
      all2=round(all2*alleleEffectMultiplCoeffg1max)
      #hom1=-effect_g1max[i]
      #hom2=+effect_g1max[i]
    } else {
      all1=0 
      all2=0
    }
    
    stringLocusEffect=paste(stringLocusEffect,"[",(i+10),",",all1,",",all2,"]",sep="")
    if (i<nbLocus_g1max) stringLocusEffect=paste(stringLocusEffect,";",sep="")
  }
  stringLocusEffect=paste(stringLocusEffect,"}\t{}\t{}\t",herit_g1max,"\t",sep="")
  if (totEnvirVar_g1max==0) {
    stringLocusEffect=paste(stringLocusEffect,"0\t0")
  } else {
    stringLocusEffect=paste(stringLocusEffect,round(totEnvirVar_g1max*alleleEffectMultiplCoeffg1max),"\t0")
  }
  cat (stringLocusEffect)  
  
  # Allele effect for species 2
  stringLocusEffect = paste0("\n", demographic_parameters$idSp[2], "\tFCRITBB\t1\t{")
  for (i in 1:nbLocus_FCRITBB)
  {
    if(herit_FCRITBB>0) {
      all1=(target_mean_FCRITBB/nbLocus_FCRITBB-pdgInventoryOptions$effect_FCRITBBdraw[i])/2
      all1=round(all1*alleleEffectMultiplCoeffFCRITBB)
      all2=(target_mean_FCRITBB/nbLocus_FCRITBB+pdgInventoryOptions$effect_FCRITBBdraw[i])/2
      all2=round(all2*alleleEffectMultiplCoeffFCRITBB)
      #hom1=-effect_FCRITBB[i]
      #hom2=+effect_FCRITBB[i]
    } else {
      all1=0 
      all2=0
    }
    
    stringLocusEffect=paste(stringLocusEffect,"[",(i),",",all1,",",all2,"]",sep="")
    if (i<nbLocus_FCRITBB) stringLocusEffect=paste(stringLocusEffect,";",sep="")
  }
  stringLocusEffect=paste(stringLocusEffect,"}\t{}\t{}\t",herit_FCRITBB,"\t",sep="")
  if (totEnvirVar_FCRITBB==0) {
    stringLocusEffect=paste(stringLocusEffect,"0\t0")
  } else {
    stringLocusEffect=paste(stringLocusEffect,round(totEnvirVar_FCRITBB*alleleEffectMultiplCoeffFCRITBB),"\t0")
  }
  cat (stringLocusEffect)                        
  
  stringLocusEffect = paste0("\n", demographic_parameters$idSp[2], "\tg1max\t1\t{")
  for (i in 1:nbLocus_g1max)
  {
    if(herit_g1max>0) {
      all1=(target_mean_g1max/nbLocus_g1max-pdgInventoryOptions$effect_g1maxdraw[i])/2
      all1=round(all1*alleleEffectMultiplCoeffg1max)
      all2=(target_mean_g1max/nbLocus_g1max+pdgInventoryOptions$effect_g1maxdraw[i])/2
      all2=round(all2*alleleEffectMultiplCoeffg1max)
      #hom1=-effect_g1max[i]
      #hom2=+effect_g1max[i]
    } else {
      all1=0 
      all2=0
    }
    
    stringLocusEffect=paste(stringLocusEffect,"[",(i+10),",",all1,",",all2,"]",sep="")
    if (i<nbLocus_g1max) stringLocusEffect=paste(stringLocusEffect,";",sep="")
  }
  stringLocusEffect=paste(stringLocusEffect,"}\t{}\t{}\t",herit_g1max,"\t",sep="")
  if (totEnvirVar_g1max==0) {
    stringLocusEffect=paste(stringLocusEffect,"0\t0")
  } else {
    stringLocusEffect=paste(stringLocusEffect,round(totEnvirVar_g1max*alleleEffectMultiplCoeffg1max),"\t0")
  }
  cat (stringLocusEffect)  
  
  
  
  
  # 4 PHI
  
  cat ("\n\n#Phi\n#speciesId	speciesPhi	DefaultPhi")
  cat("\n", demographic_parameters$idSp[1],"	{}	0", sep = "")
  cat("\n", demographic_parameters$idSp[2],"	{}	0", sep = "")
  
  # 4 CELLS
  ## Cells / cellules
  cat("\n\n## CELLS\n#cID\tclign\tccol\tsolHeight\tstoneContent\twfc\twilt\tpropMacro\tpropMacroDeep\tbulk\tSOLCLAYtop\tSOLCLAYsol\tSOLFINtop\tSOLFINsol\tSOLSANDtop\tSOLSANDsol\tdeepSoilDepth\tstoneContentDeep\tprac\tpracDeep	")
  
  for (i in 1:length(cellsTable$IdCell))
  {
    cat("\n")
    stringCell=paste(cellsTable$IdCell[i], cellsTable$clign[i], cellsTable$ccol[i], cellsTable$soilHeight[i], cellsTable$stoneContent[i], cellsTable$wfc[i], cellsTable$wilt[i], cellsTable$propMacro[i], cellsTable$propMacroDeep[i], cellsTable$bulk[i], cellsTable$SOLCLAYtop[i], cellsTable$SOLCLAYsol[i],
                     cellsTable$SOLFINtop[i], cellsTable$SOLFINsol[i], cellsTable$SOLSANDtop[i], cellsTable$SOLSANDsol[i], cellsTable$deepSoilDepth[i], cellsTable$stoneContentDeep[i], cellsTable$prac[i], cellsTable$pracDeep[i], sep="\t")
    cat(stringCell)
  }
  
  # 4 TREES
  
  ## Trees /individus
  cat("\n\n## TREES  (Individual-TREES)	\n#ID\tpop\tcID\tspeciesId\trx\try\tDBH\theight\thcb\tAge\tLAIinit\tnucDNA\tmCytDNA\tpCytDNA\tcreationDate\tmID\tpID")
  nbLocusTot= nbLocus_FCRITBB+ nbLocus_g1max+nbMsat+nbSNP
  
  i_1 = 0
  i_2 = 0
  
  # check for NAs
  NA_by_cols = apply(is.na(treesTable), FUN = sum, MARGIN = 2)
  if(sum(NA_by_cols) > 0){
    NA_columns_vec = names(which(NA_by_cols > 0))
    NA_columns = paste0(NA_columns_vec, collapse = " ")
    
    isNA_columnsImportantVector = NA_columns_vec %in% c("idTree", "pop", "sp", "cell_IDs", "xnorm", "ynorm", "dbh", "height", "hcb", "age", "LAIinit")
    
    if(TRUE %in% isNA_columnsImportantVector){
      while(sink.number() > 0){sink()}
      stop("NAs in trees for inventory: ", inventoryName, ". For columns: ", NA_columns)
    }else{
      # do nothing
    }
  }
  
  
  # rounding value
  roundingIndex = 3
  treesTable$xnorm = floor(treesTable$xnorm * 10**roundingIndex) / 10**roundingIndex
  treesTable$ynorm = floor(treesTable$ynorm * 10**roundingIndex) / 10**roundingIndex
  treesTable$dbh = round(treesTable$dbh, roundingIndex)
  treesTable$height = round(treesTable$height, roundingIndex)
  treesTable$hcb = round(treesTable$hcb, roundingIndex)
  treesTable$age = round(treesTable$age, roundingIndex)
  
  for(idTree in treesTable$idTree){
    whichTree = which(treesTable$idTree == idTree)
    tree = treesTable[whichTree, ]
    # check, errors..
    if(tree$hcb >= tree$height ){
      stop(paste0( "In ", inventoryName, ": height crown base is greater or equal to tree height" ))
    }
    
    if("dbh2013" %in% colnames(treesTable)){
      if(!is.na(tree$dbh2013) & tree$dbh > tree$dbh2013 ){
        stop(paste0( "In ", inventoryName, ": the ancient dbh estimation is greater thant actuel dbh" ))
      }
    }
    
    if("height2013" %in% colnames(treesTable)){
      if(!is.na(tree$height2013) & tree$height > tree$height2013 ){
        stop(paste0( "In ", inventoryName, ": the ancient height estimation is greater thant actuel height" ))
      }
    }
    
    if("age2013" %in% colnames(treesTable)){
      if(!is.na(tree$age2013) & tree$age2013 != -1 &  tree$age > tree$age2013 ){
        stop(paste0( "In ", inventoryName, ": the ancient age estimation is greater thant actuel age" ))
      }
    }
    
    if("hcb2013" %in% colnames(treesTable)){
      if(!is.na(tree$hcb2013) & tree$hcb > tree$hcb2013 ){
        stop(paste0( "In ", inventoryName, ": the ancient hcb estimation is greater thant actuel hcb" ))
      }
    }
    
    # for(var in colnames(tree)){
    #   val = tree[[var]]
    #   if(is.na(val)){
    #     stop(paste0( "In ", inventoryName, ": ", var, " is Na" ))
    #   }
    #   if(val < 0){
    #     stop(paste0( "In ", inventoryName, ": ", var, " is less than zero" ))
    #   }
    # }
    
    cat("\n")
    stringIndiv=paste(tree$idTree, tree$pop, tree$cell_IDs, tree$sp, tree$xnorm, tree$ynorm,
                      tree$dbh, tree$height, tree$hcb, tree$age, tree$LAIinit, "{", sep="\t")
    
    allele1 = NULL
    allele2 = NULL
    
    if(tree$sp == demographic_parameters$idSp[1]) { # if species 1
      i_1 = i_1 + 1 
      allele1 = allele1_1[i_1, ]
      allele2 = allele2_1[i_1, ]
    }else if(tree$sp == demographic_parameters$idSp[2]) { # if species 2
      i_2 = i_2 + 1 
      allele1 = allele1_2[i_2, ]
      allele2 = allele2_2[i_2, ]
    }
    
    # warning : dont use i
    stringIndiv=paste(stringIndiv, allele1[1],";",allele2[1], sep=" ")
    for (j in 2:nbLocusTot) {
      stringIndiv=paste(stringIndiv, allele1[j], allele2[j], sep="; ")
    }
    stringIndiv=paste(stringIndiv, "}", sep="")
    stringIndiv=paste(stringIndiv, "{}", "{}", "-1", "-1",	"-1", sep="\t")
    cat(stringIndiv)
  }
  
  # redirect output to console
  while(sink.number() > 0){sink()}
  
  cat("\nPDG inventory vas generated.\nFolder: ", inventoryFolder, "\nFile: ", inventoryName, "\n")
  
} # end writeInventoryPDG



# Write the inventory for CASTANEA from a demographically-regular and monospecific plot of tree cell
# FmCell line is retrieved from cellsTable and treesTable
# only the first line of these tables is used
writeInventoryCASTANEAFromRegdemoMonosp = function(inventoryFolder, inventoryName, 
                                                   pdgInventoryOptions,
                                                   castaneaSpeciesCode,
                                                   description = "",
                                                   cellsTable, treesTable, standArea, fmCellWidth = 80,
                                                   initYear,
                                                   outputInConsole, LAI){
  
  CASTANEA_simulation_parameters = pdgInventoryOptions$CASTANEA_simulation_parameters
  plotLocationParameters = pdgInventoryOptions$plotLocationParameters
  
  
  # make sure that even if error occurs, output will be redirected to console
  on.exit(while(sink.number() > 0){sink()}, add = TRUE, after = FALSE)
  
  # Open writing
  if(!outputInConsole)
    sink(paste0(inventoryFolder, inventoryName))
  
  cat(description)
  
  cat("\n\n\t#<General>");
  
  # 4 Write PARAMETERS
  
  
  cat("\n\n\t#SIMULATION parameters");
  cat(paste("\nyear= ", initYear, sep=" "));
  
  
  cat("\n\n\t#PLOT parameters\n")
  cat(paste("\nFmCellWidth =", fmCellWidth, sep=" "))
  cat("\nlatitude =", plotLocationParameters$latitude, sep=" ")
  cat("\nlongitude =", plotLocationParameters$longitude, sep=" ")
  
  
  CASTANEA_simulation_parameters = CASTANEA_simulation_parameters[-which(colnames(CASTANEA_simulation_parameters) == 'Ca')]
  CASTANEA_simulation_parameters$LAImode = "LAI_FIXED"
  CASTANEA_simulation_parameters$typeOfVegetation = "TYPE_VEG_AVG_TREE"
  
  cat("\n\n\t#CASTANEA parameters")
  cat("\n")
  cat( paste0( paste0(colnames(CASTANEA_simulation_parameters), " = ", c(CASTANEA_simulation_parameters[1, ])), collapse = "\n" ) )
  
  
  # 4 OPTIONAL VARIABLES 
  cat("\n\n#possible optional variables:  nitrogen LMA nc coefBeta MRN gfac slopePotGs Nstress Tsum woodStop GBV CRBV tronviv PotWood coefrac rootshoot aGF dateDeb rateOfSeedcell costSeedcell meanHeight")
  optionalVariables = "meanHeight"
  
  if(optionalVariables == ""){
    cat(paste0("\n", "optionalVariables = ", "-"))
  }else{
    cat(paste0("\n", "optionalVariables = ", optionalVariables))
  }
  
  # 4 FMCELL
  nha = dim(treesTable)[1] / (standArea/10000) # tree / ha
  vha = 0
  clumping = 0
  
  # stop("# TODO : vha, clumping")
  oneCell = cellsTable[1, ]
  oneTree = treesTable[1, ]
  
  
  cat("\n\n\t#Flux-model cells")
  cat("\n#Type\tName\tSpecies\txc\tyc\tzc\tsoilHeight\tstone\twfc\twilt\tpropMacro\tbulk\tSOLCLAYtop\tSOLCLAYsol\tSOLFINtop\tSOLFINsol\tSOLSANDtop\tSOLSANDsol\tprac\tpracDeep\tdeepSoilDep\tstoneDeep\tMacroDeep\tdbh\tNha\tVha\tage\tclumping\tLAI\toptional variable")
  
  roundingIndex = 3
  
  stringLine = paste("FMC", "1", castaneaSpeciesCode, "20", "20", plotLocationParameters$referenceAltitude, 
                     oneCell$soilHeight, oneCell$stoneContent, oneCell$wfc, oneCell$wilt,  
                     oneCell$propMacro, oneCell$bulk, oneCell$SOLCLAYtop, oneCell$SOLCLAYsol,  
                     oneCell$SOLFINtop, oneCell$SOLFINsol, oneCell$SOLSANDtop, oneCell$SOLSANDsol,  
                     oneCell$prac, oneCell$pracDeep, oneCell$deepSoilDepth, oneCell$stoneContentDeep, oneCell$propMacroDeep,
                     round(oneTree$dbh, roundingIndex), round(nha, roundingIndex), round(vha, roundingIndex), oneTree$age, clumping, round(LAI, 3), 
                     round(oneTree$height, roundingIndex), # OPTIONAL VARIABLES (mean height here)
                     sep = "\t")
  cat(paste0('\n', stringLine))
  
  
  # redirect output to console
  while(sink.number() > 0){sink()}
  
  cat("\nGMAP CASTANEA inventory vas generated.\nFolder : ", inventoryFolder, "\nFile : ", inventoryName, "\n")
  
} # end writeInventoryCASTANEAFromRegdemoMonosp