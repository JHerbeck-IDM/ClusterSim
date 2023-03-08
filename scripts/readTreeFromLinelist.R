# Translated from Greg's Python code:
# https://github.com/InstituteforDiseaseModeling/phyloModels/blob/master/phylomodels/trees/generate_treeFromFile.py
library(ape)
library(TreeTools)

# Read in the transmission line list ("population_summary")
# df with headers "infection_source", "ID", and "infection_day", 

linelist <- population_summary[ , c(1, 7, 8)]
#linelist$infection_day <- as.Date(transmissions$infection_day, origin = "1970-01-01")

read_treeFromLineList <- function(linelist, ...) {
  
  # Default arguments
  
  ID <- "linelist$ID"
  infectorID <- "linelist$infection_source"
  infectTime <- "linelist$infection_day"
  sampleTime <- NULL
  features <- NULL
  branchLengths <- TRUE
  
  
  # Initial empty list of trees
  
  trees <- list()
  
  
# Find IDs that infect but are never infected (the roots)
  
roots <- linelist$ID[1:samplesize]
  
# Loop through all the roots and recursively build the trees
  
  for (i in roots) {
    
    # Create a tree for each root/seed/importation
    #temp_tree <- TreeTools::ZeroTaxonTree()
    temp_tree <- TreeTools::SingleTaxonTree(label = i)
    
    # Call the recursive add function
    times <- add_childrenFromLineList(temp_tree, linelist, ID, infectorID,
                                      infectTime, sampleTime, features)
    times <- unlist(times)
    # We do not have the time of infection for the seed cases/importations
    # so here we take a list of branch lengths and use its mean to estimate
    # the time of infection for the seed cases. We are flooring this value
    # at zero.
    temp_tree$add_feature("time", max(0, temp_tree$children[[1]]$time - mean(times)))
    temp_tree$add_feature("infectTime", max(0, temp_tree$children[[1]]$time - mean(times)))
    if (length(temp_tree$children) == 1) {
      node <- temp_tree
      temp_tree <- node$children[[1]]
      node$delete()
    }
    trees <- c(trees, temp_tree)
  }
  
  if (branchLengths) {
    calculate_branch_lengths(trees, inplace = TRUE)
  }
  
  return(trees)
}