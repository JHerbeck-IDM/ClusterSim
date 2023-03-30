# read in line list
lineList <- read.csv("line_list.csv")

# Make a parenthetic pair with its source individual.
# When the source individual is also a recipient, find their source in the line list), 
# make a new parentheses with the first term being the newly identified source 
# and the second term being your original parentheses. 
# Do this iteratively for every source.

# iterate over rows
for (i in 1:nrow(lineList)) {
  
  # get source and recipient for current row
  source <- lineList$source[i]
  recipient <- lineList$recipient[i]
  
  # create new parenthetic pair with source and recipient
  Pair <- paste0("(", source, ",", recipient, ")")
  
  # check if source is also a recipient
  if (source %in% lineList$recipient) {
    # get that recipient's source the sources list
    newSource <- lineList$source[lineList$recipient == source]
    newRecipient <- source
    
    # make a new parenthetic pair
    newPair <- paste0("(", newSource, ",", newRecipient, ")")
    
    # replace recipient ID with the source's parenthetic pair
    newPair <- gsub(source, newPair, Pair)
    
    Newick_string <- newPair
  }
  
}

# print final parenthetic pair
#cat()
