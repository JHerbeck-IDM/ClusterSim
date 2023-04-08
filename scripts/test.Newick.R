# read in line list
lineList <- read.csv("line_list.csv")

# Make a parenthetic pair with its source individual.
# When the source individual is also a recipient, find their source in the line list, 
# Make a new parentheses with the first term being the newly identified source 
# and the second term being your original parentheses. 
# Do this iteratively for every source.

lineList <- lineList[order(lineList$infectionTime), ]

tree <- NULL

# iterate over rows
for (i in nrow(lineList):1) {
  
  # get source and recipient for current row
  source <- lineList$source[i]
  recipient <- lineList$recipient[i]
  
  # create new parenthetic pair with source and recipient
  Pair <- paste0("(", source, ",", recipient, ")")
  
  tree <- c(tree, Pair)
  
  # If new source is not already in tree
  if ( (lineList$source[i] %in% tree,
       
       Pair <- paste0("(", source, ",", recipient, ")")
       )
    
    # create new parenthetic pair with source and recipient
    Pair <- paste0("(", source, ",", recipient, ")")
  
    
    
  # add it to the tree
  tree <- c(tree, Pair)
  
  # check if source is also a recipient
  if (source %in% lineList$recipient) {
    # get that recipient's source from the sources list
    newSource <- lineList$source[lineList$recipient == source]
    newRecipient <- source
    
    # make a new parenthetic pair
    newPair <- paste0("(", newSource, ",", newRecipient, ")")
    
    # replace recipient ID with the source's parenthetic pair
    newPair <- gsub(source, newPair, Pair)
  }
  
    tree <- c(tree, newPair)
  
}

test = "(1, ((2, ((4, 7), (6, 8))), (3, 5)));"
tree <- ape::read.tree(text = test)
plot(tree)


# "((4,6),8)" "((2,4),7)" "((2,4),6)" "((2,3),5)" "((1,2),4)" "((1,2),3)" "((0,1),2)" "((0,1),2)"
# ((1, ((2, ((4, 7), (6, 8))), ((3, 5)))))
