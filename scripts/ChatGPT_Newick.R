# load required packages
library(ape)
library(dplyr)
library(igraph)

# read in transmission line list
line_list <- read.csv("line_list.csv")

# function to calculate leaf edge length
calculate_leaf_edge_length <- function(df, recipient_id) {
  # subset the data frame to rows where the recipient is the specified ID
  recipient_rows <- df[df$recipient == recipient_id, ]
  
  # check if the recipient is also a source in any other rows
  if (recipient_id %in% recipient_rows$source) {
    # subset to rows where the recipient is a source
    recipient_source_rows <- recipient_rows[recipient_rows$source == recipient_id, ]
    
    # find the row with the highest infection time
    highest_infection_time_row <- recipient_source_rows[which.max(recipient_source_rows$infectionTime), ]
    
    # calculate the leaf edge length using the infection time from the highest row
    leaf_edge_length <- highest_infection_time_row$sampleTime - highest_infection_time_row$infectionTime
  } else {
    # calculate the leaf edge length using the infection time from the current row
    leaf_edge_length <- recipient_rows$sampleTime - recipient_rows$infectionTime
  }
  
  # return the leaf edge length
  return(leaf_edge_length)
}

# create a new variable "leafEdgeLength" using the helper function
line_list$leafEdgeLength <- sapply(line_list$recipient, function(recipient_id) calculate_leaf_edge_length(line_list, recipient_id))

# print the resulting data frame
print(line_list)


test <- "((1:2, ((2:1, (3:1,5:1):0.5):1, ((4:2, (6:1, 8:1):1):1, 7:2):1):1):1);"
test <- ((1:2, ((2:1, (3:1,5:1):0.5):1, ((4:2, (6:1, 8:1):1):1, 7:2):1):1):1);"

tree <- ape::read.tree(text = test)
plot(tree)

transTreeFromLineList <- function(line_list) {
  
  # Convert to unweighted transmission tree graph to get edge list and node list
  infectionGraph <- graph_from_data_frame(line_list, directed = TRUE)
  
  # Add sampling time data to node/vertex data
  V(infectionGraph)$sampleTime <- c(NA,line_list$sampleTime)
  
  # Add infection time to node/vertex data
  V(infectionGraph)$timeInfected <- c(NA,line_list$infectionTime)
  
  # Identify root (first infection)
  root <- V(infectionGraph)[degree(infectionGraph, mode="in") == 0]
  
  for (i in V(infectionGraph)) {
    node <- V(infectionGraph)[i]
    childNode <- line_list$recipient[line_list$source == node]
    edge.length[i] <- line_list$timeInfected[line_list$recipient == childNode]
  }
  
  
for(i in as.numeric(V(infectionGraph)$name)){
    from <- E(infectionGraph)[.from(i)]}
  
  # Find distances between nodes / calculate edge lengths
  for (i in V(infectionGraph)) 
    locParent <- match(as.character(E(infectionGraph)$from), names)
  # Create a variable 'locParent' that matches the "from" node of each edge in 
  # the infectionGraph object to its corresponding index in the 'names' vector. 
  # This essentially finds the index of the parent node for each edge in the graph.
  locChild <- match(as.character(E(infectionGraph)$to), node_names)
  weights <- timeInfected[locChild] - timeInfected[locParent]
  
  # Time-interval weighted directed graph transmission tree
  edgeTable <- data.frame(from = E(infectionGraph)$from, to = E(infectionGraph)$to, weight = weights)
  nodeTable <- data.frame(name = V(infectionGraph)$name, timeInfected = timeInfected)
  
  transmissionTree <- graph_from_data_frame(edgeTable, directed = TRUE, vertices = nodeTable)
  
  return(transmissionTree)
}
