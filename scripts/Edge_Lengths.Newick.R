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
  if (recipient_id %in% df$source) {
    # subset to rows where the recipient is a source
    recipient_source_rows <- df[df$source == recipient_id, ]
    
    # find the row with the highest infection time
    highest_infection_time_row <- recipient_source_rows[which.max(recipient_source_rows$infectionTime), ]
    
    # calculate the leaf edge length using the infection time from the highest row
    leaf_edge_length <- recipient_rows$sampleTime - highest_infection_time_row$infectionTime
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




# Tree with dummy edge lengths
test <- "((1:1, ((2:1, ((4:2, 7:2):1, (6:1, 8:1):1):1):1, ((3:1, 5:1):1):1):1):1);"

tree <- ape::read.tree(text = test)
plot(tree)


