# Load the required library
library(dplyr)

# Load the transmission line list as a data frame
transmission_list <- read.csv("transmission_list.csv")

# Create a list of unique individuals (infected and infectors)
unique_individuals <- sort(unique(c(transmission_list$infected, transmission_list$infector)))

# Create an empty transmission matrix with dimensions equal to the number of unique individuals
transmission_matrix <- matrix(0, nrow = length(unique_individuals), ncol = length(unique_individuals))

# Loop through each row of the transmission list and update the transmission matrix
for (i in 1:nrow(transmission_list)) {
  infected <- transmission_list$infected[i]
  infector <- transmission_list$infector[i]
  transmission_matrix[which(unique_individuals == infected), which(unique_individuals == infector)] <- 1
}

# Set the row and column names of the transmission matrix to match the unique individuals
rownames(transmission_matrix) <- unique_individuals
colnames(transmission_matrix) <- unique_individuals

# Print the transmission matrix
print(transmission_matrix)





# Create a matrix of transmission events with rows as infected and columns as infectors
transmission_matrix <- matrix(c(
  0, 1, 0, 0, 0,
  1, 0, 0, 0, 0,
  0, 0, 0, 1, 1,
  0, 0, 1, 0, 0,
  0, 0, 1, 0, 0
), nrow = 5, ncol = 5, byrow = TRUE)

# Create a vector of sample dates for each infected individual
sample_dates <- c(1, 2, 2, 3, 4)

# Convert the transmission matrix to a phylogenetic tree
tree <- as.phylo(transmission_matrix, tip.label = 1:length(sample_dates), 
                 edge.length = sample_dates)

# Plot the phylogenetic tree
plot(tree)
