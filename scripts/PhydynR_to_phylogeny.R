# Install the phydynR package
install.packages("phydynR")

library(phydynR)

# Read in the transmission line list data
# df with headers "infection_source", "ID", and "infection_day", 
# "infection_day" variable is in a format that can be parsed by R's as.Date() function.

transmissions <- population_summary[ , c(1, 7, 8)]
transmissions$infection_day <- as.Date(transmissions$infection_day, origin = "1970-01-01")

# Convert the transmission data to a transmission tree object
trans_tree <- from_treetable(transmissions, 
                             transmission_from = "infection_source", 
                             transmission_to = "ID", 
                             time_column = "infection_day")

# Convert the transmission tree to a phylogenetic tree
phylo_tree <- from_transmission(trans_tree)

# Plot the phylogenetic tree
plot(phylo_tree)
