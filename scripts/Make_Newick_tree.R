##### Make line list

str(population_summary)

# Make a Newick tree from each seed infection
# "headnode" is the ID of a single seed infection

#headnode <- sampleLineList$recipient[sampleLineList$source==0]
headnode <- population_summary$recipient[population_summary$source == 0]
headnode <- population_summary$recipient[2]

test <- print(paste0(makenewickstring(headnode, population_summary),";"))

tree <- ape::read.tree(text = test)
plot(tree)
