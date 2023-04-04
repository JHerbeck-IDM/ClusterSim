lineList <- read.csv("line_list.csv")


dag <- igraph::graph_from_data_frame(lineList, directed = TRUE)
