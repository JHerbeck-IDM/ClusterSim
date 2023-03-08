add_childrenFromLineList <- function(node, lineList, ID, infectorID, infectTime,
                                     sampleTime=NULL, features=NULL) {
  
  # Use the line list to find all the people infected by the current tree node
  # and add them as children.
  #
  # Args:
  #   node (TreeNode): The node of the tree that contains the current infector.
  #   lineList (dataframe): The list line that includes who infected whom 
  #                         (edge list).
  #   ID (str): The name of the column containing the ID or name of who was
  #             infected.
  #   infectorID (str): The name of the column containing the ID or name of the
  #                      source of the infection.
  #   infectTime (str): The name of the column containing the time of infection.
  #   sampleTime (str): Optional (default NULL). The name of the column containing
  #                     the time of sampling (if sampled).
  #   features (list): Optional (default NULL). List of features to add to the tree
  #                     from the line list. An empty list (NULL) will add all features.
  #
  # Returns:
  #   times (list): This function directly modifies the tree through the node that
  #                 is passed. However it also returns a list of branch length to help
  #                 estimate the time of infection for the root.
  
  #node_name <- node$name
  node_name <- temp_tree$tip.label[1]
  children <- linelist$ID[linelist$infection_source == node_name]
  times <- list()
  for (child_name in children) {
    child <- node$add_child(name = child_name)
    idx <- linelist$ID == child_name
    child$add_feature('time', linelist[idx, infectTime][1])
    child$add_feature('infectTime', linelist[idx, infectTime][1])
    if (!is.null(sampleTime)) {
      time <- linelist[idx, sampleTime][1]
      child$add_feature('sampleTime', time)
    }
    if (!is.null(features)) {
      for (feature in features) {
        child$add_feature(feature, lineList[idx, feature][1])
      }
    }
    times <- c(times, child$time - ifelse(is.null(node$time), child$time, node$time))
    temp_times <- add_childrenFromLineList(child, linelist, ID, infectorID,
                                           infectTime, sampleTime, features)
    times <- c(times, temp_times)
  }
  return(times)
}
