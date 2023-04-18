# makenewick.R
# outputs a Newick representation of a tree based on a line list
# Dennis Chao
# for Josh Herbeck
# April 2023

lineList <- read.csv("CopyOfline_list.csv")
lineList <- read.csv("line_list.csv")

lineList <- order(lineList, decreasing = TRUE)

makenewickstring <- function(source, lineList) {
  # takes the head node(s) as an argument ("source")
  if (length(source) == 1) {
    recipients <- lineList$recipient[lineList$source == source]
    if (length(recipients) == 0) {
      # leaf
      paste0(as.character(source),
             ":",
             (lineList$sampleTime - lineList$infectionTime)[source == lineList$recipient])
    } else if (length(recipients) == 1) {
      # (a,b)
      paste0(
        "(",
        source,
        ":",
        min((
          lineList$sampleTime[source == lineList$recipient] - lineList$infectionTime
        )[source == lineList$source]),
        ",",
        makenewickstring(recipients, lineList),
        "):",
        (lineList$infectionTime[recipients == lineList$recipient] -
           lineList$infectionTime[source == lineList$recipient])
      )
    } else {
      # ((a,b), c)
      sortedrecipients <-
        rev(recipients[(order(lineList$infectionTime[match(recipients, lineList$recipient)]))]) # sort by recency when >1 recipient. Ties are not resolved.
      paste0(
        "((",
        source,
        ":",
        min(lineList$sampleTime[source == lineList$recipient] - lineList$infectionTime[source ==
                                                                                         lineList$source]),
        ",",
        makenewickstring(sortedrecipients[1], lineList),
        "):",
        (lineList$infectionTime[recipients[2] == lineList$recipient] -
           lineList$infectionTime[recipients[1] == lineList$recipient]),
        ",",
        makenewickstring(sortedrecipients[-1], lineList),
        "):",
        (lineList$infectionTime[recipients[1] == lineList$recipient] -
           lineList$infectionTime[source == lineList$recipient])
      )
    }
  } else if (length(source) > 1) {
    # recursively converts a vector of sources into binary trees. E.g., "a b c d" -> (((a,b),c),d)
    # assumes sources are sorted by time
    paste0(
      "(",
      makenewickstring(source[1], lineList),
      ",",
      makenewickstring(source[-1], lineList),
      "):?",
      (lineList$infectionTime[source[1] == lineList$recipient] -
         lineList$infectionTime[source[2] == lineList$recipient])
    ) # not sure if branch length is correct
  } else {
    # no source passed in?
    print("ERROR")
    NULL # error
  }
}

headnode <- sampleLineList$recipient[lineList$source == 0]

test <- print(paste0(makenewickstring(headnode, lineList), ";"))
test <- "(1:1,((2:-4,9:2):1,(((4:1,7:1):1,(6:0,8:1):2):2,(3:0,5:1):2):10):2):1;"
tree <- ape::read.tree(text = test)
plot(tree)
