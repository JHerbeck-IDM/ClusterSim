

lineList <- read.csv("line_list.csv")

#lineList <- tribble(
#  ~source, ~dest,
#  0,1,
#  1,2,
#  2,3,
#  2,4,
#  3,5,
#  4,6,
#  4,7,
#  6,8)

#lineList <- lineList[order(lineList$infectionTime), ]

recurhelper <- function(source) {
  # assumes sources are sorted
  if (length(source) == 1) {
    recur2(source)
  } else if (length(source) > 1) {
    paste0("(", recur2(source[1]), ",", recurhelper(source[-1]), ")")
  }
}

recur2 <- function(source) {
  if (length(source) == 1) {
    if (sum(lineList$source == source) == 0) {
      as.character(source)
    } else {
      dests <- lineList$recipient[lineList$source == source]
      if (length(dests) <= 1) {
        paste0("(", source, ",", recur2(dests), ")")
      } else {
        sorteddests <-
          rev(dests[(order(lineList$infectionTime[match(dests, lineList$dest)]))]) # sort by recency when >1 dest. Ties are not resolved.
        paste0(
          "((",
          source,
          ",",
          recur2(sorteddests[1]),
          "),",
          recurhelper(sorteddests[-1]),
          ")"
        )
      }
    }
  } else {
    print("ERROR")
    NULL # error
  }
}
headnode <- lineList$recipient[lineList$source == 0]

test <- print(paste0(recur2(headnode), ";"))

tree <- ape::read.tree(text = test)

