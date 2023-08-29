# makenewick.R
# outputs a Newick representation of a tree based on a line list
# Dennis Chao, for Josh Herbeck
# May 2023

makenewickstring <- function(source, lineList) {
    source <- source[!is.na(source)]
    # takes the head node(s) as an argument ("source")
    if (length(source)==1) {
        recipients <- lineList$recipient[!is.na(lineList$source) & lineList$source==source]
        if (length(recipients)==0) { # leaf
            paste0(as.character(source),":",(lineList$sampleTime-lineList$infectionTime)[!is.na(lineList$recipient) & lineList$recipient==source])
        } else if (length(recipients)==1) {  # (a,b)
            paste0("(",source,":",min((lineList$sampleTime[!is.na(lineList$recipient) & lineList$recipient==source]-lineList$infectionTime)[!is.na(lineList$source) & lineList$source==source]),",",
                   makenewickstring(recipients, lineList),"):",
                   (lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==recipients] -
                    lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==source]))
        } else { # ((a,b), c)
            sortedrecipients <- rev(recipients[(order(lineList$infectionTime[match(recipients, lineList$recipient)]))]) # sort by recency when >1 recipient. Ties are not resolved.
            paste0("((",source,":",min(lineList$sampleTime[!is.na(lineList$recipient) & lineList$recipient==source]-lineList$infectionTime[!is.na(lineList$source) & lineList$source==source]),",",makenewickstring(sortedrecipients[1], lineList),"):",
            (lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==recipients[2]] -
             lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==recipients[1]]),",",
            makenewickstring(sortedrecipients[-1], lineList),"):",
            (lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==recipients[1]] -
             lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==source]))
        }
    } else if (length(source)>1) {
        # recursively converts a vector of sources into binary trees. E.g., "a b c d" -> (((a,b),c),d)
        # assumes sources are sorted by time
        paste0("(",makenewickstring(source[1], lineList),",",makenewickstring(source[-1], lineList),"):",
        (lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==source[1]] -
         lineList$infectionTime[!is.na(lineList$recipient) & lineList$recipient==source[2]]))
    } else { # no source passed in?
        print("ERROR")
        NULL # error
    }
}

# working example
library(tidyr) # this is needed just to make the sample data
sampleLineList <- tribble(
    ~source, ~recipient, ~infectionTime, ~sampleTime, 
    0,1,1,3,
    1,2,2,5,
    2,3,4,6,
    2,10,4,12,
    2,4,5,9,
    3,5,6,7,
    4,6,7,9,
    4,7,8,9,
    6,8,9,10)

headnode <- sampleLineList$recipient[sampleLineList$source==0]
print(paste0(makenewickstring(headnode, sampleLineList),";"))

# second example with three-way branch
sampleLineList <- tribble(
    ~source, ~recipient, ~infectionTime, ~sampleTime, 
    0,1,1,3,
    1,2,2,5,
    2,3,4,6,
    2,4,5,9,
    2,9,9,11,
    3,5,6,7,
    4,6,7,9,
    4,7,8,9,
    6,8,9,10)

headnode <- sampleLineList$recipient[sampleLineList$source==0]
print(paste0(makenewickstring(headnode, sampleLineList),";"))

# large example
# ulimit -s 16384 # enlarge stack limit to 16 megs
# options(expressions=10000)
bigdat <- read.csv("population_summary.csv")
for (headnode in bigdat$X[bigdat$source==0]) {
    print(paste("source",headnode))
    print(paste0(makenewickstring(headnode, bigdat),";"))
}
