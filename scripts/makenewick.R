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

