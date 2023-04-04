# read in line list
lineList <- read.csv("line_list.csv")

# Gather info from the line list. Each line tells us about two nodes: 
# an internal node for the transmission event and a leaf node for the recipient.

all.the.leaf.nodes <- t( data.frame( "root" = c( name = 0, time = 0, parent = NA ) ) );
# "time" is $sampleTime
all.the.internal.nodes <- t( data.frame( "dummy" = c( source = 0, recipient = 0, time = 0, parent = NA, leftChild = NA, rightChild = NA ) ) );

add.nodes.from.line  <- function ( line ) {
    # source, recipient, infectionTime, sampleTime
    source  <- line[ "source" ];
    recipient  <- line[ "recipient" ];
    infectionTime  <- line[ "infectionTime" ];
    sampleTime  <- line[ "sampleTime" ];

    new.leaf.node  <- c( name = recipient, time = sampleTime, parent = NA );
    new.leaf.node.name  <- as.character( recipient );
    ## NOTE this is inefficient and could be optimized later:
    all.the.leaf.nodes  <<- rbind( all.the.leaf.nodes, new.leaf.node );
    rownames( all.the.leaf.nodes )[ nrow( all.the.leaf.nodes ) ]  <<- new.leaf.node.name;

    new.internal.node  <- c( source = source, recipient = recipient, time = infectionTime, parent = NA, leftChild = NA, rightChild = NA );
    new.internal.node.name  <- paste( source, "->", recipient, sep = "" );
    ## NOTE this is inefficient and could be optimized later:
    all.the.internal.nodes  <<- rbind( all.the.internal.nodes, new.internal.node );
    rownames( all.the.internal.nodes )[ nrow( all.the.internal.nodes ) ]  <<- new.internal.node.name;

    return( NULL );
} # add.nodes.from.line (..)

# ?
result.ignored <- apply( lineList, 1, add.nodes.from.line );

# all.the.leaf.nodes are not sorted in time already.
all.the.leaf.nodes <- all.the.leaf.nodes[ order( all.the.leaf.nodes[ , "time" ] ), , drop = FALSE ];

## Now figure out the parent of each leaf, and build the tree recursively leaf-first.
# start at a leaf at the end of time (one with the max observed time value)
the.tree <- t( data.frame( dummy = c( name = "dummy", parent = 0, leftChild = 0, rightChild = 0, time = 0 ) ) );
the.tree[ , c( "parent", "leftChild", "rightChild" ) ]  <- NA;

build.nodes.from.leaf <- function ( the.current.leaf.node.row.i, debug = FALSE ) {
    the.current.leaf.node.row  <- all.the.leaf.nodes[ the.current.leaf.node.row.i, ];

    leaf.name  <- the.current.leaf.node.row[ "name" ];

    # First, get all internal nodes that mention this leaf node.
    internal.nodes.mentioning.this.node.indices <-
        which( 
          ( all.the.internal.nodes[ , "source" ] == leaf.name ) |
          ( all.the.internal.nodes[ , "recipient" ] == leaf.name )
        );
    
    # Next, sort them.
    internal.nodes.mentioning.this.node <-
        all.the.internal.nodes[ internal.nodes.mentioning.this.node.indices, , drop = FALSE ];
    internal.nodes.mentioning.this.node  <-
        internal.nodes.mentioning.this.node[ order( internal.nodes.mentioning.this.node[ , "time" ], decreasing = TRUE ), , drop = FALSE ];

    # For the leaf, we will add this new node:
    new.leaf.node.for.tree <- the.tree[ "dummy", ];

    new.leaf.node.for.tree[ "name" ]  <- the.current.leaf.node.row[ "name" ];
    new.leaf.node.for.tree[ "time" ]  <- the.current.leaf.node.row[ "time" ];

    # Set the parent of the new leaf node to the latest-time internal node mentioning it.
    new.leaf.node.for.tree[ "parent" ]  <- 
        rownames( internal.nodes.mentioning.this.node )[ 1 ];
if( debug ) {
    all.the.leaf.nodes[ the.current.leaf.node.row.i, "parent" ] <-
        rownames( internal.nodes.mentioning.this.node )[ 1 ];
    the.tree <- rbind( the.tree, new.leaf.node.for.tree );
} else {
    all.the.leaf.nodes[ the.current.leaf.node.row.i, "parent" ] <<-
        rownames( internal.nodes.mentioning.this.node )[ 1 ];

    the.tree <<- rbind( the.tree, new.leaf.node.for.tree );
}

    # Start from the new leaf we just added.
    current.node.i <- nrow( the.tree );

    # For each of those internal nodes, if it is not on the tree, put it there. Either way, add a child to it.
    for( internal.node.in.reverse.order.i in 1:nrow( internal.nodes.mentioning.this.node ) ) {
        # is it already on the tree? If both children are empty, it's not on the tree.
        if( is.na( internal.nodes.mentioning.this.node[ internal.node.in.reverse.order.i, "leftChild" ] ) && is.na( internal.nodes.mentioning.this.node[ internal.node.in.reverse.order.i, "rightChild" ] ) ) {
            # Add it to the tree as a left child.
            new.internal.node.for.tree <- the.tree[ "dummy", ];

if( debug ) {
            the.tree <- rbind( the.tree, new.internal.node.for.tree );
} else {
            the.tree <<- rbind( the.tree, new.internal.node.for.tree );
}
            new.internal.node.i <- nrow( the.tree );

if( debug ) {
            the.tree[ new.internal.node.i, "name" ] <- rownames( internal.nodes.mentioning.this.node )[ internal.node.in.reverse.order.i ];
            the.tree[ new.internal.node.i, "time" ] <- internal.nodes.mentioning.this.node[ internal.node.in.reverse.order.i, "time" ];
} else {
            the.tree[ new.internal.node.i, "name" ] <<- rownames( internal.nodes.mentioning.this.node )[ internal.node.in.reverse.order.i ];
            the.tree[ new.internal.node.i, "time" ] <<- internal.nodes.mentioning.this.node[ internal.node.in.reverse.order.i, "time" ];
}
            # Set the leftChild of the new internal node to the current node.*(come back! right now this current node is the leaf node but it won't always be)
if( debug ) {
            the.tree[ new.internal.node.i, "leftChild" ] <- 
                the.tree[ current.node.i, "name" ];
} else {
            the.tree[ new.internal.node.i, "leftChild" ] <<- 
                the.tree[ current.node.i, "name" ];

}
            # always keep updated the incoming matrix too so we know what's on the tree already.
            the.current.internal.node.row.i  <-
                which( rownames( all.the.internal.nodes ) == the.tree[ new.internal.node.i, "name" ] );

if( debug ) {
            all.the.internal.nodes[ the.current.internal.node.row.i, "leftChild" ] <-
                the.tree[ current.node.i, "name" ];
            the.tree[ current.node.i, "parent" ] <-
                the.tree[ new.internal.node.i, "name" ];
} else {
            all.the.internal.nodes[ the.current.internal.node.row.i, "leftChild" ] <<-
                the.tree[ current.node.i, "name" ];
            the.tree[ current.node.i, "parent" ] <<-
                the.tree[ new.internal.node.i, "name" ];
}

            # So now for the next iteration here, use this new internal node as the current node.
            current.node.i <- new.internal.node.i;
        } else {
            # It's already on the tree.
            internal.node.i  <- which( the.tree[ , "name" ] == rownames( internal.nodes.mentioning.this.node )[ internal.node.in.reverse.order.i ] );
            stopifnot( !is.na( the.tree[ internal.node.i, "leftChild" ] ) );
            stopifnot( is.na( the.tree[ internal.node.i, "rightChild" ] ) );

            stopifnot( the.tree[ internal.node.i, "name" ] == rownames( internal.nodes.mentioning.this.node )[ internal.node.in.reverse.order.i ] );
            stopifnot( the.tree[ internal.node.i, "time" ] == internal.nodes.mentioning.this.node[ internal.node.in.reverse.order.i, "time" ] );

            # always keep updated the incoming matrix too so we know what's on the tree already.
            the.current.internal.node.row.i  <-
                which( rownames( all.the.internal.nodes ) == the.tree[ internal.node.i, "name" ] );

            # Set the rightChild of the existing internal node to the current node.
if( debug ) {
            the.tree[ internal.node.i, "rightChild" ] <- 
                the.tree[ current.node.i, "name" ];
            # always keep updated the incoming matrix too so we know what's on the tree already.
            all.the.internal.nodes[ the.current.internal.node.row.i, "rightChild" ] <-
                the.tree[ current.node.i, "name" ];
            the.tree[ current.node.i, "parent" ] <- the.tree[ internal.node.i, "name" ];
} else {
            the.tree[ internal.node.i, "rightChild" ] <<- 
                the.tree[ current.node.i, "name" ];
            # always keep updated the incoming matrix too so we know what's on the tree already.
            all.the.internal.nodes[ the.current.internal.node.row.i, "rightChild" ] <<-
                the.tree[ current.node.i, "name" ];
            the.tree[ current.node.i, "parent" ] <<- the.tree[ internal.node.i, "name" ];
}


            # So now for the next iteration here, use this existing internal node as the current node.
            current.node.i <- internal.node.i;
        } # End if the internal node is already on the tree .. else ..
        

    } # End foreach internal node in reverse order


    return( NULL );
} # build.nodes.from.leaf (..)

result.ignored  <- sapply( nrow( all.the.leaf.nodes ):1, build.nodes.from.leaf )

# Remove dummy rows
the.tree  <- the.tree[ the.tree[ , "name" ] != "dummy", , drop = FALSE ]
rownames( the.tree )  <- the.tree[ , "name" ]


# TO DO load into ape, (compute branch lengths in ape probably automatically), write as newick

# Build an ape phylo object
phylo_object <- list()
phylo_object$edge <- edge_table
phylo_object$tip.label <- unique(sample$name)
phylo_object$node.label <- c("root", unique(sample$stock))
phylo_object$Nnode <- length(phylo_object$node.label)

## Force the class to be phylo
class(phylo_object) <- "phylo"
return(phylo_object)



################

lineList <- read.csv("line_list.csv")

# Get a list of all tip nodes and their sample times
recipient.sample.times <- lineList[ ,c("recipient", "sampleTime")]
colnames(recipient.sample.times)[colnames(recipient.sample.times) == "recipient"] = "id"
source.sample.times <- lineList[ ,c("source", "sampleTime")]
colnames(source.sample.times)[colnames(source.sample.times) == "source"] = "id"
all.tip.nodes.and.sample.times <- rbind(recipient.sample.times, source.sample.times)
all.tip.nodes.and.sample.times <- all.tip.nodes.and.sample.times[unique(all.tip.nodes.and.sample.times$id),]

# Add times to this list of tip nodes (sampling times)
all.tip.nodes.sampleTime <- c

# Get a list of all internal nodes
# These are all the transmission events, which can be created using the line list
all.internal.nodes <- as.character(paste( line$source, "->", line$recipient, sep = "" ))

all.nodes <- c(all.tip.nodes, all.internal.nodes)

# Add dates to the list of all the nodes

# Make a list of all the edges
# Edges contain a parent and a child
# Tip nodes are only children
# Internal nodes are parents with 2 children (strictly bifurcating tree)

# Start with oldest (latest infection) tip node


