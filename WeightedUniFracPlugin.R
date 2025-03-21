library(microbiome)
library(ggplot2)
#library(phyloseq)
library(ape)
library(psadd)
library(picante)
library(httr)
library(foreach)
verbose()
library(sna)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


fastUniFrac <- function(physeq, weighted=FALSE, normalized=TRUE, parallel=FALSE){
	# Access the needed components. Note, will error if missing in physeq.
	OTU  <- otu_table(physeq)
	tree <- phy_tree(physeq)
	# Some important checks.
	if( is.null(tree$edge.length) ) {
	  stop("Tree has no branch lengths, cannot compute UniFrac")
	}
	if( !is.rooted(tree) ) {
	  stop("Rooted phylogeny required for UniFrac calculation")
	}
	### Some parallel-foreach housekeeping.
	# If user specifies not-parallel run (the default), register the sequential "back-end"
	#if( !parallel ){ registerDoSEQ() }
	# create N x 2 matrix of all pairwise combinations of samples.
	spn <- combn(sample_names(physeq), 2, simplify=FALSE)
	# Make sure OTU is in species-are-rows orientation
	if( !taxa_are_rows(physeq) ){OTU <- t(OTU)}
  # Convert to standard matrix
	OTU <- as(OTU, "matrix")
	# Enforce that tree and otu_table indices are the same order,
	# by re-ordering OTU, if needed
	if( !all(rownames(OTU) == taxa_names(tree)) ){
	  OTU <- OTU[taxa_names(tree), ]
	}
	########################################
	# Build the requisite matrices as defined
	# in the Fast UniFrac article.
	########################################
	## This only needs to happen once in a call to UniFrac.
	## Notice that A and B do not appear in this section.
	# Begin by building the edge descendants matrix (edge-by-sample)
  # `edge_array`
  #
	# Create a list of descendants, starting from the first internal node (root)
	ntip <- length(tree$tip.label)
	#print(str(tree))
	#print(tree$node.label)
	#print(tree$tip.label)
	#tree$edge = rbind(tree$edge, c(38, 40))
	#print(tree$edge)
	#print("**********")
	#print("EDGES")
	#print(tree$edge)
	#print(order(tree$edge[,1]))
        #print(tree$edge[order(tree$edge[,1]),][,2])
	#print("**********")
	if(ntip != ntaxa(physeq)) stop("Incompatible tree and OTU table!")
	# Create a matrix that maps each internal node to its 2 descendants
	# This matrix doesn't include the tips, so must use node#-ntip to index into it
	node.desc <- matrix(tree$edge[order(tree$edge[,1]),][,2],byrow=TRUE,ncol=2)
	# Define the edge_array object
	# Right now this is a node_array object, each row is a node (including tips)
	# It will be subset and ordered to match tree$edge later
	edge_array <- matrix(0, nrow=ntip+tree$Nnode, ncol=nsamples(physeq),
	                     dimnames=list(NULL, sample_names=sample_names(physeq)))
	# Load the tip counts in directly
	edge_array[1:ntip,] <- OTU
	# Get a list of internal nodes ordered by increasing depth
	ord.node <- order(node.depth(tree))[(ntip+1):(ntip+tree$Nnode)]
	# Loop over internal nodes, summing their descendants to get that nodes count
	#print("INTERNAL NODES")
	#print(node.desc)
	#print("INTERNAL NODES ORDERED BY INCREASING DEPTH")
	#print(ord.node)
        #print("NTIP")
	#print(ntip)
	#print("NTAXA")
	#print(ntaxa(physeq))
	for(i in ord.node){
		#print(i)
      	        #print(i-ntip)
	  edge_array[i,] <- colSums(edge_array[node.desc[i-ntip,], , drop=FALSE], na.rm = TRUE)
	}
	# Keep only those with a parental edge (drops root) and order to match tree$edge
	edge_array <- edge_array[tree$edge[,2],]
	# Remove unneeded variables.
	rm(node.desc)
	# If unweighted-UniFrac, coerce to a presence-absence contingency, occ
	if(!weighted){
		# For unweighted UniFrac, convert the edge_array to an occurrence (presence/absence binary) array
		edge_occ <- (edge_array > 0) - 0
	}
	if( weighted & normalized ){
		# This is only relevant to weighted-UniFrac.
		# For denominator in the normalized distance, we need the age of each tip.
	  # 'z' is the tree in postorder order used in calls to .C
	  # Descending order of left-hand side of edge (the ancestor to the node)
	  z = reorder.phylo(tree, order="postorder")
	  # Call phyloseq-internal function that in-turn calls ape's internal
	  # horizontal position function, in C, using the re-ordered phylo object, `z`
	  tipAges = node.depth.edgelength(tree)
	  # Keep only the tips, and add the tip labels in case `z` order differs from `tree`
	  tipAges <- tipAges[1:length(tree$tip.label)]
	  names(tipAges) <- z$tip.label
    # Explicitly re-order tipAges to match OTU
	  tipAges <- tipAges[rownames(OTU)]
	}
	########################################
  # optionally-parallel implementation with foreach
	########################################
	samplesums = sample_sums(physeq)
	distlist <- foreach( i = spn, .packages="phyloseq") %do%  {
	  A  <- i[1]
	  B  <- i[2]
	  AT <- samplesums[A]
	  BT <- samplesums[B]
	  if( weighted ){
      # weighted UniFrac
	    wUF_branchweight <- abs(edge_array[, A]/AT - edge_array[, B]/BT)
	    # calculate the w-UF numerator
	    numerator <- sum({tree$edge.length * wUF_branchweight}, na.rm = TRUE)
	    # if not-normalized weighted UniFrac, just return "numerator";
	    # the u-value in the w-UniFrac description
	    if(!normalized){
	      return(numerator)
	    } else {
	      # denominator (assumes tree-indices and otu_table indices are same order)
	      denominator <- sum({tipAges * (OTU[, A]/AT + OTU[, B]/BT)}, na.rm = TRUE)
	      # return the normalized weighted UniFrac values
	      return(numerator / denominator)
	    }
	  } else {
      # Unweighted UniFrac
	    # Subset matrix to just columns A and B
	    edge_occ_AB <- edge_occ[, c(A, B)]
      # Keep only the unique branches. Sum the lengths
      edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[rowSums(edge_occ_AB, na.rm=TRUE) < 2, ], na.rm=TRUE)
	    # Normalize this sum to the total branches among these two samples, A and B
	    uwUFpairdist <- edge_uni_AB_sum / sum(tree$edge.length[rowSums(edge_occ_AB, na.rm=TRUE) > 0])
	    return(uwUFpairdist)
	  }
	}
	# Initialize UniFracMat with NAs
	UniFracMat <- matrix(NA_real_, nsamples(physeq), nsamples(physeq))
	rownames(UniFracMat) <- colnames(UniFracMat) <- sample_names(physeq)
  # Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
  	matIndices <- do.call(rbind, spn)[, 2:1]
  	# Take care of edge case where there are two samples -> 1 pair of indices -> rbind doesn't return a matrix
  	if(!is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol=2)
	UniFracMat[matIndices] <- unlist(distlist)
	return(as.dist(UniFracMat))
}

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1]; 
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
   otu.path <<- paste(pfix, toString(parameters["otufile", 2]),sep="")
   tax.path <<- paste(pfix, toString(parameters["tax", 2]), sep="")
   map.path <<- paste(pfix, toString(parameters["mapping", 2]), sep="")
   tree.path <<- paste(pfix, toString(parameters["tree", 2]), sep="")
   userandom <<- toString(parameters["userandom", 2])
}

run <- function() {
	#print(pfix)
   physeq <<- read_csv2phyloseq(otu.file=otu.path, taxonomy.file=tax.path, metadata.file=map.path)
   if (userandom == "True") {
    print("RANDOM TREE")
	   mytree <<- rcoal(ntaxa(physeq), rooted=TRUE, br=1, tip.label=taxa_names(physeq))
           #print("BBB")
	   #OTU <- otu_table(physeq)
	   #rn <- rownames(OTU)
	   #cn <- colnames(OTU)
	   #for (x in 1:nrow(OTU)*ncol(OTU)) {
           #   x1 <- sample.int(nrow(OTU)*ncol(OTU), 1)-1#sample(1:nrow(OTU)*ncol(OTU), 1)
	   #   i1 <- as.integer(x1 / ncol(OTU))+1
	   #   j1 <- x1 %% ncol(OTU)+1
	   #   x2 <- sample.int(nrow(OTU)*ncol(OTU), 1)-1#sample(1:nrow(OTU)*ncol(OTU), 1)
	   #   i2 <- as.integer(x2 / ncol(OTU))+1
	   #   j2 <- x2 %% ncol(OTU)+1
	      #print("*******")
	      #print(nrow(OTU))
	      #print(ncol(OTU))
	      #print("X")
	      #print(x1)
	      #print("I")
	      #print(i1)
	      #print("J")
	      #print(j1)
	   #   print(x2)
	   #   print(i2)
	   #   print(j2)
	   #   print("*******")
	   #   tmp <- OTU[[i1, j1]]
	   #   print(tmp)
	   #   print("A")
	   #   OTU[[i1, j1]] <- OTU[[i2, j2]]
	   #   OTU[[i2, j2]] <- tmp
           #}
	   #print("CCC")
           #rp <- sample(OTU, nrow(OTU)*ncol(OTU))
	   #rp <- OTU[sample(nrow(OTU), nrow(OTU), replace=TRUE),]
	   #rp <- rp[,sample(ncol(OTU), ncol(OTU), replace=TRUE)]
	   #rp <- rmperm(OTU)
	   #print(rp)
	   #print(rn)
	   #print(cn)
	   #print(nrow(rp))
	   #print(ncol(rp))
	   #row.names(rp) <- rn
	   #colnames(rp) <- cn
	   #OTU <- otu_table(rp, taxa_are_rows=TRUE)
           #OTU <<- otu_table(rperm(otu_table(physeq)))
	   #print(OTU)
           #physeq <<- phyloseq(OTU, tax_table(physeq), sample_data(physeq))	   
   }
   else {
	   #print("NOT RANDOM")
	   #print(userandom)
	   #print(tree.path)
   mytree <<- read_tree(tree.path)
   }
   #print(mytree)
   physeq <<- merge_phyloseq(physeq, mytree)
   ape::write.tree(phy_tree(physeq), "output.tree")
}
output <- function(outputfile) {
  #print("Generating plot...")
  #print(rownames(otu_table(physeq)))
  #print(mytree)
  #print(rbiom::tips(mytree))
  #y <- rbiom::unifrac(otu_table(physeq), weighted=TRUE, tree=mytree)
  #y <- picante::unifrac(otu_table(physeq), mytree)
  y <- fastUniFrac(physeq, weighted=TRUE)
  #y <- distance(physeq, method="wunifrac")
  #print("Generating CSV...")
  write.csv(as.matrix(y), outputfile)
}

