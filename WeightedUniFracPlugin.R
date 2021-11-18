library(microbiome)
library(ggplot2)
#library(phyloseq)
library(ape)
library(psadd)
library(vegan)
library(multtest)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")


input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }

rownames(parameters) <<- parameters[,1]; 
   # Need to get the three files
   otu.path <<- paste(pfix, toString(parameters["otufile", 2]),sep="")
   tree.path <<- paste(pfix, toString(parameters["tree", 2]), sep="")
   map.path <<- paste(pfix, toString(parameters["mapping", 2]), sep="")
}
run <- function() {
   physeq <<- read_csv2phyloseq(otu.file=otu.path, taxonomy.file=tree.path, metadata.file=map.path)
   mytree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
   physeq <<- merge_phyloseq(physeq, mytree)
}
output <- function(outputfile) {
  print("Generating plot...")
  y <- distance(physeq, method="wunifrac")
  #meta = sample_data(physeq)
  #meta = meta[,"Category"]
  #print(str(meta))
  #print(meta@.Data[1])
  #vm = as.vector(unlist(meta@.Data[1]))
  #print(c(vm))
  #print(typeof(c(vm)))
  #res = mt.minP(physeq, c(vm))
  #print(res)
  print("Generating CSV...")
  write.csv(as.matrix(y), outputfile)
  #write.csv(res, "lalala.csv")
}

