library(ape)
library(phangorn)
library(tidyverse)


#> # clear the R environment
#> rm(list = ls())
#> 
#> # install new packages
#> install.packages("ape")
#> install.packages("phangorn")
#> install.packages("tidyverse")
#> 
#> # load packages
#> library(ape)
#> library(phangorn)
#> library(tidyverse)

#' 
#' With these packages installed, we are ready to begin!
#' 
#' ## Phylogenetics in R

# set seed to ensure the same tree is produced
set.seed(32)
# generate a tree
tree <- rtree(n = 4, tip.label = c("a", "b", "c", "d"))
tree

#' You can actually look more deeply into the data stored within the `tree` object if you want to. Try the following code and see what is inside.
#> str(tree)
#> objects(tree)
#> tree$edge
#> tree$edge.length

plot(tree)


# set seed to ensure the same tree is produced
set.seed(32)
# generate a tree
tree <- rtree(n = 5, tip.label = c("a", "b", "c", "d", "e"))


#' ### A simple example with real data - avian phylogenetics

# get bird order data
data("bird.orders")

# no.margin = TRUE gives prettier plots
plot(bird.orders, no.margin = TRUE)
segments(38, 1, 38, 5, lwd = 2)
text(39, 3, "Proaves", srt = 270)
segments(38, 6, 38, 23, lwd = 2)
text(39, 14.5, "Neoaves", srt = 270)


# Parrots and Passerines?
is.monophyletic(bird.orders, c("Passeriformes", "Psittaciformes"))
# hummingbirds and swifts?
is.monophyletic(bird.orders, c("Trochiliformes", "Apodiformes"))


plot(bird.orders, no.margin = TRUE)
segments(38, 1, 38, 5, lwd = 2)
text(39, 3, "Proaves", srt = 270)
segments(38, 6, 38, 23, lwd = 2)
text(39, 14.5, "Neoaves", srt = 270)
nodelabels()


# extract clade
neoaves <- extract.clade(bird.orders, 29)
# plot
plot(neoaves, no.margin = TRUE)


#' ### Inferring trees with R using parsimony

# get parachtes data
parachtes <- read.phyDat("ParALL153.fas", format = "fasta")


# search for the most parsimonious (MP) tree
treeRatchet  <- pratchet(parachtes, trace = 0, minit=100, all = TRUE)
parsimony(treeRatchet, parachtes)


# plot MP tree
plot(treeRatchet,no.margin = TRUE)


#> # check whether the tree is rooted
is.rooted(treeRatchet)

# plot treeRatchet rooted
treeRatchet_r <- root(treeRatchet, "Segestria_sp_k200", resolve.root = TRUE, edgelabel = TRUE)
plot(treeRatchet_r, no.margin = TRUE)

# reporting tree length (steps)
parsimony(treeRatchet_r, parachtes)

# perform character optimization
treeRatchet_r<-acctran(treeRatchet_r, parachtes)


# plot the rooted tree with branch lengths
plot(treeRatchet_r, type="phylogram", no.margin = TRUE)
add.scale.bar()


# heuristic search using random addiiton of taxa
treeRA <- random.addition(parachtes)
treeSPR  <- optim.parsimony(treeRA, parachtes)

# compare tree length
parsimony(c(treeRA, treeSPR), parachtes)

#combines both trees into a single object
obj<-c(treeRA, treeSPR)
# unrooted pratchet tree
obj_cons <- root(consensus(obj), outgroup = "Segestria_sp_k200",resolve.root = TRUE, edgelabel =TRUE)
plot(obj_cons, main="Rooted pratchet consensus tree",no.margin = TRUE)

#' # Session info
sessionInfo()