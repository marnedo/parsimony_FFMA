## ----html-doc, echo = FALSE--------------------------------------------------------------------------


## ----setup, include=FALSE----------------------------------------------------------------------------
source("setup.R")
library(tidyverse)


## ----eval = TRUE, echo = FALSE, message=FALSE, warning = FALSE---------------------------------------
library(ape)
library(phangorn)
library(tidyverse)


## ----eval = FALSE, echo = TRUE, results = "hide", message = FALSE, warning = FALSE-------------------
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


## ----html-doc, echo = FALSE--------------------------------------------------------------------------


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# set seed to ensure the same tree is produced
set.seed(32)
# generate a tree
tree <- rtree(n = 4, tip.label = c("a", "b", "c", "d"))


## ----eval = TRUE, echo = TRUE, message = FALSE-------------------------------------------------------
tree


## ----eval = FALSE, echo = TRUE, results = 'hide', message = FALSE------------------------------------
#> str(tree)
#> objects(tree)
#> tree$edge
#> tree$edge.length


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
plot(tree)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# set seed to ensure the same tree is produced
set.seed(32)
# generate a tree
tree <- rtree(n = 5, tip.label = c("a", "b", "c", "d", "e"))


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# get bird order data
data("bird.orders")


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# no.margin = TRUE gives prettier plots
plot(bird.orders, no.margin = TRUE)
segments(38, 1, 38, 5, lwd = 2)
text(39, 3, "Proaves", srt = 270)
segments(38, 6, 38, 23, lwd = 2)
text(39, 14.5, "Neoaves", srt = 270)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# Parrots and Passerines?
is.monophyletic(bird.orders, c("Passeriformes", "Psittaciformes"))
# hummingbirds and swifts?
is.monophyletic(bird.orders, c("Trochiliformes", "Apodiformes"))


## ----------------------------------------------------------------------------------------------------
plot(bird.orders, no.margin = TRUE)
segments(38, 1, 38, 5, lwd = 2)
text(39, 3, "Proaves", srt = 270)
segments(38, 6, 38, 23, lwd = 2)
text(39, 14.5, "Neoaves", srt = 270)
nodelabels()


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# extract clade
neoaves <- extract.clade(bird.orders, 29)
# plot
plot(neoaves, no.margin = TRUE)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# get parachtes data
parachtes <- read.phyDat("ParALL153.fas", format = "fasta")


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# perform model selection
treeRatchet  <- pratchet(parachtes, trace = 0, minit=100)
parsimony(treeRatchet, parachtes)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE, fig.show='hold'------------------
# plot MP tree
plot(treeRatchet,no.margin = TRUE)


## ----eval = FALSE, echo = TRUE, results = 'hidden', message = FALSE----------------------------------
#> # check whether the tree is rooted
#> is.rooted(treeRatchet)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# plot treeRatchet rooted
treeRatchet_r <- root(treeRatchet, "Segestria_sp_k200")
plot(treeRatchet_r, no.margin = TRUE)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# perform model selection
parsimony(treeRatchet_r, parachtes)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# perform character optimization
treeRatchet_r<-acctran(treeRatchet_r, parachtes)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# plot the rooted tree with branch lengths
plot(treeRatchet_r, type="phylogram", no.margin = TRUE)
add.scale.bar()


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# heuristic search using random addiiton of taxa
treeRA <- random.addition(parachtes)
treeSPR  <- optim.parsimony(treeRA, parachtes)


## ----eval = TRUE, echo = TRUE, results = 'hidden', message = FALSE-----------------------------------
# heuristic search using random addiiton of taxa
parsimony(c(treeRA, treeSPR), parachtes)

