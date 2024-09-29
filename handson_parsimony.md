---
title: "Intro and parsimony exercises FFMA"
author: "Miquel Arnedo"
date: "2024-09-29"
output:
  html_document:
    keep_md: true
---



In this tutorial, we will introduce phylogenetics  as a means to visualise the evolutionary relationships among species and will coduct our first phylogenetic analyses using parismony as method of inference

### What to expect {.unnumbered}

In this section we will:

-   learn some tools for visualising phylogenetic trees
-   learn how to create phylogenies

### Getting started {.unnumbered}

First, we need to set up our R environment. We'll load `tidyverse` a package that facilitaes dta amaipulation and visualization. aalong a few more packages today to help us handle different types of data. Chief among these is `ape` which is the basis for a lot of phylogenetic analysis in R. We will also load another phylogenetic package, `phangorn` (which has an extremely [geeky reference](https://en.wikipedia.org/wiki/Fangorn) in its name).



``` r
# clear the R environment
rm(list = ls())

# install new packages
install.packages("ape")
install.packages("phangorn")
install.packages("tidyverse")

# load packages
library(ape)
library(phangorn)
library(tidyverse)
```


