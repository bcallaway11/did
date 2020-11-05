## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
fldr <- "~/Dropbox/did/R/"
sapply(paste0(fldr,list.files(fldr)), source)
library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)

## ----echo=FALSE---------------------------------------------------------------
source("~/Dropbox/did/tests/setup_sims.R")
time.periods <- 5
reset.sim()
te <- 0
te.e <- c(-1,rep(1,time.periods-1))
bett <- betu <- rep(0,time.periods)
#set.seed(1814)

# generate data set with these parameters
data <- build_sim_dataset()
data$G <- ifelse(data$G==0, 0, data$G+1)
data$G <- ifelse(data$G==6, 0, data$G)

## -----------------------------------------------------------------------------
nrow(data)
head(data)

