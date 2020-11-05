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

## ----echo=FALSE---------------------------------------------------------------
source("~/Dropbox/did/tests/setup_sims.R")
time.periods <- 4
reset.sim()
te <- 0
set.seed(1814)

## -----------------------------------------------------------------------------
# generate dataset with 4 time periods
time.periods <- 4

# generate data set with these parameters
data <- build_sim_dataset()

nrow(data)
head(data)

