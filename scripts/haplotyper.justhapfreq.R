args <- commandArgs(TRUE)
pool <- as.character(args[1])
OutName <- as.character(args[2])
folder <- as.character(args[3])
SNPtable <- as.character(args[4])
foundername <- as.character(args[5])

library(limSolve)
source("scripts/haplotyper.justhapfreq2.code.R")
runscan(pool, OutName, folder, SNPtable, foundername)

