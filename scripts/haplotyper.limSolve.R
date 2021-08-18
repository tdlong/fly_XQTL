args <- commandArgs(TRUE)
Cpool <- as.character(args[1])
Tpool <- as.character(args[2])
OutName <- as.character(args[3])
folder <- as.character(args[4])
SNPtable <- as.character(args[5])
foundername <- as.character(args[6])

library(limSolve)
source("scripts/haplotyper.limSolve.code.R")
runscan(Cpool, Tpool, OutName, folder, SNPtable, foundername)

