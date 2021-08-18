args <- commandArgs(TRUE)
rep <- as.character(args[1])
rseed <- as.character(args[2])
source("scripts/power.code3.R")
bang_out_a_replicate(rep,rseed)

