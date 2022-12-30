
library(dplyr)
library(tidyr)
library(readxl)
library(dagLogo)
library(UniProt.ws)

snakemake@source("logo_dagLogo.R")

with(snakemake@input, {
    rds_file  <<- rds
    prot_file <<- proteome # "input/HUMAN_CRAP.fasta"
})
with(snakemake@wildcards, {
    set_name <<- set
})
output_file <- unlist(snakemake@output)

proteome <- prepareProteome(fasta = prot_file)

seq_data <- readRDS(rds_file)

bg_ztest <- buildBackgroundModel(seq_data, background = "wholeProteome", proteome = proteome, testType = "ztest")
t0 <- testDAU(seq_data, dagBackground = bg_ztest)
t0@group <- "custom"
x_labs <- c("P4", "P3", "P2", "P1", "P1'", "P2'", "P3'", "P4'")

pdf(output_file, width = 3, height = 2.5)
dagLogo2(t0, labels = x_labs, title = set_name)
dev.off()
