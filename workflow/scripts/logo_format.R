
library(dplyr)
library(tidyr)
library(readxl)
library(dagLogo)
library(UniProt.ws)

with(snakemake@input, {
    xlsx_file <<- xlsx # "input/Sequence Logo for Andrey 20220901.xlsx"
    prot_file <<- proteome # "input/HUMAN_CRAP.fasta"
})
with(snakemake@wildcards, {
    set_name <<- set
})
with(snakemake@params, {
    residues_offset <<- residues
})
output_file <- unlist(snakemake@output)

dat <- read_xlsx(xlsx_file) %>%
    gather(set, context) %>%
    filter(set == set_name, !is.na(context), nchar(context) > 1) %>%
    pull
proteome <- prepareProteome(fasta = prot_file)
seq_data <- formatSequence(seq = dat, proteome = proteome, upstreamOffset = residues_offset, downstreamOffset = residues_offset)
saveRDS(seq_data, output_file)
