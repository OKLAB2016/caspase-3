library(dplyr)
library(ontologyIndex)

input <- unlist(snakemake@input)
output <- unlist(snakemake@output)
get_ontology(input) %>%
    with(name[grepl("GO", names(name))]) %>%
    data.frame %>%
    write.table(output, sep = "\t", quote = F, col.names = F, row.names = T)
