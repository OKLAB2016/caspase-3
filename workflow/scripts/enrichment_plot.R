
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpolypath)


input_file <- unlist(snakemake@input)
output_fig <- snakemake@output[["fig"]]
output_csv <- snakemake@output[["csv"]]

go_data <- readRDS(input_file)

plot_data <- lapply(go_data, `[[`, "plot_data") %>%
    bind_rows(.id = "ontology") %>%
    filter(!is.na(qvalue))

write.csv(plot_data, output_csv)

p <- ggplot(plot_data, aes(x = V1, y = V2, color = parentTerm)) +
    geom_point(aes(size = score), alpha = 0.5) +
    geom_text_repel(data = distinct(plot_data, set, parentTerm, ontology, .keep_all = T), aes(label = parentTerm), size = 2) +
    xlab("Semantic axis 1") + ylab("Semantic axis 2") +
    scale_size(trans = "log10") +
    # scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 1)) +
    scale_color_discrete(guide = "none") +
    # scale_color_brewer(palette = "Set3") +
    theme_bw() +
    facet_grid(ontology ~ set)
ggsave(output_fig, p, width = 10, height = 7)
