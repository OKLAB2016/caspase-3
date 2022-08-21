
library(readxl)
library(topGO)
library(dplyr)
library(tidyr)
library(rrvgo)
library(ViSEAGO)
library(tibble)
library(ggplot2)
library(ggrepel)
library(GOSemSim)
library(venn)
library(ggpolypath)

with(snakemake@input, {
    data_file <<- data
    gene2GO_file <<- gene2GO
    all_proteins_file <<- all_proteins
    uniprot_file <<- uniprot
    merops_file  <<- merops
    casbah_file  <<- casbah
    aliases_file <<- aliases
    degrab_file  <<- degrab
    topfind_file <<- topfind
})

with(snakemake@output, {
    rrvgo_file <<- rrvgo
    venn_file  <<- venn
})

with(snakemake@params, {
    q_val <<- q_val # 0.01
    genome <<- genome # org.Hs.eg.db
})

read.data <- function(fname, col, aliases, entry_names, genes) {
    read.table(fname, header = T, sep = ",", fill = T, quote = '"', comment.char = "") %>%
        rename(Entry.Name = col) %>%
        separate_rows(Entry.Name, sep = "[, ]+") %>%
        separate(Entry.Name, into = c("Entry.Name", "Isoform"), sep = "-", fill = "right") %>%
        left_join(entry_names, by = "Entry.Name") %>%
        mutate(Entry2 = recode(Entry.Name, !!!aliases)) %>%
        mutate(Entry = ifelse(is.na(Entry), Entry2, Entry)) %>%
        filter(Entry %in% genes) %>%
        distinct(Entry) %>% pull
}

uniprot <- read.table(uniprot_file, sep = "\t", header = T, quote = "", fill = T)
entry_names <- select(uniprot, Entry, Entry.Name)
aliases <- read.table(aliases_file, col.names = c("Entry.Name", "Entry")) %>%
    with(setNames(Entry, Entry.Name))

gene2GO <- readMappings(gene2GO_file)
genes <- uniprot$Entry

merops  <- read.data(merops_file,  "Uniprot", aliases, entry_names, uniprot$Entry)
topfind <- read.data(topfind_file, "ac",      aliases, entry_names, uniprot$Entry)
casbah  <- read.data(casbah_file,  "UniProt", aliases, entry_names, uniprot$Entry)
degrab  <- read.data(degrab_file,  "Acc..",   aliases, entry_names, uniprot$Entry)
rawad   <- read_excel(data_file) %>%
    gather(set, gene) %>%
    filter(!is.na(gene)) %>%
    filter(set == "All D/E  n>1 Ratio<-1") %>%
    pull(gene)
all_datasets <- list(merops = merops, topfind = topfind, casbah = casbah, degrab = degrab, rawad = rawad)

p <- venn(all_datasets, ggplot = T, zcolor = "style")
ggsave(venn_file, p)

ontologies <- c("MF", "BP", "CC")
semdata <- lapply(ontologies, function(ont) godata(genome, ont = ont)) %>%
    setNames(ontologies)

gene2GO <- readMappings(gene2GO_file)
gene2GO_df <- read.table(gene2GO_file, col.names = c("gene", "GO"), sep = "\t") %>%
    separate_rows(GO, sep = ", ")
all_proteins <- pull(read.table(all_proteins_file))

topgo_test <- function(genes) {
    lapply(ontologies, function(ontology) {
        de_genes <- all_proteins %in% genes %>%
            as.numeric %>%
            as.factor %>%
            setNames(all_proteins)

        test <- new("topGOdata", ontology = ontology, allGenes = de_genes, annot = annFUN.gene2GO, gene2GO = gene2GO, nodeSize = 10) %>%
            runTest(algorithm = "classic", statistic = "fisher") %>%
            `@`("score") %>%
            data.frame(pvalue = ., ID = names(.)) %>%
            mutate(qvalue = p.adjust(pvalue, "fdr"), score = -log10(qvalue)) %>%
            filter(qvalue <= q_val)
    }) %>% setNames(ontologies) %>%
        bind_rows(.id = "ontology")
}
rrvgo_coords <- function(best_tests) {
    lapply(ontologies, function(ont) {
        test <- filter(best_tests, ontology == ont)
        simMatrix <- filter(test, ontology == ont) %>%
            with(calculateSimMatrix(ID, orgdb = genome, ont = ont, method = "Rel", semdata = semdata[[ont]]))
        scores <- with(test, setNames(score, ID))
        reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold = 0.7, orgdb = genome)
        scatterPlot(simMatrix, reducedTerms, size = "size")$data
    }) %>% setNames(ontologies) %>%
        bind_rows(.id = "ontology")
}

tests <- lapply(all_datasets, topgo_test)
best_tests <- bind_rows(tests, .id = "set") %>%
    arrange(qvalue) %>%
    distinct(ID, .keep_all = T)
coords <- rrvgo_coords(best_tests) %>%
    rownames_to_column("GO")
plot_data <- lapply(all_datasets, data.frame) %>%
    lapply(setNames, "gene") %>%
    bind_rows(.id = "set") %>%
    left_join(gene2GO_df, by = "gene") %>%
    group_by(GO, set) %>%
    summarize(n_genes = n(), .groups = "drop") %>%
    left_join(coords, by = "GO") %>%
    filter(!is.na(V1))
p <- ggplot(plot_data, aes(x = V1, y = V2, color = parentTerm)) +
    geom_point(aes(size = n_genes), alpha = 0.5) +
    geom_text_repel(data = distinct(enrichment, set, ontology, parentTerm, .keep_all = T), aes(label = parentTerm)) +
    xlab("Semantic axis 1") + ylab("Semantic axis 2") +
    scale_size(trans = "log10") +
    scale_color_discrete(guide = "none") +
    theme_bw() +
    facet_grid(ontology ~ set)
ggsave(rrvgo_file, p, width = 20, height = 12)
