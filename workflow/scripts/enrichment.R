
library(readxl)
library(topGO)
library(dplyr)
library(tidyr)
library(rrvgo)
library(ViSEAGO)
library(tibble)
library(GOSemSim)
library(venn)
library(ggplot2)

with(snakemake@input, {
    data_file <<- data
    gene2GO_file <<- gene2GO
    all_proteins_file <<- all_proteins
    idd_proteins_file <<- idd_proteins
    uniprot_file <<- uniprot
    merops_file  <<- merops
    casbah_file  <<- casbah
    aliases_file <<- aliases
    degrab_file  <<- degrab
    topfind_file <<- topfind
})

with(snakemake@output, {
    rds_file  <<- rds
    venn_file <<- venn
})

with(snakemake@params, {
    q_val <<- q_val # 0.01
    genome <<- genome # org.Hs.eg.db
    set_names <<- sets
})

read.data <- function(fname, col, aliases, entry_names, genes, my_filter = function(.)T) {
    read.table(fname, header = T, sep = ",", fill = T, quote = '"', comment.char = "") %>%
        rename(Entry.Name = col) %>%
        separate_rows(Entry.Name, sep = "[, ]+") %>%
        separate(Entry.Name, into = c("Entry.Name", "Isoform"), sep = "-", fill = "right") %>%
        left_join(entry_names, by = "Entry.Name") %>%
        mutate(Entry2 = recode(Entry.Name, !!!aliases)) %>%
        mutate(Entry = ifelse(is.na(Entry), Entry2, Entry)) %>%
        filter(Entry %in% genes) %>%
        filter(my_filter(.)) %>%
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
degrab  <- read.data(degrab_file,  "Acc..",   aliases, entry_names, uniprot$Entry, my_filter = function(.) .$P1 %in% c("D","E") )
rawad   <- read_excel(data_file) %>%
    gather(set, gene) %>%
    filter(!is.na(gene)) %>%
    mutate(gene = recode(gene, !!!aliases)) %>%
    filter(set %in% set_names) %>%
    split(f = as.factor(.$set)) %>%
    lapply(pull, "gene")

gene2GO <- readMappings(gene2GO_file)
gene2GO_df <- read.table(gene2GO_file, col.names = c("gene", "GO"), sep = "\t") %>%
    separate_rows(GO, sep = ", ")
all_proteins <- read.table(all_proteins_file, sep = "|", quote = "", fill = T) %>%
    filter(grepl("OS=Homo sapiens", V3)) %>%
    mutate(id = recode(V2, !!!aliases)) %>%
    pull(V2)
idd_proteins <- readLines(idd_proteins_file) %>%
    recode(!!!aliases)

# all_datasets <- list(merops = merops, topfind = topfind, casbah = casbah, degrab = degrab) %>%
#     c(rawad)
all_datasets <- list(topfind = topfind) %>%
    c(rawad)
all_datasets_df <- lapply(all_datasets, data.frame) %>%
    lapply(setNames, "gene") %>%
    bind_rows(.id = "set") %>%
    left_join(gene2GO_df, by = "gene") %>%
    group_by(GO, set) %>%
    summarize(n_genes = n(), .groups = "drop") %>%
    filter(GO != "")

p <- venn(all_datasets, ggplot = T, zcolor = "style")
ggsave(venn_file, p)

gene_sets <- lapply(all_datasets, function(x) as.numeric(all_proteins %in% x)) %>%
    lapply(as.factor) %>% 
    lapply(setNames, all_proteins)

ontologies <- c("MF", "BP", "CC")
# ontologies <- c("CC")

go_data <- lapply(ontologies, function(ont) {
    topgo_objects <- lapply(gene_sets, function(x) new("topGOdata", allGenes = x, ontology = ont, annot = annFUN.gene2GO, gene2GO = gene2GO, nodeSize = 10))
    topgo_results <- lapply(topgo_objects, runTest, algorithm = "classic", statistic = "fisher")
    topgo_df      <- lapply(topgo_results, slot, "score") %>%
        lapply(as.data.frame) %>%
        lapply(rownames_to_column, "GO") %>%
        bind_rows(.id = "set") %>%
        rename(pvalue = "X[[i]]") %>%
        group_by(set) %>%
        mutate(qvalue = p.adjust(pvalue, "fdr"), score = -log10(qvalue)) %>%
        filter(qvalue <= q_val)
    best_tests <- arrange(topgo_df, qvalue) %>%
        distinct(GO, .keep_all = T) %>%
        with(setNames(score, GO))
    simMatrix <- calculateSimMatrix(names(best_tests), orgdb = genome, ont = ont, method = "Rel", semdata = godata(genome, ont = ont))
    reducedTerms <- reduceSimMatrix(simMatrix, best_tests, threshold = 0.7, orgdb = genome)
    rrvgo_coords <- scatterPlot(simMatrix, reducedTerms, size = "size") %>%
        `$`("data") %>%
        rownames_to_column("GO")

    # vesiago_objs  <- mapply(function(obj, res) merge_enrich_terms(list("obj", "res")), obj = topgo_objects, res = topgo_results)

    plot_data <- left_join(all_datasets_df, rrvgo_coords, by = "GO") %>%
        left_join(topgo_df, by = c("GO", "set")) %>%
        filter(!is.na(V1))
    list(plot_data = plot_data)
}) %>% setNames(ontologies)

saveRDS(go_data, file = rds_file)
