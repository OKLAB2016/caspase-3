
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

data_file <- unlist(snakemake@input)
with(snakemake@output, {
    jit_file <<- jit
    dot_file <<- dot
})

trend_factor <- function(Trends, pos) {
    factor(substr(Trends, pos, pos), levels = c("L", "M", "H"))
}
data <- read_excel(data_file, .name_repair = "universal") %>%
    mutate(In.Vitro = as.numeric(ifelse(In.Vitro == "NA", NA, In.Vitro))) %>%
    mutate(Trends = ifelse(Trends == "NA", NA, Trends)) %>%
    mutate(In.Vitro.Trend = trend_factor(Trends, 1), In.Vivo.Trend = trend_factor(Trends, 2))

p <- ggplot(data, aes(In.Vitro, In.Vivo, color = In.Vitro > 1 & In.Vivo > 1 | In.Vitro < -1 & In.Vivo < -1)) +
    geom_hline(yintercept = 1, linetype = "dotted") + geom_hline(yintercept = -1, linetype = "dotted") +
    geom_vline(xintercept = 1, linetype = "dotted") + geom_vline(xintercept = -1, linetype = "dotted") +
    geom_point() +
    scale_color_manual(values = c("black", "red")) +
    xlab("In Vitro") + ylab("In Vivo") +
    theme_bw() +
    theme(legend.position = "none")
ggsave(dot_file, p, w = 4, h = 4)

p <- ggplot(data, aes(In.Vitro.Trend, In.Vivo.Trend)) +
    geom_jitter() +
    xlab("In Vitro") +
    ylab("In Vivo") +
    theme_bw()
ggsave(jit_file, p, w = 4, h = 4)

