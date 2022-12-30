aa = c(
    A = "#FF7979",
    C = "#FFFF00",
    D = "#C00000",
    E = "#C00000",
    F = "#F89F56",
    G = "#CC00CC",
    H = "#0070C0",
    I = "#FF7979",
    K = "#0070C0",
    L = "#FF7979",
    M = "#FF7979",
    N = "#08C81A",
    P = "#CC00CC",
    Q = "#08C81A",
    R = "#0070C0",
    S = "#08C81A",
    T = "#08C81A",
    V = "#FF7979",
    W = "#F89F56",
    Y = "#F89F56"
)
addScheme(color = aa, symbol = names(aa))

dagLogo2 <- function (testDAUresults, type = c("diff", "zscore"), pvalueCutoff = 0.05, 
    groupingSymbol = getGroupingSymbol(testDAUresults@group), 
    font = "Helvetica", fontface = "bold", fontsize = 8, title = NULL, 
    labelRelativeToAnchor = FALSE, labels = NULL, 
    alpha = 1, markers = list(), legend = FALSE, ylim = NULL) 
{
    if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") {
        stop("testDAUresults should be an object of testDAUresults\n\n         Please try ?testDAU to get help.", 
            call. = FALSE)
    }
    type <- match.arg(type)
    data <- dagLogo:::getData(type, testDAUresults)
    dat <- data$dat
    npos <- ncol(dat)
    ncha <- nrow(dat)
    if (!is.null(labels)) {
        if (length(labels) < npos) {
            stop("The length of labels is too short!", call. = FALSE)
        }
    }
    colset <- colorsets(testDAUresults@group)
    colset <- colset[rownames(dat)]
    if (any(is.na(colset))) 
        stop("Not every symbol has its color setting.", call. = FALSE)
    if (!is.null(groupingSymbol)) {
        rownames(dat) <- groupingSymbol[rownames(dat)]
        names(colset) <- groupingSymbol[names(colset)]
    }
    rname <- rownames(dat)
    if (max(nchar(rname)) > 1) {
        stop("Please using the groupingSymbol to convert", " the symbols into single letters.", 
            call. = FALSE)
    }
    key <- paste("x", ncha, font, paste(colset[rname], collapse = ""), 
        paste(rname, collapse = ""), sep = "_")
    if (exists("tmp_motifStack_symbolsCache", envir = dagLogo:::.globalEnv)) {
        symbolsCache = get("tmp_motifStack_symbolsCache", envir = dagLogo:::.globalEnv)
    }
    else {
        symbolsCache = list()
    }
    if (!is.null(symbolsCache[[key]])) {
        symbols <- symbolsCache[[key]]
    }
    else {
        symbols <- motifStack:::coloredSymbols(ncha, font, colset[rname], 
            rname, alpha = alpha, envir = dagLogo:::.globalEnv)
        symbolsCache[[key]] <- symbols
        assign("tmp_motifStack_symbolsCache", symbolsCache, envir = dagLogo:::.globalEnv)
    }
    if (legend) {
        if (is.null(groupingSymbol)) {
            groupingSymbol <- get(testDAUresults@group, envir = cachedEnv)$symbol
        }
    }
    else {
        groupingSymbol <- c(WWWW = "W")
    }
    pictureGrob <- get("pictureGrob", envir = dagLogo:::.globalEnv)
    dw <- ifelse(legend, 1/(npos + 6), 1/(npos + 2))
    grid.newpage()
    datN <- apply(dat, 2, function(.col) sum(.col[.col < 0]))
    datP <- apply(dat, 2, function(.col) sum(.col[.col > 0]))
    ylim <- c((as.integer(min(datN)/0.05) - 1) * 0.05, (as.integer(max(datP)/0.05) + 1) * 0.05)
    if (length(markers) > 0) {
        if (any(sapply(markers, function(.ele) any(nchar(.ele@label) > 
            0)))) {
            ylim[2] <- (as.integer(max(datP)/0.05) + 2) * 0.05
        }
    }
    remap <- function(x) {
        (ylim[1] - x)/(ylim[1] - ylim[2])/(1 + dw)
    }
    reheight <- function(x) {
        abs(x/(ylim[1] - ylim[2])/(1 + dw))
    }
    line1 <- as.numeric(convertUnit(unit(1, "strwidth", "W"), 
        "npc"))
    ch1 <- convertUnit(grobWidth(textGrob(label = "A", gp = gpar(fontsize = fontsize, 
        fontface = fontface))), unitTo = "npc", valueOnly = TRUE)
    cex <- dw/line1
    lwd <- cex/3
    x0 <- convertUnit(grobWidth(textGrob(label = as.character(max(ylim)), 
        gp = gpar(fontsize = fontsize, fontface = fontface))), 
        unitTo = "npc", valueOnly = TRUE)
    if (x0 > dw) {
        dw <- ifelse(legend, 1/(npos + 1 + (x0 + 2 * ch1)/dw + 
            ch1 * (max(nchar(names(groupingSymbol))) + 2)/dw), 
            1/(npos + 1 + (x0 + 2 * ch1)/dw))
    }
    else (x0 <- dw)
    x1 <- 4/5 * x0 + 2 * ch1
    x2 <- 6/5 * x0 + 2 * ch1
    x3 <- ch1 * 3/4
    x0 <- x0 + 2 * ch1
    x.pos <- 0
    grid.text(0, x0, remap(0) + dw/2, just = c(0.5, 0.5), gp = gpar(fontsize = fontsize, 
        fontface = fontface))
    grid.lines(x = x0, y = c(remap(0) + dw, 1), arrow = arrow(length = unit(0.01, 
        "npc")), gp = gpar(lwd = lwd))
    grid.lines(x = x0, y = c(remap(0), 0), arrow = arrow(length = unit(0.01, 
        "npc")), gp = gpar(lwd = lwd))
    tick <- ifelse(type == "diff", 0.1, 10^as.integer(log10(max(abs(as.numeric(dat))))))
    times <- ifelse(type == "diff", 100, 1)
    for (i in c(as.integer(min(datN)/tick):(-1), 1:as.integer(max(datP)/tick))) {
        grid.lines(x = c(x2, x0), y = remap(i * tick) + dw/2, 
            gp = gpar(lwd = lwd/2))
        grid.text(label = times * tick * i, x = x1, y = remap(i * 
            tick) + dw/2, just = c(1, 0.5), gp = gpar(fontsize = fontsize, 
            fontface = fontface))
    }
    grid.text(label = data$label, x = x3, y = remap(0) + dw/2, 
        just = "centre", rot = 90, gp = gpar(fontsize = fontsize, 
            fontface = fontface))
    x.pos <- x0 + dw/2
    y.poss <- numeric(length = npos)
    y.low.poss <- numeric(length = npos)
    for (j in 1:npos) {
        heights <- dat[, j]
        id <- order(heights)
        heights <- heights[testDAUresults@pvalue[, j] < pvalueCutoff]
        id <- id[id %in% which(testDAUresults@pvalue[, j] < pvalueCutoff)]
        id1 <- order(heights)
        y.pos <- remap(sum(heights[heights < 0]))
        y.low.poss[j] <- y.pos
        flag <- 0
        x_tick <- j
        if (labelRelativeToAnchor) {
            x_tick <- j - 1 - testDAUresults@upstreamOffset
        }
        else if (!is.null(labels)) {
            x_tick <- labels[j]
        }
        if (length(heights) > 0) {
            for (i in seq.int(length(heights))) {
                h <- reheight(heights[id1[i]])
                if (heights[id1[i]] > 0) 
                  flag <- flag + 1
                if (flag == 1) {
                  grid.text(x_tick, x.pos + dw/2, y.pos + dw/2, 
                    just = c(0.5, 0.5), gp = gpar(fontsize = fontsize * 
                      0.8, fontface = fontface))
                  y.pos <- y.pos + dw
                  flag <- flag + 1
                }
                if (h > 0) {
                  symid <- ifelse(heights[id1[i]] > 0, id[i], 
                    paste0(id[i], "_", alpha))
                  grid.draw(pictureGrob(symbols[[symid]], x.pos, 
                    y.pos, dw, h, just = c(0, 0), distort = TRUE))
                  y.pos <- y.pos + h
                }
            }
        }
        y.poss[j] <- y.pos
        if (flag == 0) {
            grid.text(x_tick, x.pos + dw/2, y.pos + dw/2, just = c(0.5, 
                0.5), gp = gpar(fontsize = fontsize * 0.8, fontface = fontface))
        }
        x.pos <- x.pos + dw
    }
    if (length(markers) > 0) {
        grid.draw(plotMarkers(markers, dw, x0 + dw/2, y.poss, 
            y.low.poss))
    }
    if (!is.null(title)) {
        grid.text(as.character(title), x = 0.5, y = 0.98, just = c(0.5, 
            0.5), gp = gpar(col = "black", fontsize = fontsize * 
            2, fontface = fontface))
    }
}
