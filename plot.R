#!/usr/bin/env Rscript
# A single renderer for audited Euclid, ranking and circle results.
binomial_interval <- function(successes, total) {
    c(if (successes == 0L) 0 else qbeta(0.025, successes, total-successes+1L),
      if (successes == total) 1 else qbeta(0.975, successes+1L, total-successes))
}

figure_data <- function(output, manifest) {
    raw <- read.csv(file.path(output, "raw.csv"), stringsAsFactors = FALSE)
    comparison <- read.csv(file.path(output, "comparisons.csv"), stringsAsFactors = FALSE)
    stopifnot(nrow(raw) == manifest$signature$expected$cutoff_rows,
              nrow(comparison) == manifest$signature$expected$comparison_rows,
              !anyDuplicated(raw[c("id", "geometry", "method")]),
              !anyDuplicated(comparison[c("id", "method")]))
    comparison$geometry <- ifelse(comparison$family == "euclidean", "euclidean", "circle")
    comparison$coverage_lower <- as.integer(comparison$decision == "IN")
    comparison$coverage_upper <- as.integer(comparison$decision != "OUT")
    comparison$attained <- NA
    fields <- c("id", "family", "parameter", "n", "replication", "geometry", "method",
                "coverage_lower", "coverage_upper", "attained")
    rows <- rbind(raw[fields], comparison[fields])
    groups <- split(rows, interaction(rows$family, rows$parameter, rows$n,
                                      rows$geometry, rows$method, drop = TRUE))
    result <- lapply(groups, function(group) {
        stopifnot(nrow(group) == manifest$signature$config$reps)
        lower <- group$coverage_lower
        upper <- group$coverage_upper
        cbind(group[1L, c("family", "parameter", "n", "geometry", "method")],
              replications = nrow(group),
              coverage_lower = mean(lower), coverage_upper = mean(upper),
              MC95_lower = binomial_interval(sum(lower), nrow(group))[1L],
              MC95_upper = binomial_interval(sum(upper), nrow(group))[2L],
              membership_unresolved = sum(lower != upper),
              cutoff_unresolved = sum(!group$attained, na.rm = TRUE))
    })
    do.call(rbind, result)
}


write_coverage <- function(data, template, path, pilot = FALSE) {
    lines <- readLines(template, warn = FALSE, encoding = "UTF-8")
    families <- c(rep("euclidean", 3), rep("ranking", 4), rep("circle", 8))
    parameters <- c(0, .5, 1, rep(NA, 4), rep(c(-.125, 0, .5, 1), 2))
    geometries <- c(rep("euclidean", 3), "euclidean", "kendall", "footrule", "graph", rep("circle", 8))
    styles <- c(G = "paperg", L = "paperl", EL_S1_RAW = "paperel",
                WALD_CHI = "paperwald", WALD_S1_KDE_CENTERED = "paperwald")
    colors <- c(G = "paperGlobal", L = "paperLocal", EL_S1_RAW = "paperEL",
                WALD_CHI = "black!75", WALD_S1_KDE_CENTERED = "black!75")
    number <- function(x) formatC(x, format = "f", digits = 6, decimal.mark = ".")
    output <- character()
    panel <- 0L
    cursor <- 1L
    while (cursor <= length(lines)) {
        line <- lines[cursor]
        output <- c(output, line)
        cursor <- cursor+1L
        if (!startsWith(line, "\\nextgroupplot")) next
        panel <- panel+1L
        family <- families[panel]
        values <- data[data$family %in% if (family == "ranking") c("ranking", "graph") else family, ]
        values <- values[if (family == "ranking") values$geometry == geometries[panel] else values$parameter == parameters[panel], ]
        methods <- switch(family, euclidean = c("G", "WALD_CHI"), ranking = c("G", "L"),
                          circle = c("G", "L", "EL_S1_RAW", "WALD_S1_KDE_CENTERED"))
        # The original paper shows identical Euclid G/L as a single Global series.
        if (family == "euclidean" && any(values$method == "L")) {
            g <- values[values$method == "G", c("n", "coverage_lower", "coverage_upper")]
            l <- values[values$method == "L", c("n", "coverage_lower", "coverage_upper")]
            if (!identical(unname(as.matrix(g)), unname(as.matrix(l)))) methods <- c("G", "L", "WALD_CHI")
        }
        if (panel <= 11L) {
            limits <- range(data$n)*c(.8, 1.25)
            output <- c(output, sprintf("\\addplot[black!55,densely dashed,line width=.45pt] coordinates {(%s,95)(%s,95)};", limits[1L], limits[2L]))
        }
        for (method in methods) {
            points <- values[values$method == method, ]
            points <- points[order(points$n), ]
            if (panel >= 8L && panel <= 11L) points <- points[100*points$MC95_upper >= 88, ]
            if (panel >= 12L) points <- points[100*points$MC95_lower <= 88, ]
            if (!nrow(points)) next
            offset <- if (family == "circle") c(G = -.09, L = -.03, EL_S1_RAW = .03, WALD_S1_KDE_CENTERED = .09)[method] else
                if (method == "G") -.0325 else .0325
            x <- points$n*2^offset
            table <- cbind(x, 100*points$coverage_lower, 100*(points$coverage_lower-points$MC95_lower),
                           100*(points$MC95_upper-points$coverage_lower))
            output <- c(output, sprintf("\\addplot[%s,error bars/y dir=both,error bars/y explicit]", styles[method]),
                        "table[x=x,y=y,y error minus=lo,y error plus=up,row sep=\\\\] {", "x y lo up\\\\",
                        apply(table, 1L, function(row) paste0(paste(number(row), collapse = " "), "\\\\")), "};")
            for (i in which(points$coverage_lower < points$coverage_upper)) {
                output <- c(output, sprintf("\\draw[%s,line width=.8pt] (axis cs:%s,%s) -- (axis cs:%s,%s);",
                    colors[method], number(x[i]), number(100*points$coverage_lower[i]), number(x[i]), number(100*points$coverage_upper[i])),
                    sprintf("\\addplot[only marks,mark=o,mark size=1.6pt,color=%s] coordinates {(%s,%s)};",
                            colors[method], number(x[i]), number(100*points$coverage_upper[i])))
            }
        }
        while (cursor <= length(lines) && !startsWith(lines[cursor], "\\nextgroupplot") &&
               !startsWith(lines[cursor], "\\end{groupplot}")) cursor <- cursor+1L
    }
    stopifnot(panel == 15L)
    count <- unique(data$replications)
    stopifnot(length(count) == 1L)
    output <- sub("1,000 datasets", paste0(format(count, big.mark = ",", trim = TRUE), " datasets"), output, fixed = TRUE)
    if (any(data$coverage_lower < data$coverage_upper)) {
        output <- sub("individual 95\\% exact binomial", "95\\% exact binomial", output, fixed = TRUE)
        output <- sub("Monte Carlo intervals", "Monte Carlo envelopes", output, fixed = TRUE)
        output <- sub("above; no percentage range is omitted.}",
                      "above; no percentage range is omitted. Numerical membership bounds are shown.}", output, fixed = TRUE)
    }
    if (pilot) output <- sub("\\caption{", "\\caption{PILOT ONLY. ", output, fixed = TRUE)
    writeLines(output, path, useBytes = TRUE)
}

main <- function(arguments = commandArgs(TRUE)) {
    pilot <- "--pilot" %in% arguments
    arguments <- arguments[arguments != "--pilot"]
    if (any(startsWith(arguments,"--"))) stop("Unsupported plotting option.")
    if (!length(arguments) || length(arguments)>2L) stop("Usage: Rscript plot.R RUN [NEW_DIRECTORY] [--pilot]")
    script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
    script <- normalizePath(script, mustWork = TRUE)
    template <- file.path(dirname(script), "coverage-layout.tex")
    source(file.path(dirname(script),"common","reproduce.R"),local=.GlobalEnv)
    output <- normalizePath(arguments[1L], mustWork = TRUE)
    manifest <- readRDS(file.path(output, "manifest.rds"))
    status <- readRDS(file.path(output, "STATUS.rds"))
    hashes <- tools::md5sum(file.path(output, "source", names(manifest$signature$source_md5)))
    stopifnot(identical(unname(hashes), unname(manifest$signature$source_md5)),
              grepl("^COMPLETE_", status$status), identical(status$signature, manifest$signature_hash))
    hashes <- tools::md5sum(file.path(output, names(status$report_md5)))
    stopifnot(length(hashes) == 4L, identical(unname(hashes), unname(status$report_md5)))
    if (!pilot) stopifnot(manifest$signature$config$reps == 1000L)
    destination <- if (length(arguments) >= 2L) arguments[2L] else file.path(output, "figures")
    data <- figure_data(output, manifest)
    provenance <- list(run_signature = manifest$signature_hash, report_md5 = status$report_md5, pilot = pilot)
    if (file.exists(destination)) stop("Choose a fresh figure directory.")
    dir.create(destination, recursive = TRUE)
    destination <- normalizePath(destination, mustWork = TRUE)
    fragment <- file.path(destination, "figure-coverage.tex")
    write_coverage(data, template, fragment, pilot)
    fields <- c("family","parameter","n","geometry","method","replications",
                "coverage_lower","coverage_upper","MC95_lower","MC95_upper","membership_unresolved")
    write.csv(data[fields], file.path(destination, "figure-data.csv"), row.names = FALSE)
    provenance$renderer_md5 <- tools::md5sum(c(script, template))
    provenance$files_md5 <- tools::md5sum(file.path(destination,
        c("figure-coverage.tex", "figure-data.csv")))
    publication <- tempfile(".publish-figure-manifest.rds-", tmpdir = destination)
    on.exit(unlink(publication), add = TRUE)
    saveRDS(provenance, publication, version = 3L)
    stopifnot(file.rename(publication, file.path(destination, "figure-manifest.rds")))
    cat("TeX figure written to", destination, "\n")
}

if (sys.nframe() == 0L) main()
