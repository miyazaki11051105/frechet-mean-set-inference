# Reproduction checks use base R; numerical fitting remains in inference.R.
sha256_file <- function(path) {
    result <- system2("sha256sum", shQuote(normalizePath(path, mustWork = TRUE)), stdout = TRUE)
    if (length(result) != 1L || !is.null(attr(result, "status"))) stop("SHA-256 check failed: ", path)
    sub(" .*", "", result)
}

check_report_hashes <- function(output, hashes) {
    # Resolve reports within the supplied run directory.
    paths <- file.path(output, basename(names(hashes)))
    stopifnot(length(paths) > 0L, !anyDuplicated(paths), all(file.exists(paths)),
              identical(unname(tools::md5sum(paths)), unname(hashes)))
    invisible(paths)
}

audit_csv <- function(path, character, numeric, logical = character(), nullable = character(),
                      extended = character()) {
    schema <- c(setNames(rep("character", length(character)), character),
        setNames(rep("numeric", length(numeric)), numeric),
        setNames(rep("logical", length(logical)), logical))
    # Character import accepts quoted numeric exports without allowing malformed
    # numeric/logical fields to disappear through type inference or NA coercion.
    rows <- read.csv(path, colClasses = "character", na.strings = NULL,
        stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(nrow(rows) > 0L, !anyDuplicated(names(rows)), setequal(names(rows), names(schema)))
    for (field in names(schema)) {
        value <- rows[[field]]
        missing <- value == "NA" & field %in% nullable
        if (schema[[field]] == "numeric") value <- suppressWarnings(as.numeric(value))
        if (schema[[field]] == "logical") {
            stopifnot(all(value %in% c("TRUE", "FALSE") | missing))
            value <- value == "TRUE"
        }
        value[missing] <- NA
        stopifnot(all(is.na(value) == missing))
        if (!field %in% nullable) stopifnot(!anyNA(value))
        if (schema[[field]] == "character") stopifnot(all(nzchar(value[!is.na(value)])))
        if (schema[[field]] == "numeric") stopifnot(!any(is.nan(value)),
            field %in% extended || all(is.finite(value[!is.na(value)])))
        rows[[field]] <- value
    }
    rows
}

check_report_inventory <- function(rows, expected, keys) {
    stopifnot(all(names(expected) %in% names(rows)))
    key <- function(x) do.call(paste, c(x[keys], sep = "|"))
    observed <- key(rows); required <- key(expected)
    stopifnot(nrow(rows) == nrow(expected), !anyDuplicated(observed), !anyDuplicated(required),
        setequal(observed, required))
    index <- match(observed, required)
    for (field in names(expected)) stopifnot(all(rows[[field]] == expected[[field]][index]))
    invisible(TRUE)
}

check_summary <- function(rows, summary, reps) {
    key <- function(x) paste(x$family, x$parameter, x$n, x$geometry, x$method, sep = "|")
    stopifnot("coverage_unresolved" %in% names(summary),
              !anyDuplicated(paste(rows$id, rows$geometry, rows$method)),
              !anyDuplicated(key(summary)),
              all(rows$coverage_lower %in% 0:1), all(rows$coverage_upper %in% 0:1),
              all(rows$coverage_lower <= rows$coverage_upper))
    groups <- split(rows, key(rows))
    stopifnot(setequal(names(groups), key(summary)))
    for (i in seq_len(nrow(summary))) {
        z <- groups[[key(summary)[i]]]; total <- nrow(z)
        lo <- sum(z$coverage_lower); hi <- sum(z$coverage_upper)
        expected <- c(lo/total, hi/total)
        fields <- c("coverage_lower", "coverage_upper")
        stopifnot(all(fields %in% names(summary)))
        observed <- as.numeric(unlist(summary[i, fields], use.names=FALSE))
        stopifnot(total == reps, summary$replications[i] == reps,
                  summary$attained[i] == sum(z$attained), summary$unresolved[i] == sum(!z$attained),
                  summary$coverage_resolved[i] == sum(z$coverage_resolved),
                  summary$coverage_unresolved[i] == sum(z$coverage_lower != z$coverage_upper),
                  summary$maximum_cutoff_gap[i] == max(z$gap),
                  all(z$evaluation_scope == summary$evaluation_scope[i]),
                  all(abs(expected-observed) <= 5e-14))
    }
    invisible(TRUE)
}

require_current_run <- function(output) {
    if (file.exists(file.path(output, "repair", "manifest.rds")))
        stop("Use the historical run's preserved source, or choose a fresh run directory for this version.")
    invisible(TRUE)
}

audit_coverage_run <- function(output) {
    output <- normalizePath(output, mustWork = TRUE)
    require_current_run(output)
    manifest <- readRDS(file.path(output, "manifest.rds"))
    status <- readRDS(file.path(output, "STATUS.rds"))
    cfg <- manifest$signature$config
    stopifnot(all(c("status","signature","expected","report_md5","attained","unresolved",
                    "coverage_resolved","coverage_unresolved") %in% names(status)),
              is.character(status$status), length(status$status) == 1L, !is.na(status$status),
              all(vapply(status[c("attained","unresolved","coverage_resolved","coverage_unresolved")],
                  function(x) is.numeric(x) && length(x)==1L && is.finite(x) && x>=0 && x==floor(x),logical(1))),
              grepl("^COMPLETE_", status$status),
              manifest$signature_hash == object_hash(manifest$signature),
              status$signature == manifest$signature_hash,
              identical(status$expected, manifest$signature$expected),
              identical(status$expected, expected_counts(cfg)),
              length(status$report_md5) == 4L,
              setequal(names(status$report_md5),
                  c("raw.csv", "summary.csv", "comparisons.csv", "comparison_summary.csv")))
    files <- names(manifest$signature$source_md5)
    stopifnot(identical(unname(tools::md5sum(file.path(output, "source", files))),
                        unname(manifest$signature$source_md5)))
    check_report_hashes(output, status$report_md5)
    plan <- manifest$signature$plan
    check_report_inventory(plan, simulation_plan(cfg)[c("id", "family", "parameter", "n", "replication")], "id")
    stopifnot(!anyDuplicated(plan$id), nrow(plan) == status$expected$tasks,
        is.character(plan$input_md5), !anyNA(plan$input_md5), all(grepl("^[0-9a-f]{32}$", plan$input_md5)),
        identical(unname(tools::md5sum(file.path(output,"inputs",paste0(plan$id,".rds")))),
                  unname(plan$input_md5)))
    raw <- audit_csv(file.path(output,"raw.csv"),
        c("id","family","geometry","method","status","evaluation_scope","input_md5","draw_hash","signature","objective"),
        c("parameter","n","replication","B","k","cutoff_lower","cutoff_upper","cutoff_used","gap",
          "coverage_lower","coverage_upper","stage","max_draw_gap","elapsed_seconds"),
        c("attained","coverage_resolved"), extended=c("cutoff_lower","cutoff_upper","cutoff_used","gap","max_draw_gap"))
    summary <- audit_csv(file.path(output,"summary.csv"),
        c("family","geometry","method","evaluation_scope"),
        c("parameter","n","replications","attained","unresolved","coverage_resolved","coverage_unresolved",
          "coverage_lower","coverage_upper","maximum_cutoff_gap"), extended="maximum_cutoff_gap")
    comparisons <- audit_csv(file.path(output,"comparisons.csv"),
        c("id","family","method","decision","assurance","evaluation_scope"),
        c("parameter","n","replication","cutoff","statistic_lower","statistic_upper","center","radius","curvature",
          "score_variance","numerical_coverage"), "numerically_certified",
        nullable=c("statistic_lower","statistic_upper","center","radius","curvature","score_variance","numerical_coverage"),
        extended=c("statistic_lower","statistic_upper"))
    expected_raw <- do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
        geometry <- if (plan$family[i] == "ranking") c("euclidean","kendall","footrule") else plan$family[i]
        methods <- expand.grid(geometry=geometry,method=c("G","L"),stringsAsFactors=FALSE)
        cbind(plan[rep(i,nrow(methods)),c("id","family","parameter","n","replication")],methods)
    }))
    expected_raw$evaluation_scope <- ifelse(expected_raw$family == "graph",
        "EXACT_RATIONAL_WHOLE_POPULATION_MEAN_SET",ifelse(expected_raw$family == "circle" & expected_raw$parameter < 0,
            "WHOLE_POPULATION_MEAN_SET_WITH_CERTIFIED_QUERY_ENCLOSURES","WHOLE_POPULATION_MEAN_SET"))
    expected_comparisons <- do.call(rbind, lapply(which(plan$family %in% c("euclidean","circle")), function(i) {
        methods <- if (plan$family[i] == "euclidean") c("WALD_CHI","WALD_F_GAUSSIAN_REFERENCE") else
            c("EL_S1_RAW","WALD_S1_KDE_CENTERED")
        cbind(plan[rep(i,2L),c("id","family","parameter","n","replication")],method=methods)
    }))
    expected_comparisons$evaluation_scope <- ifelse(expected_comparisons$family == "circle" & expected_comparisons$parameter < 0,
        ifelse(expected_comparisons$method == "EL_S1_RAW","WHOLE_POPULATION_MEAN_SET_WITH_CERTIFIED_QUERY_ENCLOSURES",
            "BINARY64_REPRESENTATIVES_OF_NONBINARY_POPULATION_MEANS"),"WHOLE_POPULATION_MEAN_SET")
    expected_comparisons$assurance <- ifelse(expected_comparisons$family == "euclidean",
        "EXACT_RATIONAL_INVERSION_FIXED_BINARY64_THRESHOLD",ifelse(expected_comparisons$method == "EL_S1_RAW",
            "TRUE_PI_INTERVAL_EL_STATISTIC_AT_ENCLOSED_QUERY_AND_FIXED_BINARY64_CUTOFF",
            "FLOAT64_KDE_WALD_REGULARITY_NOT_CERTIFIED"))
    check_report_inventory(raw, expected_raw, c("id","geometry","method"))
    check_report_inventory(summary, unique(expected_raw[c("family","parameter","n","geometry","method","evaluation_scope")]),
        c("family","parameter","n","geometry","method"))
    check_report_inventory(comparisons, expected_comparisons, c("id","method"))
    stopifnot(nrow(raw) == status$expected$cutoff_rows,
              nrow(summary) == status$expected$summary_rows,
              nrow(comparisons) == status$expected$comparison_rows,
              !anyDuplicated(paste(comparisons$id, comparisons$method)),
              setequal(raw$id,plan$id), all(raw$signature == manifest$signature_hash),
              all(raw$input_md5 == plan$input_md5[match(raw$id,plan$id)]),
              all(comparisons$decision %in% c("IN","OUT","UNRESOLVED")),
              all(raw$B == 1999L), all(raw$k == 1900L),
              all(raw$cutoff_lower >= 0), all(raw$cutoff_lower <= raw$cutoff_upper),
              all(raw$cutoff_used == raw$cutoff_upper), all(raw$gap >= 0), all(raw$max_draw_gap >= 0),
              all(raw$elapsed_seconds >= 0), all(raw$stage >= 0 & raw$stage == floor(raw$stage) & raw$stage <= cfg$policy$stages),
              identical(cfg$objective, "coverage"), all(raw$objective == "coverage"),
              identical(raw$coverage_resolved, raw$coverage_lower == raw$coverage_upper),
              all(raw$status == ifelse(raw$coverage_resolved,
                  "COVERAGE_DECISION_ATTAINED","COVERAGE_BUDGET_UNRESOLVED")),
              all(grepl("^[0-9a-f]{32}$",raw$draw_hash)), all(comparisons$cutoff >= 0),
              all(comparisons$numerically_certified == ifelse(comparisons$method == "WALD_S1_KDE_CENTERED",FALSE,
                  comparisons$decision != "UNRESOLVED")),
              identical(is.na(comparisons$numerical_coverage), comparisons$decision == "UNRESOLVED"),
              all(comparisons$numerical_coverage[comparisons$decision != "UNRESOLVED"] ==
                  as.integer(comparisons$decision[comparisons$decision != "UNRESOLVED"] == "IN")))
    check_summary(raw, summary, manifest$signature$config$reps)
    comp_summary <- audit_csv(file.path(output,"comparison_summary.csv"),
        c("family","method","assurance","evaluation_scope"),
        c("parameter","n","replications","resolved","unresolved","numerically_certified",
          "numerical_coverage_lower","numerical_coverage_upper"))
    check_report_inventory(comp_summary, unique(expected_comparisons[c("family","parameter","n","method","assurance","evaluation_scope")]),
        c("family","parameter","n","method"))
    key <- function(x) paste(x$family,x$parameter,x$n,x$method,sep="|")
    groups <- split(comparisons,key(comparisons))
    stopifnot(nrow(comp_summary)==status$expected$comparison_summary_rows,
              !anyDuplicated(key(comp_summary)),setequal(names(groups),key(comp_summary)))
    for (i in seq_len(nrow(comp_summary))) {
        z <- groups[[key(comp_summary)[i]]]; total <- nrow(z)
        stopifnot(total==manifest$signature$config$reps,comp_summary$replications[i]==total,
            comp_summary$resolved[i]==sum(z$decision!="UNRESOLVED"),
            comp_summary$unresolved[i]==sum(z$decision=="UNRESOLVED"),
            comp_summary$numerically_certified[i]==sum(z$numerically_certified),
            all(z$assurance==comp_summary$assurance[i]),all(z$evaluation_scope==comp_summary$evaluation_scope[i]),
            abs(comp_summary$numerical_coverage_lower[i]-mean(z$decision=="IN"))<5e-14,
            abs(comp_summary$numerical_coverage_upper[i]-mean(z$decision!="OUT"))<5e-14)
    }
    resolved <- sum(raw$coverage_lower == raw$coverage_upper)
    stopifnot(status$coverage_resolved == resolved,
              status$coverage_unresolved == nrow(raw)-resolved,
              status$attained == sum(raw$attained),status$unresolved == sum(!raw$attained),
              status$status == if (all(raw$coverage_resolved))
                  "COMPLETE_COVERAGE_DECISIONS_ATTAINED" else "COMPLETE_WITH_UNRESOLVED_COVERAGE_DECISIONS")
    list(status = "PASS", signature = manifest$signature_hash,
         tasks = nrow(plan), coverage_rows = nrow(raw), resolved = resolved,
         unresolved = nrow(raw)-resolved, comparison_rows=nrow(comparisons),
         comparison_unresolved=sum(comparisons$decision=='UNRESOLVED'), report_md5 = status$report_md5)
}

# Numerical completion requires every enabled decision and angle row.
certification_incomplete <- function(message) stop(structure(list(
    message=paste0('CERTIFICATION_INCOMPLETE: ',message),call=NULL),
    class=c('certification_incomplete','error','condition')))

require_certified_completion <- function(coverage, expected, ranking = NULL, angles = NULL) {
    if (!identical(coverage$status,'PASS') || coverage$tasks!=expected$tasks ||
        coverage$coverage_rows!=expected$cutoff_rows || coverage$resolved!=expected$cutoff_rows ||
        coverage$unresolved!=0L || coverage$comparison_rows!=expected$comparison_rows ||
        coverage$comparison_unresolved!=0L)
        certification_incomplete('all G/L and comparator decisions must be resolved; bounds are preserved.')
    if (!is.null(ranking) && (ranking$candidate_method_rows!=432*expected$tasks/27 ||
        ranking$unresolved_candidate_methods!=0L))
        certification_incomplete('finite-ranking candidate decisions remain unresolved.')
    if (!is.null(angles) && (angles$rows!=expected$length_rows || angles$attained!=expected$length_rows ||
        angles$wide!=0L || !identical(angles$status,'COMPLETE_ANGLE_PRECISION_ATTAINED')))
        certification_incomplete('all enabled angle rows must attain the frozen precision.')
    invisible(TRUE)
}

reproduction_file_sha256 <- function(paths) {
    stopifnot(is.character(paths), !anyDuplicated(paths), all(file.exists(paths)),
              !any(grepl("[\r\n]", paths)))
    if (!length(paths)) return(setNames(character(), character()))
    # Bound the argument list; full experiments contain tens of thousands of files.
    chunks <- split(seq_along(paths), ceiling(seq_along(paths)/128L))
    result <- unlist(lapply(chunks, function(ix) {
        lines <- system2("sha256sum", c("--", shQuote(paths[ix])), stdout = TRUE)
        stopifnot(is.null(attr(lines, "status")), length(lines) == length(ix),
                  all(grepl("^[0-9a-f]{64} [ *]", lines)))
        substr(lines, 1L, 64L)
    }), use.names = FALSE)
    setNames(result, paths)
}

reproduction_canonical_table <- function(data, keys, exclude = character()) {
    stopifnot(is.data.frame(data), !anyDuplicated(names(data)),
              all(keys %in% names(data)), !anyNA(data[, keys, drop = FALSE]),
              !anyDuplicated(data[, keys, drop = FALSE]), !any(keys %in% exclude))
    # Keep every unlisted column, including future numerical/assurance fields.
    data <- data[, setdiff(names(data), exclude), drop = FALSE]
    data <- data[do.call(order, c(unname(data[, keys, drop = FALSE]),
                                 list(na.last = TRUE, method = "radix"))), , drop = FALSE]
    rownames(data) <- NULL
    data
}

reproduction_csv_values <- function(path, keys, exclude = character()) {
    # Character import retains the exact exported decimal/hexadecimal strings;
    # no tolerance, decimal-to-double conversion, or NA coercion hides a change.
    data <- read.csv(path, colClasses = "character", na.strings = NULL,
                     stringsAsFactors = FALSE, check.names = FALSE)
    reproduction_canonical_table(data, keys, exclude)
}

reproduction_results_report <- function(output) {
    files <- c("figure-coverage.tex","table-ranking.tex","table-circle.tex","plot-style.tex")
    directory <- file.path(output,"results")
    status <- system2("python3", c(shQuote(file.path(output,"source","export_results.py")), "check", "--run",
        shQuote(output), "--output", shQuote(directory)))
    if (status != 0L) stop("Reference output verification failed.")
    paths <- file.path(directory,files)
    list(tracked=paths,tex_sha256=setNames(unname(reproduction_file_sha256(paths)),paste0("results/",files)))
}
