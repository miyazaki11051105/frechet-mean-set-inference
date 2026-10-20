# Explicit, read-only enumeration from saved finite-ranking certificates.
# This report does not calibrate, refine, or rewrite a source run.

ranking_report_paths <- function(directory) {
    names <- c("candidates", "counts", "summary", "histograms")
    c(setNames(file.path(directory, paste0(names, ".csv")), names),
      audit = file.path(directory, "audit.json"))
}

verify_ranking_candidate_report <- function(directory, identity = NULL) {
    paths <- ranking_report_paths(directory)
    manifest_path <- file.path(directory, "manifest.rds")
    if (!all(file.exists(c(paths, manifest_path)))) return(FALSE)
    tryCatch({
        manifest <- readRDS(manifest_path)
        stopifnot(identical(manifest$signature, object_hash(manifest$identity)))
        if (!is.null(identity)) stopifnot(identical(manifest$identity, identity))
        audit <- jsonlite::read_json(paths[["audit"]])
        stopifnot(identical(audit$status, "PASS"),
                  identical(audit$schema, "finite-ranking-candidates-v2"),
                  identical(audit$schema, manifest$identity$schema),
                  identical(audit$signature, manifest$signature),
                  identical(audit$parent_signature, manifest$identity$run_signature))
        hashes <- unlist(audit$report_md5, use.names = TRUE)
        artifact_paths <- paths[names(paths) != "audit"]
        stopifnot(!anyDuplicated(basename(names(hashes))),
                  setequal(basename(names(hashes)), basename(artifact_paths)))
        expected <- hashes[match(basename(artifact_paths), basename(names(hashes)))]
        stopifnot(identical(unname(tools::md5sum(artifact_paths)), unname(expected)))
        TRUE
    }, error = function(error) FALSE)
}

ranking_report_input <- function(output) {
    require_current_run(output)
    manifest <- readRDS(file.path(output, "manifest.rds"))
    status <- readRDS(file.path(output, "STATUS.rds"))
    stopifnot(identical(manifest$signature_hash, object_hash(manifest$signature)),
              identical(status$signature, manifest$signature_hash),
              identical(status$expected, manifest$signature$expected),
              grepl("^COMPLETE_", status$status))
    source <- manifest$signature$source_md5
    stopifnot(identical(unname(tools::md5sum(file.path(output, "source", names(source)))),
                        unname(source)))
    sources <- list(source)
    reports <- check_report_hashes(output, status$report_md5)
    stopifnot(setequal(basename(reports),
                      c("raw.csv", "summary.csv", "comparisons.csv", "comparison_summary.csv")))
    plan <- subset(manifest$signature$plan, family == "ranking")
    stopifnot(nrow(plan) > 0L, !anyDuplicated(plan$id),
              all(grepl("^[A-Za-z0-9_.-]+$", plan$id)))
    metrics <- manifest$signature$config$metrics
    stopifnot(identical(metrics, c("euclidean", "kendall", "footrule")),
              nrow(plan) == length(manifest$signature$config$n) * manifest$signature$config$reps)
    raw <- read.csv(file.path(output, "raw.csv"), stringsAsFactors = FALSE)
    raw <- subset(raw, family == "ranking")
    stopifnot(nrow(raw) == nrow(plan) * length(metrics) * 2L,
              !anyDuplicated(paste(raw$id, raw$geometry, raw$method)),
              all(raw$method %in% c("G", "L")), setequal(raw$id, plan$id),
              all(raw$signature == manifest$signature_hash),
              all(raw$input_md5 == plan$input_md5[match(raw$id, plan$id)]),
              all(raw$B == 1999L), all(raw$k == 1900L),
              all(raw$cutoff_used == raw$cutoff_upper))
    summary <- read.csv(file.path(output, "summary.csv"), stringsAsFactors = FALSE)
    check_summary(raw, subset(summary, family == "ranking"), manifest$signature$config$reps)
    list(manifest = manifest, status = status, plan = plan, raw = split(raw, raw$id),
         metrics = metrics, sources = sources)
}

ranking_report_part <- function(task, output, reference, part_path, identity_hash) {
    input_path <- file.path(output, "inputs", paste0(task$id, ".rds"))
    paths <- file.path(output, "certificates", paste0(task$id, "-", reference$metrics, ".rds"))
    hashes <- unname(tools::md5sum(c(input_path, paths)))
    stopifnot(!anyNA(hashes), identical(hashes[1L], task$input_md5))
    if (file.exists(part_path)) {
        saved <- readRDS(part_path)
        stopifnot(identical(saved$identity_hash, identity_hash),
                  identical(saved$source_hashes, hashes),
                  identical(saved$payload_hash, object_hash(saved$payload)))
        return(saved$payload)
    }
    input <- readRDS(input_path)
    stopifnot(identical(input$family, "ranking"), length(input$observations) == task$n,
              input$B == 1999L, input$alpha == .05)
    original_draws <- regenerate_draws(input)
    original <- reference$raw[[task$id]]
    candidates <- counts <- vector("list", length(reference$metrics))
    for (i in seq_along(reference$metrics)) {
        metric <- reference$metrics[i]
        checkpoint <- readRDS(paths[i])
        payload <- verify_checkpoint(checkpoint, reference$manifest$signature_hash, task$input_md5)
        fit <- payload$fit
        stopifnot(inherits(fit, "gl_certificate"), identical(fit$domain, "finite_ranking"),
                  identical(fit$backend$metric, metric),
                  identical(fit$backend$observations, input$observations),
                  identical(fit$original_draw_hash, input$original_draw_hash),
                  identical(fit$draw_hash, object_hash(fit$draws)),
                  identical(fit$draws, original_draws[fit$backend$seen, , drop = FALSE]),
                  fit$B == 1999L, fit$k == 1900L, fit$alpha == .05,
                  identical(fit$eps, input$tuning[[metric]]$eps),
                  identical(fit$tau, input$tuning[[metric]]$tau),
                  any(vapply(reference$sources, function(source)
                      identical(unname(fit$source_md5), unname(source[names(fit$source_md5)])), logical(1))))
        expected_backend <- ranking_discrete_backend(input$observations, metric)
        stopifnot(identical(fit$backend, expected_backend))
        rank <- cutoff_order_statistic(fit$bounds$lower, fit$bounds$upper, fit$k, fit$policy$cutoff_tolerance)
        stopifnot(nrow(fit$bounds$lower) == fit$B,
                  all(is.finite(fit$bounds$lower)), all(is.finite(fit$bounds$upper)),
                  identical(rank$lower, fit$cutoff$lower), identical(rank$upper, fit$cutoff$upper),
                  identical(rank$used, fit$cutoff$used))
        rows <- original[original$geometry == metric, , drop = FALSE]
        rows <- rows[match(c("G", "L"), rows$method), , drop = FALSE]
        stopifnot(nrow(rows) == 2L, !anyNA(rows$method),
                  identical(rows$draw_hash, rep(fit$draw_hash, 2L)),
                  all(rows$cutoff_lower == fit$cutoff$lower),
                  all(rows$cutoff_upper == fit$cutoff$upper),
                  all(rows$coverage_lower == payload$evaluation$coverage_lower),
                  all(rows$coverage_upper == payload$evaluation$coverage_upper))
        result <- ranking_candidate_sets(fit)
        z <- result$candidates
        stopifnot(nrow(z) == 24L, identical(z$candidate_id, seq_len(24L)),
                  identical(z$ranking_key, ranking_discrete_geometry(metric)$keys),
                  identical(result$counts$set,
                    c("G", "L", "G_only", "L_only", "intersection", "symmetric_difference")),
                  all(result$counts$lower >= 0L), all(result$counts$upper <= 24L),
                  all(result$counts$lower <= result$counts$upper),
                  all(result$counts$unresolved == result$counts$upper - result$counts$lower))
        base <- data.frame(id = task$id, n = task$n, replication = task$replication,
                           geometry = metric, input_md5 = task$input_md5,
                           draw_hash = fit$draw_hash, original_draw_hash = fit$original_draw_hash,
                           certificate_md5 = hashes[i + 1L], payload_hash = checkpoint$payload_hash)
        candidates[[i]] <- do.call(rbind, lapply(c("G", "L"), function(method) {
            fields <- c("lower_cutoff_decision", "upper_cutoff_decision",
                        "lower_cutoff_raw_decision", "upper_cutoff_raw_decision",
                        "lower_cutoff_certified", "upper_cutoff_certified", "decision", "lower", "upper")
            zz <- z[, paste0(method, "_", fields), drop = FALSE]
            names(zz) <- fields
            cbind(base[rep(1L, 24L), ], method = method,
                  candidate_id = z$candidate_id, ranking_key = z$ranking_key,
                  cutoff_lower = unname(fit$cutoff$lower[match(method, c("G", "L"))]),
                  cutoff_upper = unname(fit$cutoff$upper[match(method, c("G", "L"))]), zz, row.names = NULL)
        }))
        count <- base
        for (j in seq_len(nrow(result$counts))) {
            set <- result$counts$set[j]
            for (field in c("lower", "upper", "unresolved", "resolved"))
                count[[paste0(set, "_", field)]] <- result$counts[[field]][j]
        }
        count$G_minus_L_lower <- count$G_only_lower - count$L_only_upper
        count$G_minus_L_upper <- count$G_only_upper - count$L_only_lower
        count$equal_certified <- count$symmetric_difference_upper == 0L
        count$equal_possible <- count$symmetric_difference_lower == 0L
        count$L_strict_subset_G_certified <- count$L_only_upper == 0L && count$G_only_lower > 0L
        count$L_strict_subset_G_possible <- count$L_only_lower == 0L && count$G_only_upper > 0L
        calibration <- result$calibration_candidates
        if (!is.null(calibration)) stopifnot(identical(calibration$method, c("G", "L")),
            all(calibration$total_ordered_pairs == 552L),
            all(calibration$feasible_ordered_pairs >= 0L),
            all(calibration$feasible_ordered_pairs <= 552L),
            all(calibration$first_endpoint_candidates >= 0L),
            all(calibration$first_endpoint_candidates <= 24L))
        for (method in c("G", "L")) {
            j <- match(method, c("G", "L"))
            count[[paste0(method, "_calibration_pairs")]] <-
                if (is.null(calibration)) NA_integer_ else calibration$feasible_ordered_pairs[j]
            count[[paste0(method, "_calibration_first_endpoints")]] <-
                if (is.null(calibration)) NA_integer_ else calibration$first_endpoint_candidates[j]
        }
        # Endpoint results describe the implemented upper-cutoff set separately.
        # Joint G/L nesting above applies to the ideal enclosed cutoffs.
        for (method in c("G", "L")) {
            decision <- z[[paste0(method, "_upper_cutoff_decision")]]
            certified <- z[[paste0(method, "_upper_cutoff_certified")]]
            count[[paste0(method, "_upper_cutoff_count_lower")]] <- sum(certified & decision == "IN")
            count[[paste0(method, "_upper_cutoff_count_upper")]] <- sum(!(certified & decision == "OUT"))
        }
        count$nesting <- result$nesting
        count$assurance <- result$assurance
        counts[[i]] <- count
    }
    stopifnot(identical(unname(tools::md5sum(c(input_path, paths))), hashes))
    result <- list(candidates = do.call(rbind, candidates), counts = do.call(rbind, counts))
    save_atomic(list(identity_hash = identity_hash, source_hashes = hashes,
                     payload = result, payload_hash = object_hash(result)), part_path)
    result
}

summarize_ranking_candidates <- function(counts) {
    groups <- split(counts, paste(counts$geometry, counts$n, sep = "|"))
    sets <- c("G", "L", "G_only", "L_only", "intersection", "symmetric_difference", "G_minus_L")
    summary <- histogram <- list()
    for (key in names(groups)) {
        z <- groups[[key]]
        row <- data.frame(geometry = z$geometry[1L], n = z$n[1L], replications = nrow(z))
        for (set in sets) {
            lower <- z[[paste0(set, "_lower")]]
            upper <- z[[paste0(set, "_upper")]]
            row[[paste0(set, "_mean_lower")]] <- mean(lower)
            row[[paste0(set, "_mean_upper")]] <- mean(upper)
            row[[paste0(set, "_minimum_lower")]] <- min(lower)
            row[[paste0(set, "_minimum_upper")]] <- min(upper)
            row[[paste0(set, "_maximum_lower")]] <- max(lower)
            row[[paste0(set, "_maximum_upper")]] <- max(upper)
            row[[paste0(set, "_count_unresolved")]] <- sum(lower != upper)
            values <- if (set == "G_minus_L") -24:24 else 0:24
            for (value in values) {
                certain <- sum(lower == value & upper == value)
                possible <- sum(lower <= value & value <= upper)
                histogram[[length(histogram) + 1L]] <- data.frame(
                    geometry = z$geometry[1L], n = z$n[1L], set = set, count = value,
                    frequency_lower = certain, frequency_upper = possible,
                    replications = nrow(z))
            }
        }
        for (field in c("equal_certified", "equal_possible", "L_strict_subset_G_certified",
                        "L_strict_subset_G_possible")) row[[field]] <- sum(z[[field]])
        for (method in c("G", "L")) {
            lo <- z[[paste0(method, "_upper_cutoff_count_lower")]]
            hi <- z[[paste0(method, "_upper_cutoff_count_upper")]]
            row[[paste0(method, "_upper_cutoff_mean_lower")]] <- mean(lo)
            row[[paste0(method, "_upper_cutoff_mean_upper")]] <- mean(hi)
            row[[paste0(method, "_upper_cutoff_unresolved")]] <- sum(lo != hi)
            for (quantity in c("calibration_pairs", "calibration_first_endpoints")) {
                value <- z[[paste0(method, "_", quantity)]]
                present <- value[!is.na(value)]
                row[[paste0(method, "_", quantity, "_missing")]] <- sum(is.na(value))
                row[[paste0(method, "_", quantity, "_mean")]] <- if (length(present)) mean(present) else NA_real_
                row[[paste0(method, "_", quantity, "_minimum")]] <- if (length(present)) min(present) else NA_real_
                row[[paste0(method, "_", quantity, "_maximum")]] <- if (length(present)) max(present) else NA_real_
            }
        }
        summary[[key]] <- row
    }
    list(summary = do.call(rbind, summary), histograms = do.call(rbind, histogram))
}

write_ranking_candidate_report <- function(output, destination, workers = 1L, batch_size = 24L) {
    stopifnot(length(workers) == 1L, is.finite(workers), workers == as.integer(workers),
              workers >= 1L, workers <= 4L, length(batch_size) == 1L,
              is.finite(batch_size), batch_size >= 1L, batch_size == as.integer(batch_size))
    if (.Platform$OS.type != "unix") stop("Use the tested Linux/WSL runtime for saved-certificate analysis")
    output <- normalizePath(output, mustWork = TRUE)
    destination <- path.expand(destination)
    if (!startsWith(destination, "/")) destination <- file.path(getwd(), destination)
    tail <- character()
    parent <- destination
    while (!dir.exists(parent)) {
        stopifnot(!file.exists(parent), dirname(parent) != parent)
        tail <- c(basename(parent), tail)
        parent <- dirname(parent)
    }
    destination <- do.call(file.path, as.list(c(normalizePath(parent, mustWork = TRUE), tail)))
    if (any(strsplit(destination, "/", fixed = TRUE)[[1L]] == ".."))
        stop("Use a report destination without parent-directory components")
    if (identical(destination, output) || startsWith(destination, paste0(output, "/")))
        stop("Use a new report directory outside the read-only source run")
    if (!dir.exists(destination)) dir.create(destination, recursive = TRUE)
    destination <- normalizePath(destination, mustWork = TRUE)
    lock <- acquire_lock(file.path(destination, "analysis.lock"), seconds = 0)
    on.exit(try(lock$kill(), silent = TRUE), add = TRUE)
    reference <- ranking_report_input(output)
    source_files <- c("ranking/report.R", "ranking/geometry.R", "ranking/numerics.hpp",
                      "common/inference.R", "common/reproduce.R", "common/numerics.cpp",
                      "common/interval.hpp", "common/exact.py")
    source_hashes <- tools::md5sum(source_files)
    stopifnot(!anyNA(source_hashes))
    identity <- list(schema = "finite-ranking-candidates-v2",
                     run_signature = reference$manifest$signature_hash,
                     report_md5 = reference$status$report_md5,
                     plan = reference$plan, metrics = reference$metrics,
                     analysis_source_md5 = source_hashes, runtime = runtime_versions())
    identity_hash <- object_hash(identity)
    manifest_path <- file.path(destination, "manifest.rds")
    if (file.exists(manifest_path)) {
        saved <- readRDS(manifest_path)
        stopifnot(identical(saved$identity, identity), identical(saved$signature, identity_hash))
    } else {
        if (length(setdiff(list.files(destination, all.files = TRUE, no.. = TRUE), "analysis.lock")))
            stop("Destination is not an empty or compatible ranking analysis directory")
        save_atomic(list(identity = identity, signature = identity_hash), manifest_path)
    }
    if (verify_ranking_candidate_report(destination, identity)) {
        return(invisible(ranking_report_paths(destination)))
    }
    load_numerics()
    directory <- file.path(destination, "parts")
    dir.create(directory, showWarnings = FALSE)
    plan <- reference$plan
    all_parts <- vector("list", nrow(plan))
    started <- proc.time()[["elapsed"]]
    for (first in seq.int(1L, nrow(plan), by = batch_size)) {
        indexes <- seq.int(first, min(nrow(plan), first + batch_size - 1L))
        results <- parallel::mclapply(indexes, function(i) {
            task <- lapply(plan[i, , drop = FALSE], `[[`, 1L)
            ranking_report_part(task, output, reference,
                file.path(directory, paste0(task$id, ".rds")), identity_hash)
        }, mc.cores = min(as.integer(workers), length(indexes)), mc.preschedule = FALSE,
           mc.set.seed = FALSE)
        failed <- vapply(results, inherits, logical(1), "try-error")
        if (any(failed)) stop("Ranking candidate analysis failed; completed report parts retained: ",
                              paste(as.character(results[failed]), collapse = "\n"))
        all_parts[indexes] <- results
    }
    candidates <- do.call(rbind, lapply(all_parts, `[[`, "candidates"))
    counts <- do.call(rbind, lapply(all_parts, `[[`, "counts"))
    stopifnot(nrow(candidates) == nrow(plan) * length(reference$metrics) * 48L,
              nrow(counts) == nrow(plan) * length(reference$metrics),
              !anyDuplicated(paste(candidates$id, candidates$geometry, candidates$method, candidates$candidate_id)),
              !anyDuplicated(paste(counts$id, counts$geometry)),
              identical(source_hashes, tools::md5sum(source_files)))
    # Recheck the historical reports after analysis; no output below enters that run.
    check_report_hashes(output, reference$status$report_md5)
    summaries <- summarize_ranking_candidates(counts)
    paths <- ranking_report_paths(destination)
    reports <- c(list(candidates = candidates, counts = counts), summaries)
    for (name in names(reports)) save_atomic(reports[[name]], paths[[name]])
    audit <- list(status = "PASS", schema = identity$schema, signature = identity_hash,
        parent_signature = identity$run_signature, datasets = nrow(plan), fits = nrow(counts),
        candidate_method_rows = nrow(candidates), summary_groups = nrow(summaries$summary),
        unresolved_candidate_methods = sum(candidates$lower != candidates$upper),
        unresolved_G_counts = sum(counts$G_lower != counts$G_upper),
        unresolved_L_counts = sum(counts$L_lower != counts$L_upper),
        missing_calibration_selector_fits = sum(is.na(counts$G_calibration_pairs) | is.na(counts$L_calibration_pairs)),
        original_report_md5 = as.list(reference$status$report_md5),
        analysis_source_md5 = as.list(source_hashes),
        report_md5 = as.list(tools::md5sum(paths[names(paths) != "audit"])),
        assurance = "EXHAUSTIVE_24_CANDIDATES_FROM_SAVED_CALIBRATION_ENCLOSURES_NO_NEW_CALIBRATION",
        uncertainty = "Deterministic enclosures conditional on the saved binary64 inputs and B draws; no simultaneous statistical confidence claim.",
        histograms = "Frequency bounds are pointwise across count values; upper frequencies need not sum to the number of replications.",
        elapsed_seconds = proc.time()[["elapsed"]] - started)
    save_atomic(audit, paths[["audit"]])
    invisible(paths)
}
