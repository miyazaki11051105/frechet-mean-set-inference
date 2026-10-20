# Shared calibration, confidence-set interface and reproducible execution.
elapsed_time <- function() nanonext::mclock() / 1000

# Probe completed files before targets hashes them on the Windows-mounted drive.
# A failed probe is an execution error, never a numerical decision.
ready_files <- function(paths, attempts = 5L) {
    stopifnot(is.character(paths), length(paths) > 0L, attempts >= 1L)
    for (i in seq_len(attempts)) {
        hashes <- suppressWarnings(tools::md5sum(paths))
        if (all(file.exists(paths)) && !anyNA(hashes)) return(paths)
        Sys.sleep(.1 * i)
    }
    stop("Completed output files are not readable: ", paste(paths[is.na(hashes)], collapse = ", "))
}

task_attempt <- function(task, output, signature, complete = FALSE) {
    controller <- readRDS(file.path(output, "execution-attempts", "controller.rds"))
    stopifnot(is.character(controller$epoch), length(controller$epoch) == 1L, nzchar(controller$epoch),
              identical(controller$signature, signature))
    path <- file.path(output, "execution-attempts", paste0("coverage-", task$id, ".rds"))
    if (complete) {
        save_atomic(list(epoch = controller$epoch, fingerprint = NULL, count = 0L, complete = TRUE), path)
        return(invisible(path))
    }
    fits <- if (task$family == "ranking") c("euclidean", "kendall", "footrule") else task$family
    files <- file.path(output, "certificates", paste0(task$id, "-", fits, ".rds"))
    files <- files[file.exists(files)]
    fingerprint <- object_hash(list(input = task$input_md5,
        source = tools::md5sum(c("_targets.R", "common/inference.R")), checkpoints = tools::md5sum(files)))
    prior <- if (file.exists(path)) readRDS(path) else NULL
    count <- if (!is.null(prior) && identical(prior$epoch, controller$epoch) &&
        !isTRUE(prior$complete) && identical(prior$fingerprint, fingerprint)) prior$count else 0L
    if (!complete) count <- count + 1L
    save_atomic(list(epoch = controller$epoch, fingerprint = fingerprint, count = count, complete = complete), path)
    if (count > 3L) stop("Repeated worker loss without checkpoint progress: ", task$id,
        "; repair the cause before resuming. Existing results are retained.")
    invisible(path)
}

begin_execution_epoch <- function(output, signature) {
    record <- list(epoch = basename(tempfile("controller-")), signature = signature,
                   pid = Sys.getpid(), started = as.character(Sys.time()))
    save_atomic(record, file.path(output, "execution-attempts", "controller.rds"))
    invisible(record)
}

reproduction_files <- function() {
    c("run.R", "setup.R", "_targets.R", "plot.R", "coverage-layout.tex", "export_results.py", "verify_release.py",
      "reference/coverage.csv", "reference/ranking.csv", "reference/circle.csv", "reference/figure-coverage.template.tex",
      "reference/table-ranking.template.tex", "reference/table-circle.template.tex", "reference/manifest.json", "reference/plot-style.tex", "settings.dcf",
      "environment/runtime.dcf", "environment/packages.dcf", "common/reproduce.R", "ranking/report.R", "circle/length.R", "circle/length.hpp",
      "common/inference.R", "common/numerics.cpp", "common/interval.hpp",
      "common/exact.py",
      unlist(lapply(c("euclid", "ranking", "circle"), function(geometry)
          file.path(geometry, c("geometry.R", "numerics.hpp"))), use.names = FALSE))
}

calibrate_gl <- function(backend = NULL, Y = NULL, draws, alpha = 0.05, eps = 0.1, tau = 2, policy = refinement_policy(),
    witnesses = NULL, coarse = FALSE) {
    started <- elapsed_time()
    if (is.null(backend) == is.null(Y))
        stop("supply exactly one observed backend or Y")
    check_alpha(alpha)
    if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps < 0 || !is.numeric(tau) || length(tau) !=
        1L || !is.finite(tau) || tau <= 0 || !is.matrix(draws) || !is.numeric(draws) || ncol(draws) < 1L || any(!is.finite(draws)) ||
        !identical(policy, do.call(refinement_policy, policy)) || !is.logical(coarse) || length(coarse) != 1L || is.na(coarse))
        stop("invalid fixed-input G/L controls")
    load_numerics()
    B <- ncol(draws)
    load_seconds <- elapsed_time() - started
    k <- as.integer(cutoff_index_cpp(B, alpha))
    if (is.na(k))
        k <- exact_cutoff_index(B, alpha)
    if (k > B)
        stop("add-one rank exceeds B; increase B or use an explicit infinite cutoff")
    if (!is.null(Y)) {
        euclid <- euclid_bounds(Y, draws, witnesses)
        Y <- euclid$Y
        bounds <- euclid$bounds
        metadata <- euclid$metadata
        domain <- euclid$domain
        identity <- euclid$identity
        trace <- euclid$trace
    }
    else {
        domain <- backend$domain
        if (!domain %in% c("full_circle", "continuous_ranking_graph", "finite_ranking") || nrow(draws) != length(backend$counts))
            stop("unsupported backend or draw row mismatch")
        identity <- backend[intersect(c("domain", "counts", "seen", "raw_observations", "observations", "pieces",
            "loss", "pairs", "metric"), names(backend))]
        if (is.function(identity$loss))
            identity$loss <- NULL
        if (domain == "finite_ranking") {
            if (!is.null(witnesses))
                stop("finite ranking enumerates all comparisons")
            raw <- ranking_finite_bounds(backend, draws, eps, tau)
            bounds <- raw[c("lower", "upper")]
            metadata <- raw[setdiff(names(raw), c("lower", "upper"))]
            trace <- list()
        }
        else {
            input <- if (domain == "full_circle")
                list(pieces = circle_charts_cpp(backend$raw_observations), counts = backend$counts)
            else piecewise_input(backend)
            if (coarse && domain == "full_circle")
                witnesses <- unique(rbind(witnesses, circle_witnesses(backend, eps)))
            points <- witness_coordinates(backend, input, witnesses)
            seeds <- matrix(0, B, 2L)
            seed_metadata <- NULL
            if (nrow(points) || coarse) {
                model <- build_piecewise_model_cpp(input$pieces, input$counts, eps, tau, TRUE, TRUE, TRUE)
                tryCatch({
                  seed_metadata <- evaluate_witnesses_cpp(model, input$pieces, input$counts, points, draws, eps,
                    tau)
                  seeds <- seed_metadata$lower[, 1:2, drop = FALSE]
                }, finally = release_piecewise_model_cpp(model))
            }
            trace <- list()
            bounds <- NULL
            if (coarse) {
                bounds <- list(lower = seeds, upper = cbind(seed_metadata$global_cap, seed_metadata$global_cap))
                rank <- cutoff_order_statistic(bounds$lower, bounds$upper, k, policy$cutoff_tolerance)
                cover <- list(fmin = seed_metadata$fmin, solver_state = NULL)
                trace <- list(list(stage = 0L, cutoff_lower = rank$lower, cutoff_upper = rank$upper,
                    gap = rank$gap, scope = "FEASIBLE_WITNESSES_AND_FULL_DOMAIN_GAUSSIAN_NORM_CAP"))
            } else {
                workspace <- build_piecewise_model_cpp(input$pieces, input$counts, eps, tau,
                    TRUE, TRUE, TRUE, draws)
                on.exit(release_piecewise_model_cpp(workspace), add = TRUE)
            }
            for (stage in if (coarse) integer() else seq_len(policy$stages)) {
                if (stage == 1L) {
                  active <- rep(TRUE, B)
                  methods <- c(TRUE, TRUE)
                }
                else {
                  rank <- cutoff_order_statistic(bounds$lower, bounds$upper, k, policy$cutoff_tolerance)
                  methods <- !rank$precision_met
                  active <- rep(FALSE, B)
                  for (j in which(!rank$precision_met)) active <- active | (bounds$upper[, j] > rank$lower[j] &
                    bounds$lower[, j] <= rank$upper[j])
                  if (!policy$rank_refinement)
                    active[] <- TRUE
                  if (!any(active))
                    break
                }
                lo <- hi <- matrix(0, B, 2L)
                if (!is.null(bounds)) {
                  lo[, 1:2] <- bounds$lower
                  hi[, 1:2] <- bounds$upper
                  lo[active, which(methods)] <- hi[active, which(methods)] <- 0
                }
                args <- list(g = draws, eps = eps, tau = tau, max_cells = as.integer(policy$max_cells *
                  2^(stage - 1L)), max_depth = min(30L, policy$max_depth + stage - 1L), direction_tolerance = policy$direction_tolerance/2^(stage -
                  1L), continuous_cones = TRUE, circle_prefix = TRUE, radial_feasibility = TRUE, sharp_rectangles = TRUE,
                  global_incumbent = TRUE, skip_inactive_chunks = TRUE, pair_batch = policy$pair_batch,
                  cutoff_rank = if (policy$rank_refinement) k else NULL, cutoff_tolerance = policy$cutoff_tolerance,
                  model = workspace)
                cover <- piecewise_bounds(input, args, B, policy$draw_batch, refinement = list(methods = methods,
                  active = active, global_lower = if (is.null(bounds)) seeds else bounds$lower, bounds = list(lower = lo, upper = hi, active_methods = integer(2L))))
                next_lower <- pmax(seeds, cover$lower[, 1:2, drop = FALSE])
                next_upper <- cover$upper[, 1:2, drop = FALSE]
                if (!is.null(bounds)) {
                  next_lower <- pmax(bounds$lower, next_lower)
                  next_upper <- pmin(bounds$upper, next_upper)
                }
                bounds <- intersect_nested_bounds(next_lower, next_upper)
                rank <- cutoff_order_statistic(bounds$lower, bounds$upper, k, policy$cutoff_tolerance)
                trace[[stage]] <- list(stage = stage, active_draws = sum(active), cutoff_lower = rank$lower, cutoff_upper = rank$upper,
                  gap = rank$gap, max_cells = args$max_cells, max_depth = args$max_depth, metadata = cover[c("visited_cells",
                    "budget_leaves", "blocks", "total_pairs", "active_methods")])
                if (all(rank$precision_met))
                  break
            }
            metadata <- list(witnesses = points, seed_metadata = seed_metadata, total_pairs = length(input$pieces)^2,
                minimum = cover$fmin, solver_state = if (!all(rank$precision_met)) cover$solver_state else NULL)
        }
    }
    bounds <- intersect_nested_bounds(bounds$lower, bounds$upper)
    colnames(bounds$lower) <- colnames(bounds$upper) <- c("G", "L")
    cutoff <- cutoff_order_statistic(bounds$lower, bounds$upper, k, policy$cutoff_tolerance)
    cutoff$status <- ifelse(cutoff$precision_met, "CUTOFF_GAP_ATTAINED", "CUTOFF_BUDGET_UNRESOLVED")
    files <- reproduction_files()
    structure(list(schema = "fixed-input-certified-gl-v1", domain = domain, Y = Y, backend = backend, draws = draws,
        input_hash = object_hash(identity), draw_hash = object_hash(draws), B = B, k = k, alpha = alpha, eps = eps,
        tau = tau, policy = policy, bounds = bounds, cutoff = cutoff, trace = trace, metadata = metadata, source_md5 = tools::md5sum(files),
        timing = list(load_seconds = load_seconds, calibration_seconds = elapsed_time() - started - load_seconds,
            total_seconds = elapsed_time() - started), assurance = "DIRECTED_INTERVAL_FIXED_BINARY64_INPUT_AND_FIXED_B_ORDER_STATISTIC",
        mc_assurance = "SUPPLIED_DRAWS_NO_IID_CLAIM_NO_BOOTSTRAP_APPROXIMATION_GUARANTEE"), class = "gl_certificate")
}

check_alpha <- function(alpha) {
    if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 0.5)
        stop("alpha must lie in (0, 0.5)")
    invisible(alpha)
}

check_counts <- function(counts) {
    if (!is.numeric(counts) || !length(counts) || any(!is.finite(counts)) || any(counts <= 0) || any(counts !=
        floor(counts)))
        stop("counts must be positive integer-valued numbers for observed objects")
    if (sum(counts) > 2^52)
        stop("count total exceeds exact integer envelope")
    invisible(counts)
}

check_sample_size <- function(n) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n > .Machine$integer.max || n != floor(n))
        stop("n must be a positive integer")
    invisible(n)
}

save_atomic <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(paste0(".publish-", basename(path), "-"), tmpdir = dirname(path))
    on.exit(unlink(tmp), add = TRUE)
    if (grepl("[.]csv$", path)) {
        value[] <- lapply(value, function(x) if (is.numeric(x))
            sprintf("%.17g", x)
        else x)
        write.csv(value, tmp, row.names = FALSE)
    }
    else if (grepl("[.]json$", path)) jsonlite::write_json(value, tmp, auto_unbox = TRUE, null = "null", digits = 17)
    else saveRDS(value, tmp, version = 3L)
    stopifnot(file.rename(tmp, path))
    invisible(path)
}

# Call only while holding the destination's controller lock. Preserve abandoned
# publications; legacy anonymous files require an exact match to the intended RDS.
recover_atomic_publications <- function(directory, targets, legacy_values = list(), archive_root = directory) {
    if (!dir.exists(directory)) return(invisible(NULL))
    paths <- list.files(directory, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    paths <- paths[!dir.exists(paths)]
    owned <- vapply(paths, function(path) {
        name <- basename(path)
        for (target in targets) {
            prefix <- paste0(".publish-", target, "-")
            if (startsWith(name, prefix) && grepl("^[[:xdigit:]]+$", substring(name, nchar(prefix) + 1L))) return(TRUE)
        }
        if (!grepl("^[.]publish-[[:xdigit:]]+$", name)) return(FALSE)
        value <- tryCatch(readRDS(path), error = function(e) NULL)
        !is.null(value) && any(vapply(legacy_values, identical, logical(1), y = value))
    }, logical(1))
    if (!any(owned)) return(invisible(NULL))
    archive <- tempfile(paste0(basename(archive_root), "-publication-incomplete-"), tmpdir = dirname(archive_root))
    stopifnot(dir.create(archive), all(file.rename(paths[owned], file.path(archive, basename(paths[owned])))))
    invisible(archive)
}

freeze_source <- function(files, hashes, directory, archive_root) {
    for (f in files) {
        to <- file.path(directory, f)
        dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
        recover_atomic_publications(dirname(to), basename(to), archive_root = archive_root)
        if (file.exists(to)) stopifnot(identical(unname(tools::md5sum(to)), unname(hashes[f])))
        else {
            tmp <- tempfile(paste0(".publish-", basename(to), "-"), tmpdir = dirname(to))
            tryCatch({
                stopifnot(file.copy(f, tmp, overwrite = FALSE),
                    identical(unname(tools::md5sum(tmp)), unname(hashes[f])), file.rename(tmp, to))
            }, finally = unlink(tmp))
        }
    }
    stopifnot(identical(unname(tools::md5sum(file.path(directory, files))), unname(hashes)))
    invisible(directory)
}

verify_cutoff_gap <- function(f, decimal_gap) {
    fmt <- function(x) sprintf("%.17g", x)
    encode <- function(x) matrix(fmt(x), nrow(x), ncol(x))
    proof <- list(lower = encode(f$bounds$lower), upper = encode(f$bounds$upper),
        cutoff = lapply(f$cutoff[c("lower", "upper", "used", "gap")], fmt), B = f$B, k = f$k, alpha = fmt(f$alpha),
        requested_decimal_gap = decimal_gap)
    path <- tempfile(fileext = ".json")
    log <- tempfile()
    on.exit(unlink(c(path, log)), add = TRUE)
    jsonlite::write_json(proof, path, digits = 17, auto_unbox = TRUE, matrix = "rowmajor")
    status <- system2("python3", c("common/exact.py", "cutoff", shQuote(path)), stdout = log, stderr = log)
    if (!status %in% c(0L, 2L))
        stop("independent all-draw rank/decimal audit failed: ", paste(readLines(log), collapse = " "))
    list(attained = status == 0L, log = readLines(log, warn = FALSE))
}

verify_checkpoint <- function(z, signature, input_md5) {
    stopifnot(identical(z$signature, signature), identical(z$input_md5, input_md5), identical(z$payload_hash, object_hash(z$payload)))
    z$payload
}

evaluate_comparisons <- function(input, task, f) {
    if (!input$family %in% c("euclidean", "circle"))
        return(NULL)
    if (input$family == "euclidean") {
        r <- f$metadata$exact_rank
        n <- nrow(f$Y)
        cuts <- if (r == 0L)
            c(0, 0)
        else c(sqrt(qchisq(1 - input$alpha, r)), sqrt(r * n/(n - r) * qf(1 - input$alpha, r, n - r)))
        decisions <- vapply(cuts, function(cut) euclid_exact(f$Y, "query", rep(0, ncol(f$Y)), cut)$decision, character(1))
        methods <- c("WALD_CHI", "WALD_F_GAUSSIAN_REFERENCE")
        certified <- rep(TRUE, 2)
        lower <- upper <- rep(NA_real_, 2)
        assurance <- rep("EXACT_RATIONAL_INVERSION_FIXED_BINARY64_THRESHOLD", 2)
        center <- radius <- curvature <- score_variance <- rep(NA_real_, 2)
    }
    else {
        load_numerics()
        points <- circle_population_means(task$parameter)$means
        intervals <- circle_population_mean_bounds_cpp(task$parameter)
        el <- circle_el_bounds_cpp(input$observations, intervals[, 1L], qchisq(1 - input$alpha, 1), intervals[, 2L])
        e <- if (any(el$decision == "OUT"))
            "OUT"
        else if (all(el$decision == "IN"))
            "IN"
        else "UNRESOLVED"
        w <- lapply(points, function(q) circle_wald(input$observations, q, input$alpha))
        wd <- vapply(w, `[[`, character(1), "decision")
        decisions <- c(e, if (any(wd == "OUT")) "OUT" else if (all(wd == "IN")) "IN" else "UNRESOLVED")
        methods <- c("EL_S1_RAW", "WALD_S1_KDE_CENTERED")
        cuts <- rep(qchisq(1 - input$alpha, 1), 2)
        certified <- c(e != "UNRESOLVED", FALSE)
        lower <- c(max(el$lower), NA_real_)
        upper <- c(max(el$upper), NA_real_)
        assurance <- c(el$assurance, "FLOAT64_KDE_WALD_REGULARITY_NOT_CERTIFIED")
        center <- radius <- curvature <- score_variance <- rep(NA_real_, 2)
        if (w[[1]]$status == "WALD_COMPUTED") {
            center[2] <- w[[1]]$mean
            curvature[2] <- w[[1]]$curvature
            score_variance[2] <- w[[1]]$score_variance
            radius[2] <- min(pi, sqrt(cuts[2] * score_variance[2]/(length(input$observations) * curvature[2]^2)))
        }
    }
    data.frame(id = task$id, family = task$family, parameter = task$parameter, n = task$n, replication = task$replication,
        method = methods, cutoff = cuts, decision = decisions, statistic_lower = lower, statistic_upper = upper,
        numerically_certified = certified, center = center, radius = radius, curvature = curvature,
        score_variance = score_variance, numerical_coverage = ifelse(decisions == "UNRESOLVED", NA, as.integer(decisions ==
            "IN")), assurance = assurance, evaluation_scope = if (task$family == "circle" && task$parameter < 0)
            c("WHOLE_POPULATION_MEAN_SET_WITH_CERTIFIED_QUERY_ENCLOSURES",
              "BINARY64_REPRESENTATIVES_OF_NONBINARY_POPULATION_MEANS")
        else "WHOLE_POPULATION_MEAN_SET")
}

decision_refinement_controls <- function(budgets = c(262144,1048576,4194304,16777216),
    max_depth = 60L, max_draws = 32L, after_stage = 1L, fraction = .25) {
    integer <- function(x,lo,hi) is.numeric(x) && length(x)==1L && is.finite(x) &&
        x==floor(x) && x>=lo && x<=hi
    if (!is.numeric(budgets) || !length(budgets) || any(!is.finite(budgets)) ||
        any(budgets!=floor(budgets)) || any(budgets<1 | budgets>.Machine$integer.max) ||
        is.unsorted(budgets,strictly=TRUE) || !integer(max_depth,0,60) ||
        !integer(max_draws,1,1999) || !integer(after_stage,1,20) ||
        !is.numeric(fraction) || length(fraction)!=1L || !is.finite(fraction) ||
        fraction<=0 || fraction>=1) stop('Invalid fixed decision-refinement controls')
    list(schema='OBSERVED_ONLY_LOCAL_RANK_V1',budgets=as.integer(budgets),
        max_depth=as.integer(max_depth),max_draws=as.integer(max_draws),
        after_stage=as.integer(after_stage),fraction=fraction)
}

# No task, population parameter or mean-set oracle enters this numerical step.
# Each root pair is covered; the kernel retains all discarded/unvisited uppers.
refine_local_decision_rank <- function(fit, controls, budget) {
    expected <- do.call(decision_refinement_controls,controls[setdiff(names(controls),'schema')])
    if (!identical(controls,expected) || !inherits(fit,'gl_certificate') ||
        !fit$domain %in% c('full_circle','continuous_ranking_graph') ||
        !identical(fit$assurance,'DIRECTED_INTERVAL_FIXED_BINARY64_INPUT_AND_FIXED_B_ORDER_STATISTIC') ||
        !identical(object_hash(fit$draws),fit$draw_hash) || !budget %in% controls$budgets)
        stop('Trusted fixed-input certificate and frozen rank controls required')
    b <- fit$backend
    identity <- b[intersect(c('domain','counts','seen','raw_observations','observations','pieces',
        'loss','pairs','metric'),names(b))]
    if (is.function(identity$loss)) identity$loss <- NULL
    stopifnot(identical(object_hash(identity),fit$input_hash))
    lo <- fit$cutoff$lower[2L]; hi <- fit$cutoff$upper[2L]
    target <- unname(lo+(hi-lo)*controls$fraction)
    # No representable interior threshold: preserve unresolved, never round a tie.
    if (!is.finite(target) || target<=lo || target>=hi) return(NULL)
    ids <- which(fit$bounds$lower[,2L]<=target & fit$bounds$upper[,2L]>target)
    roots <- if (fit$domain=='full_circle') (length(b$raw_observations)+1L)^2 else length(b$pieces)^2
    if (!length(ids) || length(ids)>controls$max_draws || budget<roots) return(NULL)
    before <- fit
    calculate <- function() {
        load_numerics()
        input <- if(fit$domain=='full_circle')
            list(pieces=circle_charts_cpp(b$raw_observations),counts=b$counts) else piecewise_input(b)
        workspace <- build_piecewise_model_cpp(input$pieces,input$counts,fit$eps,fit$tau,
            TRUE,TRUE,TRUE,fit$draws)
        on.exit(release_piecewise_model_cpp(workspace),add=TRUE)
        refine_rank_draws_cpp(workspace,fit$draws,fit$bounds$lower,fit$bounds$upper,
            as.integer(ids),target,as.integer(budget),controls$max_depth)
    }
    # One normal heavy slot, 1 GiB cap; existing bounded exclusive allocation
    # retry frees this workspace before changing caps. No new light admission.
    answer <- with_refinement_resources(calculate,memory_bytes=1024^3)
    stopifnot(all(answer$lower<=answer$upper),all(answer$lower>=fit$bounds$lower),
        all(answer$upper<=fit$bounds$upper),
        identical(answer$lower[,1L],fit$bounds$lower[,1L]),
        identical(answer$upper[,1L],fit$bounds$upper[,1L]))
    fit$bounds <- intersect_nested_bounds(pmax(fit$bounds$lower,answer$lower),
        pmin(fit$bounds$upper,answer$upper))
    fit$cutoff <- cutoff_order_statistic(fit$bounds$lower,fit$bounds$upper,fit$k,fit$policy$cutoff_tolerance)
    fit$cutoff$status <- ifelse(fit$cutoff$precision_met,'CUTOFF_GAP_ATTAINED','CUTOFF_BUDGET_UNRESOLVED')
    record <- list(engine='OBSERVED_ONLY_LOCAL_RANK_V1',draw_ids=ids,target=target,
        max_cells=as.integer(budget),max_depth=controls$max_depth,visited_cells=answer$visited_cells,
        remaining_peaks=answer$remaining_peaks,cutoff_before=before$cutoff,cutoff_after=fit$cutoff)
    fit$trace[[length(fit$trace)+1L]] <- record
    fit$source_md5 <- tools::md5sum(reproduction_files())
    stopifnot(identical(fit$draws,before$draws),identical(fit$backend,before$backend),
        identical(fit$eps,before$eps),identical(fit$tau,before$tau),
        all(fit$bounds$lower>=before$bounds$lower),all(fit$bounds$upper<=before$bounds$upper))
    list(fit=fit,record=record)
}

continue_coverage_decisions <- function(payload, task, name, controls, checkpoint) {
    if (is.null(controls) || payload$stage<controls$after_stage ||
        !task$family %in% c('circle','graph') ||
        payload$evaluation$coverage_lower[2L]==payload$evaluation$coverage_upper[2L]) return(payload)
    if (is.null(payload$decision_refinement))
        payload$decision_refinement <- list(controls=controls,passes=list())
    stopifnot(identical(payload$decision_refinement$controls,controls))
    completed <- length(payload$decision_refinement$passes)
    if (completed>=length(controls$budgets)) return(payload)
    for (pass in seq.int(completed+1L,length(controls$budgets))) {
        answer <- refine_local_decision_rank(payload$fit,controls,controls$budgets[pass])
        if (is.null(answer)) break
        payload$fit <- answer$fit
        payload$decision_refinement$passes[[pass]] <- answer$record
        # Truth only evaluates the already certified numerical interval.
        payload$evaluation <- evaluate_mean_set(payload$fit,task,name)
        payload$completion <- NULL
        checkpoint(payload)
        if (payload$evaluation$coverage_lower[2L]==payload$evaluation$coverage_upper[2L]) break
    }
    payload
}

compute_task <- function(task, input_file, signature, output, cfg = NULL) {
    started <- elapsed_time()
    input <- readRDS(input_file)
    stopifnot(unname(tools::md5sum(input_file)) == task$input_md5)
    if (is.null(cfg)) {
        m <- readRDS(file.path(output, "manifest.rds"))
        cfg <- m$signature$config
        stopifnot(identical(signature, m$signature_hash))
    }
    stopifnot(identical(cfg$objective, "coverage"))
    names_fit <- if (input$family == "ranking")
        cfg$metrics
    else input$family
    row_path <- file.path(output, "rows", paste0(task$id, ".csv"))
    previous_rows <- if (file.exists(row_path)) read.csv(row_path) else NULL
    completion <- list(objective = cfg$objective, policy = cfg$policy, decimal_gap = cfg$decimal_gap)
    if (!is.null(cfg$decision_refinement)) completion$decision_refinement <- cfg$decision_refinement
    certificates <- file.path(output, "certificates", paste0(task$id, "-", names_fit, ".rds"))
    comparison <- if (input$family %in% c("euclidean", "circle"))
        file.path(output, "comparisons", paste0(task$id, ".csv")) else character()
    if (!is.null(previous_rows) && all(file.exists(c(certificates, comparison)))) {
        cached <- lapply(certificates, function(path)
            verify_checkpoint(readRDS(path), signature, task$input_md5))
        ready <- vapply(cached, function(payload) {
            compatible <- !is.null(payload$completion) &&
                identical(payload$completion, completion)
            stopped <- !is.null(payload$evaluation) && all(payload$evaluation$coverage_lower == payload$evaluation$coverage_upper)
            compatible || (is.null(payload$completion) && stopped &&
                !is.null(payload$audit) && !is.null(payload$evaluation))
        }, logical(1))
        # A newly committed rank pass can precede publication of its CSV row.
        # Legacy completion flags alone must not reuse that stale row.
        coherent <- vapply(seq_along(cached),function(i) {
            z <- previous_rows[previous_rows$geometry==names_fit[i],,drop=FALSE]
            if(nrow(z)!=2L || !setequal(z$method,c('G','L')))return(FALSE)
            z <- z[match(c('G','L'),z$method),,drop=FALSE]
            f <- cached[[i]]$fit;ev <- cached[[i]]$evaluation
            !is.null(ev) && identical(as.numeric(z$cutoff_lower),unname(f$cutoff$lower)) &&
                identical(as.numeric(z$cutoff_upper),unname(f$cutoff$upper)) &&
                identical(as.numeric(z$coverage_lower),as.numeric(ev$coverage_lower)) &&
                identical(as.numeric(z$coverage_upper),as.numeric(ev$coverage_upper))
        },logical(1))
        if (all(ready) && all(coherent)) {
            stopifnot(nrow(previous_rows) == 2L * length(names_fit),
                      all(previous_rows$signature == signature), all(previous_rows$input_md5 == task$input_md5),
                      all(previous_rows$objective == cfg$objective),
                      identical(previous_rows$draw_hash,
                          rep(vapply(cached, function(payload) payload$fit$draw_hash, character(1)), each = 2L)))
            return(c(certificates, comparison, row_path))
        }
    }
    g <- regenerate_draws(input)
    rows <- list()
    paths <- character()
    for (name in names_fit) {
        path <- file.path(output, "certificates", paste0(task$id, "-", name, ".rds"))
        payload <- if (file.exists(path))
            verify_checkpoint(readRDS(path), signature, task$input_md5)
        else NULL
        if (is.null(payload))
            payload <- list(fit = NULL, stage = 0L, audit = NULL, evaluation = NULL)
        checkpoint <- function(value) save_atomic(list(signature=signature,input_md5=task$input_md5,
            payload=value,payload_hash=object_hash(value)),path)
        if (!is.null(payload$fit))
            payload <- continue_coverage_decisions(payload,task,name,cfg$decision_refinement,checkpoint)
        for (stage in seq_len(cfg$policy$stages)) {
            if (stage > 1L && !input$family %in% c("circle", "graph"))
                break
            if (stage <= payload$stage)
                next
            if (!is.null(payload$evaluation) &&
                all(payload$evaluation$coverage_lower == payload$evaluation$coverage_upper))
                break
            p <- cfg$policy
            p$stages <- 1L
            p$max_cells <- as.integer(p$max_cells * 2^(stage - 1L))
            p$cutoff_tolerance <- cfg$policy$cutoff_tolerance * (1 - .Machine$double.eps)
            requested <- c(TRUE, TRUE)
            if (!is.null(payload$fit)) {
                requested <- payload$evaluation$coverage_lower != payload$evaluation$coverage_upper
                gaps <- payload$fit$cutoff$gap[requested]
                if (all(gaps == 0)) break
                p$cutoff_tolerance <- min(p$cutoff_tolerance, min(gaps[gaps > 0]) / 2)
            }
            p$max_depth <- p$max_depth + stage - 1L
            p$direction_tolerance <- p$direction_tolerance/2^(stage - 1L)
            f <- if (stage == 1L)
                fit_observed_data(input, name, g, p, coarse = TRUE)
            else refine_cutoff(payload$fit, p, requested_methods = requested)
            payload$fit <- f
            payload$completion <- NULL
            payload$stage <- stage
            # Truth is evaluated only AFTER fitting. It can stop further work,
            # but cannot change any draw, feasible witness or global upper bound.
            payload$evaluation <- evaluate_mean_set(f, task, name)
            record <- list(signature = signature, input_md5 = task$input_md5, payload = payload, payload_hash = object_hash(payload))
            save_atomic(record, path)
            payload <- continue_coverage_decisions(payload,task,name,cfg$decision_refinement,checkpoint)
            f <- payload$fit
            if (!input$family %in% c("circle", "graph"))
                break
        }
        f <- payload$fit
        payload$audit <- verify_cutoff_gap(f, cfg$decimal_gap)
        if (is.null(payload$evaluation)) payload$evaluation <- evaluate_mean_set(f, task, name)
        stopifnot(all(payload$evaluation$coverage_lower <= payload$evaluation$coverage_upper))
        payload$completion <- completion
        save_atomic(list(signature = signature, input_md5 = task$input_md5, payload = payload, payload_hash = object_hash(payload)),
            path)
        attained <- vapply(payload$audit$log, function(s) grepl("attained=True", s, fixed = TRUE), logical(1))
        stopifnot(length(attained) == 2L)
        ev <- payload$evaluation
        rows[[name]] <- data.frame(id = task$id, family = task$family, parameter = task$parameter, n = task$n,
            replication = task$replication, geometry = name, method = c("G", "L"), B = f$B, k = f$k, cutoff_lower = f$cutoff$lower,
            cutoff_upper = f$cutoff$upper, cutoff_used = f$cutoff$used, gap = f$cutoff$gap, attained = attained,
            status = ifelse(ev$coverage_lower == ev$coverage_upper,
                "COVERAGE_DECISION_ATTAINED", "COVERAGE_BUDGET_UNRESOLVED"), coverage_lower = ev$coverage_lower,
            coverage_upper = ev$coverage_upper, evaluation_scope = ev$evaluation_scope, stage = payload$stage,
            max_draw_gap = apply(f$bounds$upper - f$bounds$lower, 2, max), input_md5 = task$input_md5, draw_hash = f$draw_hash,
            signature = signature)
        paths <- c(paths, path)
        rows[[name]]$objective <- cfg$objective
        rows[[name]]$coverage_resolved <- ev$coverage_lower == ev$coverage_upper
    }
    rows <- do.call(rbind, rows)
    comparisons <- evaluate_comparisons(input, task, payload$fit)
    if (!is.null(comparisons)) {
        comp_path <- file.path(output, "comparisons", paste0(task$id, ".csv"))
        save_atomic(comparisons, comp_path)
        paths <- c(paths, comp_path)
    }
    rows$elapsed_seconds <- elapsed_time() - started
    fields <- c("id", "method", "objective", "input_md5", "draw_hash")
    if (!is.null(previous_rows) && identical(unname(as.matrix(rows[fields])), unname(as.matrix(previous_rows[fields]))))
        rows$elapsed_seconds <- pmax(rows$elapsed_seconds, previous_rows$elapsed_seconds)
    save_atomic(rows, row_path)
    c(paths, row_path)
}

read_settings <- function(profile = "settings.dcf", reps = NULL, workers = NULL) {
    p <- as.list(read.dcf(profile)[1, ])
    num <- function(key) as.numeric(strsplit(p[[key]], " ", fixed = TRUE)[[1]])
    policy <- refinement_policy(as.numeric(p[["Decimal-cutoff-gap"]]), num("Stages"), num("Initial-cells"), num("Initial-depth"),
        num("Pair-batch"), num("Draw-batch"), num("Initial-direction-tolerance"))
    cfg <- list(schema = p$Schema, objective = p$Objective, reps = as.integer(if (is.null(reps)) num("Replications") else reps), workers = as.integer(if (is.null(workers)) num("Workers") else workers),
        execution_attempts = if (is.null(p[["Execution-attempts"]])) 1L else as.integer(num("Execution-attempts")),
        n = as.integer(num("Sample-sizes")), delta = num("Euclidean-delta"), lambda = num("Circle-lambda"), metrics = c("euclidean",
            "kendall", "footrule"), B = as.integer(num("B")), alpha = num("Alpha"), decimal_gap = p[["Decimal-cutoff-gap"]],
        policy = policy)
    stopifnot(cfg$schema == "certified-gl-targets-v1", cfg$objective == "coverage",
        cfg$reps >= 1L, cfg$reps <= 1000L, cfg$workers >= 1L, cfg$workers <=
        32L, cfg$execution_attempts >= 1L, cfg$execution_attempts <= 3L,
        identical(cfg$n, c(50L, 100L, 200L)), identical(cfg$delta, c(0, 0.5, 1)), identical(cfg$lambda, c(-0.125,
        0, 0.5, 1)), cfg$B == 1999L, cfg$alpha == 0.05, cfg$decimal_gap == ".01", policy$max_cells * 2^(policy$stages -
        1) == num("Maximum-cells"), policy$max_depth + policy$stages - 1 == num("Maximum-depth"))
    if (!is.null(p[["Circle-angles"]])) cfg$lengths <- circle_length_controls(
        enabled=identical(p[["Circle-angles"]],"Global Local EL"), tolerance=num("Angle-tolerance-degrees"),
        stages=num("Angle-stages"),initial_cells=num("Angle-initial-cells"),
        maximum_cells=num("Angle-maximum-cells"),calibration_cells=num("Angle-calibration-cells"))
    if (!is.null(p[["Circle-angles"]])) stopifnot(p[["Circle-angles"]] %in% c("Global Local EL","off"))
    if (!is.null(p[['Decision-refinement']])) {
        stopifnot(identical(p[['Decision-refinement']],'observed-only-local-rank-v1'))
        cfg$decision_refinement <- decision_refinement_controls(num('Decision-rank-cells'),
            num('Decision-rank-depth'),num('Decision-rank-draws'),num('Decision-rank-start-stage'),
            num('Decision-rank-fraction'))
    }
    if (!is.null(reps))
        stopifnot(reps == cfg$reps)
    if (!is.null(workers))
        stopifnot(workers == cfg$workers)
    cfg
}

regenerate_draws <- function(input) {
    dimension <- if (input$family %in% c("ranking", "graph"))
        24L
    else nrow(as.matrix(input$observations))
    g <- with_seed(input$validation_seed, matrix(rnorm(dimension * input$B), dimension, input$B))
    stopifnot(identical(object_hash(input$observations), input$observations_hash), identical(object_hash(g),
        input$original_draw_hash))
    g
}

evaluate_mean_set <- function(f, task, name) {
    if (task$family == "graph") {
        oracle <- ranking_graph_population_means(ranking_distribution("contamination", epsilon = c(1, 2)))
        query <- function(cuts) evaluate_graph_mean_set(f$backend, oracle, c(G = unname(cuts[1]), L = unname(cuts[2])))$decision
        lo <- query(f$cutoff$lower)
        hi <- query(f$cutoff$upper)
        return(data.frame(method = c("G", "L"), coverage_lower = as.integer(lo == "IN"), coverage_upper = as.integer(hi !=
            "OUT"), evaluation_scope = "EXACT_RATIONAL_WHOLE_POPULATION_MEAN_SET"))
    }
    points <- if (task$family == "euclidean")
        list(rep(0, ncol(f$Y)))
    else if (task$family == "circle") {
        load_numerics()
        intervals <- circle_population_mean_bounds_cpp(task$parameter)
        lapply(seq_len(nrow(intervals)), function(i) intervals[i, ])
    }
    else as.list(ranking_discrete_population_means(name)$mean_ids)
    rows <- mean_set_membership(f, points)$rows
    rows$evaluation_scope <- if (task$family == "circle" && task$parameter < 0)
        "WHOLE_POPULATION_MEAN_SET_WITH_CERTIFIED_QUERY_ENCLOSURES"
    else "WHOLE_POPULATION_MEAN_SET"
    rows
}

expected_counts <- function(cfg) {
    result <- list(tasks=27L*cfg$reps,unique_datasets=24L*cfg$reps,fits=33L*cfg$reps,
        cutoff_rows=66L*cfg$reps,summary_rows=66L,comparison_rows=42L*cfg$reps,comparison_summary_rows=42L)
    if(isTRUE(cfg$lengths$enabled))result<-c(result,list(length_datasets=12L*cfg$reps,
        length_rows=36L*cfg$reps,length_summary_rows=36L,length_paired_groups=36L))
    result
}

fit_observed_data <- function(input, name, g, policy, coarse = FALSE) {
    family <- input$family
    t <- if (family == "ranking")
        input$tuning[[name]]
    else input$tuning
    args <- list(draws = g, alpha = input$alpha, eps = t$eps, tau = t$tau, policy = policy, coarse = coarse)
    if (family == "euclidean")
        args$Y <- input$observations
    else if (family == "ranking") {
        return(calibrate_ranking(input$observations, name, draws = g, alpha = input$alpha, eps = t$eps, tau = t$tau,
            policy = policy))
    }
    else {
        args$backend <- if (family == "circle")
            circle_backend(input$observations)
        else ranking_graph_backend(ranking_graph(4L), input$observations)
        if (family == "graph")
            args$draws <- g[args$backend$seen, , drop = FALSE]
        args$witnesses <- input$witnesses[[name]]
    }
    do.call(calibrate_gl, args)
}

simulation_plan <- function(cfg) {
    paper <- graph <- "GENERATED"
    make <- function(family, parameters, root) {
        p <- expand.grid(replication = seq_len(cfg$reps), parameter = parameters, n = cfg$n, KEEP.OUT.ATTRS = FALSE)
        source_family <- if (family == "graph")
            "ranking"
        else family
        source_id <- sprintf("%s-p%s-n%d-r%06d", source_family, sprintf("%.17g", p$parameter), p$n, p$replication)
        p$family <- family
        p$id <- if (family == "graph")
            paste0("graph-", source_id)
        else source_id
        p$source <- file.path(root, "checkpoints", paste0(source_id, ".rds"))
        p
    }
    p <- rbind(make("euclidean", cfg$delta, paper), make("ranking", 0, paper), make("circle", cfg$lambda, paper),
        make("graph", 0, graph))
    p <- p[order(p$replication, p$n, p$family, p$parameter), ]
    rownames(p) <- NULL
    stopifnot(nrow(p) == 27L * cfg$reps, !anyDuplicated(p$id))
    p
}

prepare_run <- function(output, cfg = read_settings(), freeze_only = FALSE) {
    stopifnot(!file.exists(file.path(output, "manifest.rds")))
    plan <- simulation_plan(cfg)
    files <- reproduction_files()
    hashes <- tools::md5sum(files)
    preparation <- list(config = cfg, plan = plan, source_md5 = hashes, runtime = runtime_versions())
    recover_atomic_publications(output, "PREPARATION.rds", list(preparation))
    if (dir.exists(output) && !file.exists(file.path(output, "PREPARATION.rds")) &&
        any(!list.files(output, all.files = TRUE, no.. = TRUE) %in% "source"))
        stop("Output already contains unrelated or unrecorded files; choose a fresh run directory.")
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output, mustWork = TRUE)
    for (file in files[grepl("[.]R$", files)]) parse(file)
    freeze_source(files, hashes, file.path(output, "source"), output)
    preparation_path <- file.path(output, "PREPARATION.rds")
    if (file.exists(preparation_path)) stopifnot(identical(readRDS(preparation_path), preparation))
    else save_atomic(preparation, preparation_path)
    if (freeze_only) return(invisible(preparation))
    streams <- simulation_streams()
    load_numerics()
    chunks <- split(seq_len(nrow(plan)), ceiling(seq_len(nrow(plan))/100L))
    inputs <- parallel::mclapply(chunks, function(ix) do.call(rbind, lapply(ix, function(i) {
        task <- as.list(plan[i, ])
        generate_input(task, output, streams)
    })), mc.cores = cfg$workers, mc.set.seed = FALSE)
    stopifnot(!any(vapply(inputs, inherits, logical(1), "try-error")))
    inputs <- do.call(rbind, inputs)
    plan$input_md5 <- inputs$input_md5[match(plan$id, inputs$id)]
    stopifnot(!anyNA(plan$input_md5), identical(hashes, tools::md5sum(files)))
    signature <- list(config = cfg, plan = plan, inputs = inputs, source_md5 = hashes, runtime = runtime_versions(),
        expected = expected_counts(cfg), private_library = normalizePath(.libPaths()[1], mustWork = TRUE))
    m <- list(signature = signature, signature_hash = object_hash(signature), created = as.character(Sys.time()))
    write.csv(plan, file.path(output, "plan.csv"), row.names = FALSE)
    write.csv(data.frame(file = files, md5 = unname(hashes)), file.path(output, "source_manifest.csv"), row.names = FALSE)
    save_atomic(list(status = "PREPARED", expected = signature$expected), file.path(output, "STATUS.rds"))
    save_atomic(m, file.path(output, "manifest.rds"))
    invisible(m)
}

summarize_run <- function(paths, output) {
    m <- readRDS(file.path(output, "manifest.rds"))
    cfg <- m$signature$config
    stopifnot(identical(cfg$objective, "coverage"))
    expected <- m$signature$expected
    input_paths <- file.path(output,"inputs",paste0(m$signature$plan$id,".rds"))
    source_hash <- tools::md5sum(file.path(output,"source",names(m$signature$source_md5)))
    names(source_hash) <- names(m$signature$source_md5)
    stopifnot(all(unname(tools::md5sum(input_paths))==m$signature$plan$input_md5),
      identical(source_hash,m$signature$source_md5),identical(runtime_versions(),m$signature$runtime))
    files <- sort(paths[grepl("[.]csv$", paths) & basename(dirname(paths)) == "rows"])
    stopifnot(length(files) == expected$tasks)
    rows <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
    stopifnot(nrow(rows) == expected$cutoff_rows, !anyDuplicated(rows[c("id", "geometry", "method")]), all(rows$B ==
        1999L), all(rows$k == 1900L), all(rows$cutoff_used == rows$cutoff_upper), all(rows$cutoff_lower <= rows$cutoff_upper),
        all(rows$coverage_lower <= rows$coverage_upper), all(rows$objective == cfg$objective),
        identical(rows$coverage_resolved, rows$coverage_lower == rows$coverage_upper),
        all(rows$signature == m$signature_hash), all(table(rows$id) %in%
            c(2L, 6L)))
    groups <- split(rows, interaction(rows$family, rows$parameter, rows$n, rows$geometry, rows$method, drop = TRUE))
    summary <- do.call(rbind, lapply(groups, function(z) {
        stopifnot(nrow(z) == cfg$reps)
        cbind(z[1, c("family", "parameter", "n", "geometry", "method", "evaluation_scope")], replications = nrow(z),
            attained = sum(z$attained), unresolved = sum(!z$attained), coverage_resolved = sum(z$coverage_resolved),
            coverage_unresolved = sum(!z$coverage_resolved), coverage_lower = mean(z$coverage_lower),
            coverage_upper = mean(z$coverage_upper), maximum_cutoff_gap = max(z$gap))
    }))
    stopifnot(nrow(summary) == expected$summary_rows)
    raw <- file.path(output, "raw.csv")
    summ <- file.path(output, "summary.csv")
    save_atomic(rows, raw)
    save_atomic(summary, summ)
    comparisons <- paths[basename(dirname(paths)) == "comparisons"]
    comparison_rows <- do.call(rbind, lapply(comparisons, read.csv, stringsAsFactors = FALSE))
    stopifnot(nrow(comparison_rows) == expected$comparison_rows, !anyDuplicated(comparison_rows[c("id", "method")]),
        all(comparison_rows$decision %in% c("IN", "OUT", "UNRESOLVED")))
    comp_path <- file.path(output, "comparisons.csv")
    save_atomic(comparison_rows, comp_path)
    groups <- split(comparison_rows, interaction(comparison_rows$family, comparison_rows$parameter, comparison_rows$n,
        comparison_rows$method, drop = TRUE))
    comp_summary <- do.call(rbind, lapply(groups, function(z) {
        stopifnot(nrow(z) == cfg$reps)
        cbind(z[1, c("family", "parameter", "n", "method", "assurance", "evaluation_scope")], replications = nrow(z),
            resolved = sum(z$decision != "UNRESOLVED"), unresolved = sum(z$decision == "UNRESOLVED"), numerically_certified = sum(z$numerically_certified),
            numerical_coverage_lower = mean(z$decision == "IN"), numerical_coverage_upper = mean(z$decision !=
                "OUT"))
    }))
    stopifnot(nrow(comp_summary) == expected$comparison_summary_rows)
    comp_summ_path <- file.path(output, "comparison_summary.csv")
    save_atomic(comp_summary, comp_summ_path)
    report_paths <- c(raw, summ, comp_path, comp_summ_path)
    report_hashes <- tools::md5sum(report_paths)
    names(report_hashes) <- basename(report_paths)
    status <- list(status = if (all(rows$coverage_resolved))
            "COMPLETE_COVERAGE_DECISIONS_ATTAINED" else "COMPLETE_WITH_UNRESOLVED_COVERAGE_DECISIONS",
        expected = expected, attained = sum(rows$attained), unresolved = sum(!rows$attained),
        coverage_resolved = sum(rows$coverage_resolved), coverage_unresolved = sum(!rows$coverage_resolved),
        signature = m$signature_hash, report_md5 = report_hashes)
    save_atomic(status, file.path(output, "STATUS.rds"))
    c(raw, summ, comp_path, comp_summ_path, file.path(output, "STATUS.rds"))
}

numerics_compiler_contract <- function() {
    setup <- new.env(parent = baseenv())
    sys.source("setup.R", envir = setup)
    setup$compiler_contract(normalizePath("."))
}

runtime_versions <- function() list(R = R.version.string, platform = R.version$platform, software = extSoftVersion(),
    packages = utils::installed.packages()[, c("Package", "Version", "Built"), drop = FALSE], python = system2("python3",
        "--version", stdout = TRUE), compiler = system2("g++", "--version", stdout = TRUE)[1L],
    compiler_contract = numerics_compiler_contract(), threads = Sys.getenv(c("OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")))

piecewise_input <- function(backend) {
    if (backend$domain != "continuous_ranking_graph")
        stop("graph backend required")
    point <- function(x) cbind(x, x)
    pieces <- lapply(seq_along(backend$pieces), function(j) {
        x <- backend$pieces[[j]]
        list(left = rep(x$left, 2L), right = rep(x$right, 2L), b = point(x$b), c = point(x$c), left_id = backend$graph$edges[j,
            1L], right_id = backend$graph$edges[j, 2L])
    })
    list(pieces = pieces, counts = backend$counts)
}

piecewise_bounds <- function(input, args, B = ncol(args$g), batch = 1024L, rng_before = NULL, refinement = NULL,
                             save_state = TRUE) {
    stopifnot(is.logical(save_state), length(save_state) == 1L, !is.na(save_state))
    total <- length(input$pieces)^2
    if (total > .Machine$integer.max)
        stop("pair index limit exceeded")
    if (is.null(batch))
        batch <- 1024L
    check_sample_size(B)
    check_sample_size(batch)
    dimension <- length(input$counts)
    supplied <- !is.null(args$g)
    if (!supplied && is.null(rng_before))
        stop("RNG state required for chunk replay")
    controls <- args[intersect(c("eps", "tau", "max_cells", "max_depth", "direction_tolerance", "continuous_cones"),
        names(args))]
    requested_methods <- if (is.null(refinement) || is.null(refinement$methods))
        rep(TRUE, 2)
    else refinement$methods
    model <- args$model
    if (is.null(model)) {
        model <- do.call(build_piecewise_model_cpp, c(input, controls[c("eps", "tau")], list(circle_prefix = !identical(args$circle_prefix,
            FALSE), radial_feasibility = !identical(args$radial_feasibility, FALSE), sharp_rectangles = isTRUE(args$sharp_rectangles), draws = if (supplied) args$g else NULL)))
        on.exit(release_piecewise_model_cpp(model), add = TRUE)
    }
    block <- if (total <= args$max_cells)
        total
    else max(1L, as.integer(args$max_cells%/%4L))
    if (!is.null(args$pair_batch))
        block <- min(block, args$pair_batch)
    answer <- list(lower = matrix(0, B, 2L), upper = matrix(0, B, 2L))
    blocks <- 0L
    seeded <- isTRUE(args$global_incumbent)
    if (seeded && !is.null(refinement$global_lower)) {
        seed <- refinement$global_lower
        if (!is.matrix(seed) || !identical(dim(seed), c(as.integer(B), 2L)) || any(!is.finite(seed)) || any(seed <
            0))
            stop("invalid certified global lower seeds")
        answer$lower[] <- seed
        answer$upper[] <- seed
    }
    skip_inactive <- isTRUE(args$skip_inactive_chunks) && supplied && !is.null(refinement) && any(refinement$active)
    if (skip_inactive) {
        inactive <- which(!refinement$active)
        answer$lower[inactive, ] <- pmax(answer$lower[inactive, , drop = FALSE], refinement$bounds$lower[inactive,
            , drop = FALSE])
        answer$upper[inactive, ] <- pmax(answer$upper[inactive, , drop = FALSE], refinement$bounds$upper[inactive,
            , drop = FALSE])
    }
    counters <- c("active_methods", "representatives", "outer_directions", "visited_cells", "leaves", "budget_leaves",
        "possible_cells", "represented_cells", "exact_selection_points")
    for (start in seq.int(0, total - 1, by = block)) {
        targets <- if (!is.null(args$cutoff_rank))
            matrix(rep(apply(answer$lower, 2L, function(x) sort(x, partial = args$cutoff_rank)[args$cutoff_rank]) +
                args$cutoff_tolerance * 0.99, each = B), B, 2L) else NULL
        geometry <- do.call(build_piecewise_cover_cpp, c(list(model_ptr = model), input, controls, list(pair_start = as.integer(start),
            pair_count = as.integer(min(block, total - start)), requested_methods = as.integer(sum(c(1L, 2L)[requested_methods])))))
        if (!supplied)
            assign(".Random.seed", rng_before, envir = .GlobalEnv)
        tryCatch({
            first <- TRUE
            for (lo in seq.int(1L, B, by = batch)) {
                ix <- seq.int(lo, min(B, lo + batch - 1L))
                active <- if (is.null(refinement))
                  seq_along(ix)
                else which(refinement$active[ix])
                if (skip_inactive && !length(active))
                  next
                g <- if (supplied)
                  args$g[, ix, drop = FALSE]
                else matrix(rnorm(dimension * length(ix)), dimension)
                chunk_lower <- chunk_upper <- matrix(0, length(ix), 2L)
                if (!is.null(refinement)) {
                  chunk_lower <- refinement$bounds$lower[ix, , drop = FALSE]
                  chunk_upper <- refinement$bounds$upper[ix, , drop = FALSE]
                }
                prior <- if (seeded && length(active))
                  list(lower = answer$lower[ix[active], , drop = FALSE], upper = answer$upper[ix[active], , drop = FALSE])
                else NULL
                current <- evaluate_piecewise_draws_cpp(geometry, if (length(active))
                  g[, active, drop = FALSE]
                else matrix(0, dimension, 1L), prior = prior,
                  draw_ids = if (supplied && !is.null(attr(input$pieces, "circle_observations")) &&
                    !identical(args$circle_prefix, FALSE) && length(active)) as.integer(ix[active]) else NULL,
                  cutoff_targets = if (!is.null(targets) && length(active)) targets[ix[active], , drop = FALSE] else NULL)
                if (length(active)) {
                  chunk_lower[active, ] <- current$lower
                  chunk_upper[active, ] <- current$upper
                }
                if (!is.null(refinement) && any(!requested_methods)) {
                  chunk_lower[, !requested_methods] <- refinement$bounds$lower[ix, !requested_methods, drop = FALSE]
                  chunk_upper[, !requested_methods] <- refinement$bounds$upper[ix, !requested_methods, drop = FALSE]
                }
                merge_piecewise_bounds_cpp(answer$lower, answer$upper, chunk_lower, chunk_upper, as.integer(lo - 1L))
                if (first) {
                  for (key in counters) answer[[key]] <- if (is.null(answer[[key]]))
                    current[[key]]
                  else answer[[key]] + current[[key]]
                  answer$fmin <- current$fmin
                  first <- FALSE
                }
            }
        }, finally = release_piecewise_cover_cpp(geometry))
        blocks <- blocks + 1L
    }
    answer$pair_start <- 0L
    answer$pair_count <- total
    answer$total_pairs <- total
    answer$blocks <- blocks
    answer$cell_budget_per_block <- args$max_cells
    answer$draw_batch <- min(batch, B)
    if (save_state) answer$solver_state <- piecewise_state_cpp(model)
    if (!is.null(refinement) && any(!requested_methods))
        answer$active_methods[!requested_methods] <- refinement$bounds$active_methods[!requested_methods]
    answer
}

compile_numerics <- function(file, dependencies, label, compiler = numerics_compiler_contract()) {
    key <- object_hash(list(tools::md5sum(dependencies), R.version.string, R.version$platform,
        as.character(utils::packageVersion("Rcpp")), compiler))
    cache <- file.path(".cache", paste0("rcpp-gl-", label), key)
    dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
    lock <- acquire_lock(paste0(cache, ".lock"), 180)
    on.exit(lock$kill(), add = TRUE)
    Rcpp::sourceCpp(file, cacheDir = cache, env = .GlobalEnv, rebuild = FALSE)
}

acquire_lock <- function(path, seconds = 300, poll = .05, shared = FALSE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    # A shell builtin waits on the pipe. Parent death closes it and releases flock.
    lock <- processx::process$new("flock", c("--no-fork", if (shared) "-s", "-w", as.character(seconds), path,
        "sh", "-c", "printf 'LOCKED\\n'; read -r release"), stdin = "|", stdout = "|", stderr = "|")
    started <- elapsed_time()
    repeat {
        if (length(lock$read_output_lines())) return(lock)
        if (!lock$is_alive() || elapsed_time() - started > seconds + 1) {
            lock$kill()
            stop("Could not acquire numerical initialization lock: ", path)
        }
        Sys.sleep(poll)
    }
}

mean_set_membership <- function(fit, points) {
    if (!is.list(points) || !length(points))
        stop("nonempty list of evaluation points required")
    queries <- lapply(points, function(q) confidence_set_contains(fit, q))
    rows <- lapply(seq_len(2L), function(j) {
        z <- do.call(rbind, lapply(queries, function(x) x$rows[j, ]))
        lo <- min(z$coverage_lower)
        hi <- min(z$coverage_upper)
        data.frame(method = c("G", "L")[j], coverage_lower = lo, coverage_upper = hi, decision = if (lo == 1L)
            "IN"
        else if (hi == 0L)
            "OUT"
        else "UNRESOLVED")
    })
    list(rows = do.call(rbind, rows), queries = queries, scope = "ALL_SUPPLIED_EVALUATION_POINTS_NO_TRUTH_IN_CALIBRATION")
}

exact_cutoff_index <- function(B, alpha) {
    path <- tempfile()
    on.exit(unlink(path))
    status <- system2("python3", c("common/exact.py", "index", B, sprintf("%.17g", alpha)), stdout = path)
    if (status != 0L)
        stop("exact cutoff index failed")
    k <- scan(path, quiet = TRUE)
    stopifnot(length(k) == 1L, k >= 1L, k <= B + 1L, k == floor(k))
    k
}

intersect_nested_bounds <- function(lower, upper) {
    lower[, 1L] <- pmax(lower[, 1L], lower[, 2L])
    upper[, 2L] <- pmin(upper[, 2L], upper[, 1L])
    if (any(lower > upper))
        stop("nested G/L enclosure conflict")
    list(lower = lower, upper = upper)
}

refinement_policy <- function(cutoff_tolerance = 0.01, stages = 4L, max_cells = 256L, max_depth = 2L, pair_batch = 64L,
    draw_batch = 1024L, direction_tolerance = 0.1, rank_refinement = TRUE) {
    scalar <- function(x, lo, hi) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lo && x <= hi
    integer <- function(x, lo, hi) scalar(x, lo, hi) && x == floor(x)
    if (!scalar(cutoff_tolerance, 0, Inf) || !integer(stages, 1, 20) || !integer(max_cells, 1, .Machine$integer.max) ||
        !integer(max_depth, 0, 30) || !integer(pair_batch, 1, max_cells) || !integer(draw_batch, 1, .Machine$integer.max) ||
        !scalar(direction_tolerance, .Machine$double.xmin, 2) || max_cells * 2^(stages - 1) > .Machine$integer.max ||
        !is.logical(rank_refinement) || length(rank_refinement) != 1L || is.na(rank_refinement))
        stop("invalid fixed-B enclosure policy")
    list(cutoff_tolerance = cutoff_tolerance, stages = as.integer(stages), max_cells = as.integer(max_cells), max_depth = as.integer(max_depth),
        pair_batch = as.integer(pair_batch), draw_batch = as.integer(draw_batch), direction_tolerance = direction_tolerance,
        rank_refinement = rank_refinement)
}

acquire_refinement_resources <- function(exclusive = FALSE,
    directory = file.path("/tmp", paste0("frechet-refinement-", system2("id", "-u", stdout = TRUE))),
    slots = 4L, memory_bytes = 8 * 1024^3, headroom_bytes = 8 * 1024^3, root_global = FALSE,
    resource_kind = NULL) {
    stopifnot(.Platform$OS.type == "unix", nzchar(Sys.which("prlimit")), slots >= 1L,
        memory_bytes > 0, headroom_bytes >= 0, is.logical(root_global), length(root_global) == 1L,
        !is.na(root_global))
    prefix <- if (root_global) "root-slot-" else "slot-"
    other_prefix <- if (root_global) "slot-" else "root-slot-"
    other_slots <- if (root_global) 4L else 32L
    other_memory <- if (root_global) 8 * 1024^3 else 1024^3
    if (is.null(resource_kind)) resource_kind <- if (root_global) "GLOBAL_ROOT" else "REFINEMENT"
    stopifnot(is.character(resource_kind), length(resource_kind) == 1L, !is.na(resource_kind),
        if (root_global) resource_kind %in% c("GLOBAL_ROOT", "CIRCLE_RANK_LOCAL") else resource_kind == "REFINEMENT")
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    pid <- Sys.getpid()
    status <- file.path(directory, paste0("worker-", pid, ".json"))
    record <- function(phase, ...) save_atomic(list(pid = pid, phase = phase, exclusive = exclusive,
        resource_kind = resource_kind,
        updated = as.character(Sys.time()), ...), status)
    rss <- function(pid) {
        path <- file.path("/proc", pid, "status")
        if (!file.exists(path)) return(0)
        x <- tryCatch(grep("^VmRSS:", readLines(path, warn = FALSE), value = TRUE), error = function(e) character())
        if (length(x)) as.numeric(strsplit(trimws(x), "[[:space:]]+")[[1L]][2L]) * 1024 else 0
    }
    available <- function() {
        x <- grep("^MemAvailable:", readLines("/proc/meminfo"), value = TRUE)
        stopifnot(length(x) == 1L)
        as.numeric(strsplit(trimws(x), "[[:space:]]+")[[1L]][2L]) * 1024
    }
    try_lock <- function(path) {
        lock <- processx::process$new("flock", c("--no-fork", "-n", "-E", "75", path,
            "sh", "-c", "printf 'LOCKED\\n'; read -r release"), stdin = "|", stdout = "|", stderr = "|")
        repeat {
            if (length(lock$read_output_lines())) return(lock)
            if (!lock$is_alive()) {
                if (identical(lock$get_exit_status(), 75L)) return(NULL)
                stop("Resource slot lock failed: ", lock$read_error())
            }
            Sys.sleep(.005)
        }
    }
    if (!exists("address_space_limits_cpp", mode = "function")) load_numerics()
    original_limits <- address_space_limits_cpp()
    original <- ifelse(is.infinite(original_limits), "unlimited",
        format(original_limits, scientific = FALSE, trim = TRUE))
    cap <- if (original[1L] == "unlimited") Inf else as.numeric(original[1L])
    admission <- NULL; locks <- list(); paths <- character(); changed <- FALSE; published <- FALSE
    on.exit(if (!published) {
        if (changed) address_space_limits_cpp(original_limits)
        unlink(paths)
        for (lock in locks) lock$kill()
        if (!is.null(admission)) admission$kill()
        try(record("ERROR"), silent = TRUE)
    })
    record("WAITING")
    repeat {
        if (is.null(admission) || !admission$is_alive())
            admission <- acquire_lock(file.path(directory, "admission.lock"), 604800, 1)
        locks <- list(); outstanding <- 0; other_busy <- 0L; foreign_exclusive <- FALSE
        for (i in seq_len(slots)) {
            path <- file.path(directory, paste0(prefix, i))
            lock <- try_lock(paste0(path, ".lock"))
            if (!is.null(lock)) locks[[as.character(i)]] <- lock else {
                metadata <- tryCatch(jsonlite::read_json(paste0(path, ".json")), error = function(e) NULL)
                foreign_exclusive <- foreign_exclusive || isTRUE(metadata$exclusive)
                outstanding <- outstanding + if (is.null(metadata)) memory_bytes else
                    max(0, metadata$limit_bytes - rss(metadata$pid))
            }
        }
        # Root-only G and bounded Circle rank work share a small-memory queue.
        # Both queues use
        # this admission lock and subtract each other's still-unused reservations.
        # An exclusive retry waits for the other queue to drain; its published
        # exclusive lease prevents new admissions in that queue until release.
        for (i in seq_len(other_slots)) {
            path <- file.path(directory, paste0(other_prefix, i))
            probe <- try_lock(paste0(path, ".lock"))
            if (!is.null(probe)) probe$kill() else {
                other_busy <- other_busy + 1L
                metadata <- tryCatch(jsonlite::read_json(paste0(path, ".json")), error = function(e) NULL)
                foreign_exclusive <- foreign_exclusive || isTRUE(metadata$exclusive)
                outstanding <- outstanding + if (is.null(metadata)) other_memory else
                    max(0, metadata$limit_bytes - rss(metadata$pid))
            }
        }
        limit <- min(cap, if (exclusive) min(32 * 1024^3,
            available() + rss(pid) - headroom_bytes - outstanding) else memory_bytes)
        ready <- !foreign_exclusive && if (exclusive)
            length(locks) == slots && other_busy == 0L && limit >= min(cap, 2 * memory_bytes) else
            length(locks) > 0L && available() >= headroom_bytes + outstanding + max(0, limit - rss(pid))
        if (ready) {
            if (!exclusive && length(locks) > 1L) {
                for (lock in locks[-1L]) lock$kill()
                locks <- locks[1L]
            }
            virtual_line <- grep("^VmSize:", readLines(file.path("/proc",pid,"status")), value=TRUE)
            virtual_bytes <- as.numeric(strsplit(trimws(virtual_line), "[[:space:]]+")[[1L]][2L])*1024
            if (limit < virtual_bytes + 16*1024^2)
                stop("cannot allocate refinement workspace below current virtual address usage")
            address_space_limits_cpp(c(floor(limit),original_limits[2L]))
            changed <- TRUE
            paths <- file.path(directory, paste0(prefix, names(locks), ".json"))
            for (path in paths) save_atomic(list(pid = pid, limit_bytes = limit, exclusive = exclusive,
                resource_kind = resource_kind), path)
            admission$kill()
            record("ACTIVE", limit_bytes = limit, slots = length(locks))
            published <- TRUE
            return(list(locks = locks, paths = paths, original = original, original_limits = original_limits, status = status, pid = pid))
        }
        for (lock in locks) lock$kill()
        if (!exclusive) admission$kill()
        Sys.sleep(1)
    }
}

release_refinement_resources <- function(lease) {
    stopifnot(identical(lease$pid, Sys.getpid()))
    # Restore before any filesystem/process work, without spawning a child.
    address_space_limits_cpp(lease$original_limits)
    unlink(lease$paths)
    for (lock in lease$locks) lock$kill()
    save_atomic(list(pid = lease$pid, phase = "IDLE", updated = as.character(Sys.time())), lease$status)
    invisible(NULL)
}

with_refinement_resources <- function(compute, retry = function() NULL, ...) {
    attempt <- function(exclusive) {
        lease <- acquire_refinement_resources(exclusive = exclusive, ...)
        on.exit(release_refinement_resources(lease))
        compute()
    }
    tryCatch(attempt(FALSE), error = function(error) {
        if (!grepl("cannot allocate|bad_alloc|bad allocation|memory exhausted",
            conditionMessage(error), ignore.case = TRUE)) stop(error)
        gc()
        retry()
        message("Memory allocation limit reached; retrying once with exclusive refinement resources.")
        attempt(TRUE)
    })
}

refine_cutoff <- function(fit, policy = refinement_policy(cutoff_tolerance = 0.01, stages = 8L, max_cells = 512L, max_depth = 3L,
    direction_tolerance = 0.05, draw_batch = 1024L), checkpoint_directory = NULL,
    requested_methods = c(TRUE, TRUE)) {
    compute <- function() refine_cutoff_enclosed(fit, policy, checkpoint_directory, requested_methods)
    # Scheduling only: retain the same observations, draws, target and cell
    # budget. The light class is restricted to the measured n<=200, B=1999
    # Circle routes; a real allocation failure retries once with all queues idle.
    small_circle <- inherits(fit, "gl_certificate") && identical(fit$domain, "full_circle") &&
        isTRUE(fit$B == 1999L) && isTRUE(fit$backend$n <= 200L)
    light_global <- small_circle && identical(unname(requested_methods), c(TRUE, FALSE))
    light_local <- small_circle && identical(unname(requested_methods), c(FALSE, TRUE)) &&
        isTRUE(policy$max_cells >= (length(fit$backend$raw_observations) + 1L)^2) &&
        isTRUE(policy$max_depth >= 13L) && isTRUE(policy$max_cells * 2^(policy$stages - 1L) <= 262144L)
    if (light_global || light_local)
        return(with_refinement_resources(compute, root_global = TRUE, slots = 32L, memory_bytes = 1024^3,
            resource_kind = if (light_global) "GLOBAL_ROOT" else "CIRCLE_RANK_LOCAL",
            retry = function() {
                if (!is.null(checkpoint_directory) && dir.exists(checkpoint_directory))
                    stopifnot(file.rename(checkpoint_directory, paste0(checkpoint_directory, ".allocation-", Sys.getpid())))
            }))
    if (policy$max_cells * 2^(policy$stages - 1L) <= 8192L) return(compute())
    with_refinement_resources(compute, retry = function() {
        if (!is.null(checkpoint_directory) && dir.exists(checkpoint_directory))
            stopifnot(file.rename(checkpoint_directory, paste0(checkpoint_directory, ".allocation-", Sys.getpid())))
    })
}

refine_cutoff_enclosed <- function(fit, policy = refinement_policy(cutoff_tolerance = 0.01, stages = 8L, max_cells = 512L, max_depth = 3L,
    direction_tolerance = 0.05, draw_batch = 1024L), checkpoint_directory = NULL,
    requested_methods = c(TRUE, TRUE)) {
    if (!inherits(fit, "gl_certificate") || !fit$domain %in% c("full_circle", "continuous_ranking_graph") || fit$schema !=
        "fixed-input-certified-gl-v1" || fit$assurance != "DIRECTED_INTERVAL_FIXED_BINARY64_INPUT_AND_FIXED_B_ORDER_STATISTIC" ||
        !identical(policy, do.call(refinement_policy, policy)) || !identical(object_hash(fit$draws), fit$draw_hash) ||
        nrow(fit$bounds$lower) != fit$B || ncol(fit$draws) != fit$B || !is.logical(requested_methods) || length(requested_methods) != 2L || anyNA(requested_methods))
        stop("trusted fixed-input continuous certificate required")
    b <- fit$backend
    identity <- b[intersect(c("domain", "counts", "seen", "raw_observations", "observations", "pieces", "loss",
        "pairs", "metric"), names(b))]
    if (is.function(identity$loss))
        identity$loss <- NULL
    if (!identical(object_hash(identity), fit$input_hash))
        stop("continuous input hash mismatch")
    if (!is.null(checkpoint_directory)) {
        if (file.exists(checkpoint_directory) || dir.exists(checkpoint_directory))
            stop("fresh checkpoint directory required")
        dir.create(checkpoint_directory, recursive = TRUE)
    }
    started <- elapsed_time()
    parent <- list(bounds_hash = object_hash(fit$bounds), source_md5 = fit$source_md5, input_hash = fit$input_hash,
        draw_hash = fit$draw_hash, cutoff = fit$cutoff)
    load_numerics()
    if (fit$k != cutoff_index_cpp(fit$B, fit$alpha))
        stop("continuous rank mismatch")
    input <- if (fit$domain == "full_circle")
        list(pieces = circle_charts_cpp(b$raw_observations), counts = b$counts)
    else piecewise_input(b)
    workspace <- build_piecewise_model_cpp(input$pieces, input$counts, fit$eps, fit$tau,
        TRUE, TRUE, TRUE, fit$draws)
    on.exit(release_piecewise_model_cpp(workspace), add = TRUE)
    # The threshold solver budgets ALL root pairs per draw; the block solver
    # budgets each block separately. Keep the complete block cover below this gate.
    global_root_only <- fit$domain == "full_circle" && identical(unname(requested_methods), c(TRUE, FALSE))
    priority <- global_root_only || (fit$domain == "full_circle" && policy$max_cells >= length(input$pieces)^2 &&
        (policy$max_depth >= 13L || isTRUE(fit$cutoff$precision_met[2])))
    if (!priority && length(fit$metadata$solver_state)) piecewise_state_cpp(workspace, fit$metadata$solver_state)
    fit$policy <- policy
    fit$cutoff <- cutoff_order_statistic(fit$bounds$lower, fit$bounds$upper, fit$k, policy$cutoff_tolerance)
    fit$cutoff$status <- ifelse(fit$cutoff$precision_met, "CUTOFF_GAP_ATTAINED", "CUTOFF_BUDGET_UNRESOLVED")
    fit$polish <- list(parent = parent, policy = policy, passes = 0L, elapsed_seconds = 0)
    for (pass in seq_len(policy$stages)) {
        rank <- cutoff_order_statistic(fit$bounds$lower, fit$bounds$upper, fit$k, policy$cutoff_tolerance)
        if (all(rank$precision_met | !requested_methods))
            break
        methods <- !rank$precision_met & requested_methods
        active <- rep(FALSE, fit$B)
        for (j in which(methods[1:2])) active <- active | (fit$bounds$upper[, j] > rank$lower[j] & fit$bounds$lower[,
            j] <= rank$upper[j])
        if (!policy$rank_refinement)
            active[] <- TRUE
        lo <- hi <- matrix(0, fit$B, 2L)
        lo[, 1:2] <- fit$bounds$lower
        hi[, 1:2] <- fit$bounds$upper
        lo[active, which(methods)] <- hi[active, which(methods)] <- 0
        args <- list(g = fit$draws, eps = fit$eps, tau = fit$tau, max_cells = as.integer(policy$max_cells *
            2^(pass - 1L)), max_depth = min(30L, policy$max_depth + pass - 1L), direction_tolerance = policy$direction_tolerance/2^(pass -
            1L), continuous_cones = TRUE, circle_prefix = TRUE, radial_feasibility = TRUE, sharp_rectangles = TRUE,
            global_incumbent = TRUE, skip_inactive_chunks = TRUE, pair_batch = policy$pair_batch,
            cutoff_rank = if (policy$rank_refinement) fit$k else NULL, cutoff_tolerance = policy$cutoff_tolerance,
            model = workspace)
        if (priority) {
            cover <- list(lower = fit$bounds$lower, upper = fit$bounds$upper)
            if (methods[2]) {
                cover <- refine_rank_draws_cpp(workspace, fit$draws, fit$bounds$lower, fit$bounds$upper,
                    which(active), rank$lower[2] + policy$cutoff_tolerance*.99, args$max_cells)
            }
            if (methods[1]) {
                # G already uses the independent ROOT cone cover. Rebuilding
                # millions of Local leaves cannot sharpen that G upper bound.
                root_args <- args; root_args$max_depth <- 0L
                global <- piecewise_bounds(input, root_args, fit$B, policy$draw_batch,
                    refinement = list(methods = c(TRUE, FALSE), active = active,
                        global_lower = fit$bounds$lower,
                        bounds = list(lower = lo, upper = hi, active_methods = integer(2L))), save_state = FALSE)
                cover$lower[, 1] <- pmax(cover$lower[, 1], global$lower[, 1])
                cover$upper[, 1] <- pmin(cover$upper[, 1], global$upper[, 1])
            }
            cover$engine <- if (all(methods[1:2])) "ROOT_G_AND_RANK_LOCAL" else
                if (methods[1]) "ROOT_CONE_G" else "RANK_THRESHOLD_BRANCH_AND_BOUND"
            cover$draw_ids <- which(active)
            cover$target <- setNames(rank$lower[which(methods[1:2])] + policy$cutoff_tolerance*.99,
                c("G", "L")[which(methods[1:2])])
        } else cover <- piecewise_bounds(input, args, fit$B, policy$draw_batch, refinement = list(methods = methods,
            active = active, global_lower = fit$bounds$lower, bounds = list(lower = lo, upper = hi, active_methods = integer(2L))))
        fit$bounds <- intersect_nested_bounds(pmax(fit$bounds$lower, cover$lower[, 1:2]),
            pmin(fit$bounds$upper, cover$upper[, 1:2]))
        fit$cutoff <- cutoff_order_statistic(fit$bounds$lower, fit$bounds$upper, fit$k, policy$cutoff_tolerance)
        fit$cutoff$status <- ifelse(fit$cutoff$precision_met, "CUTOFF_GAP_ATTAINED", "CUTOFF_BUDGET_UNRESOLVED")
        fit$metadata$solver_state <- if (!all(fit$cutoff$precision_met)) cover$solver_state else NULL
        fit$trace[[length(fit$trace) + 1L]] <- list(pass = pass, active_draws = sum(active), cutoff_lower = fit$cutoff$lower,
            cutoff_upper = fit$cutoff$upper, gap = fit$cutoff$gap, max_cells = args$max_cells,
            max_depth = if (priority) { if (methods[2]) 60L else 0L } else args$max_depth,
            metadata = cover[c("engine", "draw_ids", "target", "visited_cells", "remaining_peaks", "budget_leaves", "blocks", "total_pairs", "active_methods")])
        fit$policy <- policy
        fit$source_md5 <- tools::md5sum(reproduction_files())
        fit$polish <- list(parent = parent, policy = policy, passes = pass, elapsed_seconds = elapsed_time() -
            started)
        if (!is.null(checkpoint_directory))
            saveRDS(fit, file.path(checkpoint_directory, paste0("pass", pass, ".rds")), version = 3)
    }
    fit$polish$elapsed_seconds <- elapsed_time() - started
    fit
}

confidence_set_contains <- function(fit, q, timeout = 60) {
    stopifnot(inherits(fit, "gl_certificate"),
        fit$domain %in% c("full_euclidean", "full_circle", "finite_ranking"))
    load_numerics()
    query <- function(cutoff) {
        if (fit$domain == "full_euclidean") {
            ans <- euclid_exact(fit$Y, "query", q, cutoff)$decision
            return(list(decision = ans, certified = TRUE))
        }
        b <- fit$backend
        if (fit$domain == "full_circle")
            return(circle_membership(b, q, cutoff, timeout))
        if (length(q) != 1L || !is.numeric(q) || !q %in% seq_len(ncol(b$loss)))
            stop("invalid ranking vertex query")
        ranking_finite_membership(b$loss[, q], b$loss[, q], b$loss, b$loss, b$counts, cutoff)
    }
    endpoints <- lapply(seq_len(2L), function(j) list(lower = query(fit$cutoff$lower[j]), upper = query(fit$cutoff$upper[j])))
    rows <- lapply(seq_len(2L), function(j) {
        low <- endpoints[[j]]$lower
        high <- endpoints[[j]]$upper
        decision <- if (isTRUE(low$certified) && low$decision == "IN")
            "IN"
        else if (isTRUE(high$certified) && high$decision == "OUT")
            "OUT"
        else "UNRESOLVED"
        if (low$decision == "IN" && high$decision == "OUT")
            stop("cutoff-monotonicity conflict")
        data.frame(method = c("G", "L")[j], lower_cutoff_decision = low$decision, upper_cutoff_decision = high$decision,
            ideal_finite_B_decision = decision, coverage_lower = as.integer(isTRUE(low$certified) && low$decision ==
                "IN"), coverage_upper = as.integer(!(isTRUE(high$certified) && high$decision == "OUT")))
    })
    list(rows = do.call(rbind, rows), endpoints = endpoints, assurance = "FULL_COMPARATOR_INVERSION_AND_CUTOFF_ENCLOSURE")
}

cutoff_order_statistic <- function(lower, upper, k, tolerance) {
    if (!is.matrix(lower) || !identical(dim(lower), dim(upper)) || ncol(lower) != 2L || anyNA(lower) || anyNA(upper) ||
        any(lower < 0) || any(lower > upper) || k < 1L || k > nrow(lower))
        stop("invalid G/L draw enclosure")
    lo <- apply(lower, 2L, function(x) sort(x, partial = k)[k])
    hi <- apply(upper, 2L, function(x) sort(x, partial = k)[k])
    names(lo) <- names(hi) <- c("G", "L")
    if (!exists("interval_width_cpp", mode = "function"))
        load_numerics()
    gap <- interval_width_cpp(lo, hi)
    names(gap) <- names(lo)
    gap[lo == hi] <- 0
    list(lower = lo, upper = hi, used = hi, gap = gap, precision_met = gap <= tolerance)
}

witness_coordinates <- function(backend, input, witnesses) {
    if (is.null(witnesses))
        return(matrix(numeric(), 0L, 4L))
    if (!is.matrix(witnesses) || any(!is.finite(witnesses)))
        stop("finite witness matrix required")
    if (backend$domain == "continuous_ranking_graph") {
        if (ncol(witnesses) != 4L)
            stop("graph witnesses: piece_q,t_q,piece_r,t_r")
        return(unique(witnesses))
    }
    if (ncol(witnesses) != 2L || any(abs(witnesses) > pi))
        stop("circle witnesses: q,r in [-pi,pi]")
    locate <- function(x) {
        ids <- which(vapply(input$pieces, function(p) p$left[2L] < x && x < p$right[1L], logical(1L)))
        if (length(ids) != 1L)
            return(NULL)
        p <- input$pieces[[ids]]
        left <- mean(p$left)
        right <- mean(p$right)
        u <- (x - left)/(right - left)
        if (!is.finite(u) || u < 0 || u > 1)
            return(NULL)
        c(ids, u)
    }
    result <- lapply(seq_len(nrow(witnesses)), function(i) {
        a <- locate(witnesses[i, 1L])
        b <- locate(witnesses[i, 2L])
        if (is.null(a) || is.null(b))
            NULL
        else c(a, b)
    })
    if (!any(lengths(result)))
        matrix(numeric(), 0L, 4L)
    else unique(do.call(rbind, result))
}

load_numerics <- function() {
    compiler <- numerics_compiler_contract()
    files <- c("common/numerics.cpp", "common/interval.hpp",
               "euclid/numerics.hpp", "ranking/numerics.hpp", "circle/numerics.hpp", "circle/length.hpp")
    hash <- object_hash(list(source = tools::md5sum(files), compiler = compiler))
    if (!exists("euclid_projection_cpp", mode = "function") ||
        !identical(attr(euclid_projection_cpp, "source_md5"), hash)) {
        flags <- Sys.getenv("PKG_CPPFLAGS", unset = NA_character_)
        on.exit(if (is.na(flags)) Sys.unsetenv("PKG_CPPFLAGS") else Sys.setenv(PKG_CPPFLAGS = flags))
        Sys.setenv(PKG_CPPFLAGS = paste(if (is.na(flags)) "" else flags,
                                      "-I", shQuote(normalizePath("."))))
        compile_numerics("common/numerics.cpp", files, "numerics", compiler)
        attr(euclid_projection_cpp, "source_md5") <<- hash
    }
    invisible(NULL)
}

object_hash <- function(object) {
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    saveRDS(object, path, version = 3L, compress = FALSE)
    unname(tools::md5sum(path))
}

membership_result <- function(decision, lower = NA_real_, upper = NA_real_, witness = NULL, certified = FALSE, detail = NULL) {
    decision <- match.arg(decision, c("IN", "OUT", "UNRESOLVED"))
    structure(list(decision = decision, lower = lower, upper = upper, witness = witness, certified = certified,
        assurance = if (certified) "CERTIFIED" else "FLOAT64_REFERENCE", detail = detail), class = "membership_result")
}

exact_selector_masks <- function(pieces, counts, points, eps, tau) {
    input <- tempfile()
    queries <- tempfile()
    output <- tempfile()
    errors <- tempfile()
    on.exit(unlink(c(input, queries, output, errors)), add = TRUE)
    fmt <- function(x) sprintf("%.17g", x)
    s <- length(counts)
    J <- length(pieces)
    rows <- data.frame(piece = rep(seq_len(J) - 1L, each = s), left = rep(fmt(vapply(pieces, function(a) a$left[1],
        numeric(1))), each = s), right = rep(fmt(vapply(pieces, function(a) a$right[1], numeric(1))), each = s),
        b = fmt(unlist(lapply(pieces, function(a) a$b[, 1]), use.names = FALSE)), c = fmt(unlist(lapply(pieces,
            function(a) a$c[, 1]), use.names = FALSE)), count = rep(fmt(counts), times = J))
    write.csv(rows, input, row.names = FALSE)
    colnames(points) <- c("j", "k", "u", "v")
    write.csv(as.data.frame(lapply(as.data.frame(points), fmt)), queries, row.names = FALSE)
    status <- system2("python3", c("common/exact.py", "selection", shQuote(input), shQuote(queries), fmt(eps),
        fmt(tau)), stdout = output, stderr = errors, timeout = 120)
    if (status != 0)
        stop("exact selection fallback failed: ", paste(readLines(errors), collapse = " "))
    result <- read.csv(output)$mask
    if (length(result) != nrow(points) || anyNA(result) || any(!result %in% c(0, 1, 3)))
        stop("invalid exact masks")
    as.integer(result)
}

run_exact_solver <- function(command, args, output, errors, timeout) {
    for (attempt in 1:3) {
        status <- system2(command, args, stdout = output, stderr = errors, timeout = timeout)
        if (status != 0L)
            return(status)
        if (file.exists(output) && is.finite(file.info(output)$size) && file.info(output)$size > 0)
            return(status)
        if (attempt < 3L) {
            options(exact_solver_retries = getOption("exact_solver_retries", 0L) + 1L)
            Sys.sleep(0.05 * attempt)
        }
    }
    cat("empty solver output after three identical-input attempts\n", file = errors, append = TRUE)
    1L
}

with_seed <- function(seed, expr) {
    previous_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed)
        previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(previous_kind))
        if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv,
            inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("L'Ecuyer-CMRG", "Inversion")
    assign(".Random.seed", seed, envir = .GlobalEnv)
    eval.parent(substitute(expr))
}

write_exact_input <- function(rows, path) {
    rows[] <- lapply(rows, function(x) if (is.numeric(x))
        sprintf("%.17g", x)
    else x)
    write.csv(rows, path, row.names = FALSE)
}

simulation_streams <- function(master_seed = 202610074L) {
    reference <- simulation_plan(read_settings())
    ids <- sort(unique(sub("^graph-", "", reference$id)), method = "radix")
    previous_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(previous_kind))
        if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv) else
            rm(list = intersect(".Random.seed", ls(.GlobalEnv, all.names = TRUE)), envir = .GlobalEnv)
    })
    RNGkind("L'Ecuyer-CMRG", "Inversion")
    set.seed(master_seed)
    seed <- .Random.seed
    streams <- setNames(vector("list", length(ids)), ids)
    for (id in ids) {
        validation <- parallel::nextRNGSubStream(parallel::nextRNGSubStream(seed))
        streams[[id]] <- list(data = seed, validation = validation)
        seed <- parallel::nextRNGStream(seed)
    }
    streams
}

observed_tuning <- function(minimum, n) {
    stopifnot(length(minimum) == 1L, is.finite(minimum), minimum >= 0, n >= 3L)
    # v=0 is retained only in the frozen input schema to preserve reference hashes.
    list(eps = minimum * sqrt(log(log(n))/n), tau = n^(1/5), v = 0,
         rule = "OBSERVED_FSTAR_LOGLOG_SQRTN_TAU_FIFTH_G_L_ONLY")
}

generate_input <- function(task, output, streams) {
    seed <- streams[[sub("^graph-", "", task$id)]]
    observations <- with_seed(seed$data, switch(task$family,
        euclidean = sample_euclid(task$n, task$parameter),
        circle = sample_circle(task$n, task$parameter),
        ranking_distribution()$sample(task$n)))
    tuning <- switch(task$family,
        euclidean = {
            centered <- sweep(observations, 2, colMeans(observations))
            observed_tuning(mean(rowSums(centered^2)), task$n)
        },
        circle = observed_tuning(circle_empirical_mean(observations, tolerance = 0)$risk, task$n),
        ranking = setNames(lapply(c("euclidean", "kendall", "footrule"), function(metric) {
            risk <- drop(crossprod(tabulate(observations, 24L), ranking_discrete_geometry(metric)$loss))/task$n
            observed_tuning(min(risk), task$n)
        }), c("euclidean", "kendall", "footrule")),
        graph = {
            backend <- ranking_graph_backend(ranking_graph(), observations)
            minima <- vapply(backend$pieces, function(piece) {
                b <- sum(backend$p * piece$b)
                c <- sum(backend$p * piece$c)
                t <- min(1, max(0, -b/2))
                t^2 + b*t + c
            }, numeric(1))
            observed_tuning(min(minima), task$n)
        })
    dimension <- if (task$family %in% c("ranking", "graph")) 24L else task$n
    draws <- with_seed(seed$validation, matrix(rnorm(dimension * 1999L), dimension, 1999L))
    input <- list(family = task$family, observations = observations,
                  validation_seed = seed$validation, data_seed = seed$data,
                  observations_hash = object_hash(observations), original_draw_hash = object_hash(draws),
                  tuning = tuning, witnesses = NULL, alpha = 0.05, B = 1999L)
    path <- file.path(output, "inputs", paste0(task$id, ".rds"))
    save_frozen_input(input, path)
    data.frame(id = task$id, source_md5 = "GENERATED_FROM_RESERVED_RNG_STREAM",
               input_md5 = unname(tools::md5sum(path)))
}

save_frozen_input <- function(input, path) {
    if (file.exists(path)) stopifnot(identical(object_hash(readRDS(path)), object_hash(input)))
    else save_atomic(input, path)
    invisible(path)
}

for (geometry in c("euclid", "ranking", "circle")) {
    source(file.path(geometry, "geometry.R"))
}
source("common/reproduce.R")
source("ranking/report.R")

source("circle/length.R")
