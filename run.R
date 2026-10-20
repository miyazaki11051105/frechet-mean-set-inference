#!/usr/bin/env Rscript
# One command prepares, runs and exports; the same command resumes frozen work.
configure_runtime <- function() {
    library_path <- Sys.getenv("FRECHET_PRIVATE_LIBRARY",
        unset = if (dir.exists(".R-library")) ".R-library" else ".private-r-library")
    if (dir.exists(library_path)) {
        .libPaths(c(normalizePath(library_path, mustWork = TRUE), .Library.site, .Library))
        Sys.setenv(FRECHET_PRIVATE_LIBRARY = normalizePath(library_path, mustWork = TRUE))
    }
    Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep),
               OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    setup <- new.env(parent = baseenv())
    sys.source("setup.R", envir = setup)
    setup$check_runtime_contract(normalizePath("."), library_path)
}

write_figures <- function(project, output, options = character()) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(file.path(project, "plot.R")), shQuote(output), shQuote(options)))
    if (status != 0L) stop("Figure export failed; numerical results are preserved.")
    invisible(NULL)
}

verify_results_inputs <- function(output) {
    manifest <- readRDS(file.path(output,"manifest.rds"))
    coverage <- audit_coverage_run(output)
    require_certified_completion(coverage,manifest$signature$expected)
    directory <- paste0(output,"-ranking-candidates")
    stopifnot(verify_ranking_candidate_report(directory),
        identical(readRDS(file.path(directory,"manifest.rds"))$identity$run_signature,manifest$signature_hash))
    ranking <- jsonlite::read_json(ranking_report_paths(directory)[["audit"]])
    require_certified_completion(coverage,manifest$signature$expected,ranking)
    if (!isTRUE(manifest$signature$config$lengths$enabled))
        stop("Reference results require the enabled Circle angle analysis.")
    angles <- verify_circle_lengths(output)
    require_certified_completion(coverage,manifest$signature$expected,ranking,angles)
    invisible(TRUE)
}

verify_frozen_run <- function(output) {
    manifest <- readRDS(file.path(output, "manifest.rds"))
    signature <- manifest$signature
    stopifnot(identical(manifest$signature_hash, object_hash(signature)),
              identical(signature$runtime, runtime_versions()),
              identical(signature$expected, expected_counts(signature$config)))
    stopifnot(identical(unname(tools::md5sum(file.path(output, "source", names(signature$source_md5)))),
                        unname(signature$source_md5)))
    load_numerics()
    cat("Frozen source, runtime and numerical engines verified.\n")
}

verify_reports <- function(output, manifest, status) {
    stopifnot(identical(status$signature, manifest$signature_hash),
              identical(status$expected, manifest$signature$expected), length(status$report_md5) == 4L,
              identical(unname(tools::md5sum(file.path(output, names(status$report_md5)))), unname(status$report_md5)))
    invisible(status)
}

publish_execution <- function(output, phase, ...) {
    save_atomic(c(list(phase = phase, updated = as.character(Sys.time()), pid = Sys.getpid()), list(...)),
                file.path(output, "execution.json"))
}

recover_checkpoints <- function(output, manifest) {
    cfg <- manifest$signature$config
    plan <- manifest$signature$plan
    missing <- which(!file.exists(file.path(output, "rows", paste0(plan$id, ".csv"))))
    # Recover saved tasks before traversing the full targets graph; new tasks stay in the pool.
    for (i in head(missing, cfg$workers)) {
        task <- as.list(plan[i, ])
        fits <- if (task$family == "ranking") cfg$metrics else task$family
        if (!all(file.exists(file.path(output, "certificates", paste0(task$id, "-", fits, ".rds"))))) next
        publish_execution(output, "RECOVERING_CHECKPOINTS", task = task$id, workers = cfg$workers)
        paths <- compute_task(task, file.path(output, "inputs", paste0(task$id, ".rds")),
                              manifest$signature_hash, output, cfg)
        save_atomic(list(id = task$id, files_md5 = tools::md5sum(paths), updated = as.character(Sys.time())),
                    file.path(output, "execution-recovery", paste0(task$id, ".rds")))
    }
    invisible(NULL)
}

run_pipeline <- function(output, project = file.path(output, "source")) {
    manifest <- readRDS(file.path(output, "manifest.rds"))
    status_path <- file.path(output, "STATUS.rds")
    verify_frozen_run(output)
    if (file.exists(status_path)) {
        status <- readRDS(status_path)
        if (grepl("^COMPLETE_", status$status)) return(verify_reports(output, manifest, status))
    }
    .libPaths(c(manifest$signature$private_library, .libPaths()))
    Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep),
               FRECHET_CERTIFIED_RUN = output)
    setwd(project)
    warm_status <- system2(file.path(R.home("bin"), "Rscript"),
                           c("--vanilla", "run.R", "warm", shQuote(output)))
    stopifnot(warm_status == 0L)

    # The entry point holds the Linux advisory lock for this controller.
    begin_execution_epoch(output, manifest$signature_hash)
    save_atomic(list(status = "RUNNING", started = as.character(Sys.time()),
                 pid = Sys.getpid(), workers = manifest$signature$config$workers,
                 signature = manifest$signature_hash), status_path)
    recover_checkpoints(output, manifest)
    # Resolve at runtime so targets sees a definition file without a nested run.
    make <- getExportedValue("targets", "tar_make")
    attempts <- manifest$signature$config$execution_attempts
    if (is.null(attempts)) attempts <- 1L
    errors <- list()
    for (attempt in seq_len(attempts)) {
        publish_execution(output, "SOLVING", attempt = attempt, workers = manifest$signature$config$workers)
        error <- tryCatch({
            make(script = "_targets.R", store = file.path(output, "_targets"),
             reporter = "timestamp", seconds_meta_append = 5,
             envir = new.env(parent = baseenv()))
            NULL
        }, error = conditionMessage)
        metadata <- targets::tar_meta(fields = c("name", "error"),
                                  store = file.path(output, "_targets"))
        failed <- !is.na(metadata$error) & nzchar(metadata$error)
        if (is.null(error) && !any(failed)) return(invisible(readRDS(status_path)))
        errors[[attempt]] <- list(attempt = attempt, time = as.character(Sys.time()),
                                  error = error, targets = metadata[failed, ])
        save_atomic(errors, file.path(output, "execution-errors.rds"))
        message("Execution attempt ", attempt, " failed; completed checkpoints are retained.")
    }
    save_atomic(list(status = "INCOMPLETE_EXECUTION_ERRORS", errors = errors,
                     signature = manifest$signature_hash), status_path)
    stop("Execution attempts exhausted; rerun the same command to resume. See execution-errors.rds.")
}

figure_publication_ready <- function(output, manifest, status) {
    directory <- file.path(output, "figures")
    link <- Sys.readlink(directory)
    if (!is.na(link) && nzchar(link)) stop("The figure directory must not be a symbolic link.")
    path <- file.path(output, "figures", "figure-manifest.rds")
    figure <- if (file.exists(path)) suppressWarnings(tryCatch(readRDS(path), error = function(e) NULL)) else NULL
    if (!is.list(figure)) return(FALSE)
    if (!is.null(figure$run_signature) && !identical(figure$run_signature, manifest$signature_hash))
        stop("Existing figure manifest belongs to a different run; preserve it and choose a fresh figure directory.")
    if (!is.null(figure$report_md5) && !identical(figure$report_md5, status$report_md5))
        stop("Existing figure manifest refers to different numerical reports; no figure files were changed.")
    expected <- file.path(output, "figures", c("figure-coverage.tex", "figure-data.csv"))
    links <- Sys.readlink(expected)
    is.character(figure$files_md5) && length(figure$files_md5) == length(expected) &&
        !anyNA(figure$files_md5) && all(grepl("^[a-f0-9]{32}$", figure$files_md5)) &&
        !anyDuplicated(names(figure$files_md5)) && setequal(names(figure$files_md5), expected) &&
        all(file.exists(expected)) && !any(dir.exists(expected)) && !any(!is.na(links) & nzchar(links)) &&
        identical(figure$run_signature, manifest$signature_hash) && identical(figure$report_md5, status$report_md5) &&
        identical(unname(tools::md5sum(names(figure$files_md5))), unname(figure$files_md5))
}

ensure_run_figures <- function(output, project, manifest, status) {
    if (figure_publication_ready(output, manifest, status)) return(invisible(NULL))
    directory <- file.path(output, "figures")
    if (file.exists(directory)) stopifnot(file.rename(directory, tempfile("figures-incomplete-", tmpdir = output)))
    write_figures(project, output, if (manifest$signature$config$reps == 1000L) character() else "--pilot")
    stopifnot(figure_publication_ready(output, manifest, status))
    invisible(NULL)
}

finish_run <- function(output, project) {
    verify_frozen_run(output)
    manifest <- readRDS(file.path(output, "manifest.rds"))
    status <- readRDS(file.path(output, "STATUS.rds"))
    verify_reports(output, manifest, status)
    coverage_audit <- audit_coverage_run(output)
    require_certified_completion(coverage_audit,manifest$signature$expected)
    publish_execution(output, "EXPORTING", workers = manifest$signature$config$workers)
    ensure_run_figures(output, project, manifest, status)
    candidate_directory <- paste0(output, "-ranking-candidates")
    publish_execution(output, "EXPORTING", workers = manifest$signature$config$workers,
        ranking_candidates = list(status = "ANALYZING", directory = candidate_directory))
    candidate_paths <- write_ranking_candidate_report(output, candidate_directory,
        workers = min(4L, manifest$signature$config$workers))
    candidate_audit <- jsonlite::read_json(candidate_paths[["audit"]])
    require_certified_completion(coverage_audit,manifest$signature$expected,ranking=candidate_audit)
    candidate_record <- list(status = "COMPLETE", signature = candidate_audit$signature,
        paths = as.list(candidate_paths), report_md5 = candidate_audit$report_md5,
        candidate_method_rows = candidate_audit$candidate_method_rows,
        unresolved_candidate_methods = candidate_audit$unresolved_candidate_methods)
    length_record <- NULL
    if(isTRUE(manifest$signature$config$lengths$enabled)) {
        publish_execution(output,"ANGLES",workers=manifest$signature$config$workers,ranking_candidates=candidate_record)
        length_status <- run_circle_lengths(output,project)
        require_certified_completion(coverage_audit,manifest$signature$expected,candidate_audit,length_status)
        length_record <- list(status=length_status$status,rows=length_status$rows,
            precision_attained=length_status$attained,valid_wide=length_status$wide,
            reports=as.list(length_status$reports),signature=length_status$signature)
    }
    audit_coverage_run(output)
    if (manifest$signature$config$reps == 1000L) {
        export_status <- system2("python3", c(shQuote(file.path(project,"export_results.py")), "export", "--run",
            shQuote(output), "--output", shQuote(file.path(output,"results"))))
        if (export_status != 0L) stop("Results export failed; certified numerical results are preserved.")
    }
    publish_execution(output, "COMPLETE",
                      workers = manifest$signature$config$workers,
                      coverage_resolved = status$coverage_resolved, coverage_unresolved = status$coverage_unresolved,
                      circle_angles = length_record, ranking_candidates = candidate_record)
    invisible(status)
}


main <- function(arguments = commandArgs(TRUE)) {
    locked <- "--controller" %in% arguments
    arguments <- arguments[arguments != "--controller"]
    commands <- c("prepare", "run", "warm", "plot", "audit", "export-results", "ranking-counts", "length-worker")
    if (length(arguments) < 2L || !arguments[1L] %in% commands) {
        stop("Usage: Rscript run.R run OUTPUT [REPS WORKERS], or prepare|plot OUTPUT [options]")
    }
    script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
    # Rscript encodes spaces in --file when invoked by a frozen-run subprocess.
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
    script <- normalizePath(script, mustWork = TRUE)
    project <- dirname(script)
    output <- arguments[2L]
    if (!grepl("^(/|[A-Za-z]:)", output)) output <- file.path(getwd(), output)
    path_arguments <- if (arguments[1L] %in% c("ranking-counts","export-results","plot")) 3L else integer()
    for (j in path_arguments[path_arguments<=length(arguments)])
        if (!startsWith(arguments[j],"--") && !grepl("^(/|[A-Za-z]:)",arguments[j]))
            arguments[j] <- file.path(getwd(),arguments[j])
    setwd(project)
    configure_runtime()
    command <- arguments[1L]
    manifest_path <- file.path(output,"manifest.rds")
    if (command %in% c("prepare","run") && file.exists(manifest_path)) {
        historical <- readRDS(manifest_path)
        if (!is.null(historical$parent) || isTRUE(historical$signature$config$size_reps > 0L))
            stop("Historical length-enabled runs are preserved read-only. Use audit or a fresh Coverage output.")
    }
    if (command == "length-worker") {
        if (file.exists(file.path(output,"repair","manifest.rds")))
            stop("Historical repaired runs are read-only.")
        source("common/inference.R")
        return(circle_length_worker(output))
    }
    if (command == "audit") {
        source("common/inference.R")
        result<-audit_coverage_run(output)
        if(isTRUE(readRDS(file.path(output,"manifest.rds"))$signature$config$lengths$enabled)) {
            lengths<-verify_circle_lengths(output)
            result$circle_angles<-lengths[c("status","rows","attained","wide","signature")]
        }
        manifest <- readRDS(file.path(output,"manifest.rds"))
        if ("export_results.py" %in% names(manifest$signature$source_md5) && manifest$signature$config$reps == 1000L) {
            verify_results_inputs(output)
            exports <- reproduction_results_report(output)
            result$reference_results <- list(status="PASS",tex_sha256=as.list(exports$tex_sha256))
        }
        return(print(result))
    }
    if (command == "ranking-counts") {
        if (!length(arguments) %in% 3:4)
            stop("Usage: Rscript run.R ranking-counts RUN NEW_DIRECTORY [WORKERS_1_TO_4]")
        source("common/inference.R")
        return(invisible(write_ranking_candidate_report(output,arguments[3L],
            workers=if (length(arguments)==4L) as.integer(arguments[4L]) else 1L)))
    }
    if (command == "export-results") {
        if (length(arguments) != 3L) stop("Usage: Rscript run.R export-results SAVED_RUN NEW_OUTPUT_DIRECTORY")
        source("common/inference.R")
        verify_results_inputs(output)
        result <- system2("python3", c("export_results.py", "export", "--run", shQuote(output),
            "--output", shQuote(arguments[3L])))
        if (result != 0L) stop("Results export failed; saved numerical results are preserved.")
        return(invisible(NULL))
    }
    if (command %in% c("prepare", "run") && !locked) {
        if (.Platform$OS.type != "unix" || !nzchar(Sys.which("flock"))) stop("Run under Linux/WSL with flock installed.")
        dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
        output <- if (dir.exists(output)) normalizePath(output, mustWork = TRUE) else
            file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
        result <- system2("flock", c("-n", "-E", "75", shQuote(paste0(output, ".lock")),
            shQuote(file.path(R.home("bin"), "Rscript")), "--vanilla", shQuote(script),
            shQuote(c(command, output, arguments[-c(1L, 2L)], "--controller"))))
        if (result == 75L) stop("This experiment already has an active controller or worker; no duplicate was started.")
        if (result != 0L) stop("Run incomplete; existing inputs and checkpoints are preserved. See the preceding error.")
        return(invisible(NULL))
    }
    if (command %in% c("prepare", "run", "warm") &&
        file.exists(file.path(output, "repair", "manifest.rds")))
        stop("Use the historical run's preserved source, or choose a fresh run directory for this version.")
    frozen <- file.path(output, "source")
    if (command %in% c("run", "prepare") &&
        (file.exists(file.path(output, "manifest.rds")) || file.exists(file.path(output, "PREPARATION.rds"))) &&
        normalizePath(project) != normalizePath(frozen)) {
        result <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(file.path(frozen, "run.R")),
            shQuote(c(command, output, arguments[-c(1L, 2L)], "--controller"))))
        if (result != 0L) stop("Frozen execution failed; inputs and checkpoints are preserved.")
        return(invisible(NULL))
    }
    if (command == "warm") setwd(frozen)
    source("common/inference.R")
    if (command == "run") {
        stopifnot(length(arguments) <= 4L)
        if (!file.exists(file.path(output, "manifest.rds"))) {
            if (file.exists(file.path(output, "PREPARATION.rds"))) {
                preparation <- readRDS(file.path(output, "PREPARATION.rds"))
                if (length(arguments) >= 3L) stopifnot(as.integer(arguments[3L]) == preparation$config$reps)
                if (length(arguments) >= 4L) stopifnot(as.integer(arguments[4L]) == preparation$config$workers)
                prepare_run(output, cfg = preparation$config)
            } else prepare_run(output, cfg = read_settings(
                reps = if (length(arguments) >= 3L) as.integer(arguments[3L]) else NULL,
                workers = if (length(arguments) >= 4L) as.integer(arguments[4L]) else NULL), freeze_only = TRUE)
            return(main(c("run", output, "--controller")))
        }
        output <- normalizePath(output, mustWork = TRUE)
        manifest <- readRDS(file.path(output, "manifest.rds"))
        if (length(arguments) >= 3L) stopifnot(as.integer(arguments[3L]) == manifest$signature$config$reps)
        if (length(arguments) >= 4L) stopifnot(as.integer(arguments[4L]) == manifest$signature$config$workers)
        tryCatch({ run_pipeline(output, project); finish_run(output, project) }, error = function(error) {
            state <- tryCatch(jsonlite::read_json(file.path(output, "execution.json")),
                error=function(e)list(phase='UNKNOWN'))
            if (inherits(error,'certification_incomplete'))
                publish_execution(output,'CERTIFICATION_INCOMPLETE',error=conditionMessage(error))
            else if (state$phase != "BUDGET_EXHAUSTED") publish_execution(output, "EXECUTION_ERROR", error = conditionMessage(error))
            stop(error)
        })
        return(invisible(NULL))
    }
    if (command == "plot") {
        options <- arguments[-c(1L, 2L)]
        return(write_figures(project, output, options))
    }
    if (command == "prepare") {
        if (length(arguments) != 2L) stop("Usage: Rscript run.R prepare OUTPUT")
        return(prepare_run(output))
    }
    if (command == "warm") return(verify_frozen_run(output))

}

if (sys.nframe() == 0L) {
    main()
    quit(save = "no", status = 0L)
}
