library(targets)

run <- normalizePath(Sys.getenv("FRECHET_CERTIFIED_RUN"), mustWork = TRUE)
record <- readRDS(file.path(run, "manifest.rds"))
workers <- record$signature$config$workers

compute_in_worker <- function(task, input_file, settings, output) {
    if (!exists(".numerics_ready", envir = .GlobalEnv, inherits = FALSE)) {
        # Rcpp updates its cache index on load. Only initialization is serialized.
        source("common/inference.R")
        lock <- acquire_lock(file.path(".cache", "worker-load.lock"))
        tryCatch({
            load_numerics()
            assign(".numerics_ready", TRUE, envir = .GlobalEnv)
        }, finally = lock$kill())
    }
    task_attempt(task, output, settings$signature)
    paths <- get("compute_task", envir = .GlobalEnv)(
        task, input_file, settings$signature, output, settings$config
    )
    paths <- ready_files(paths)
    task_attempt(task, output, settings$signature, complete = TRUE)
    paths
}

write_report <- function(paths, output) {
    source("common/inference.R")
    summarize_run(paths, output)
}

tar_option_set(
    packages = character(), library = .libPaths(), format = "rds", error = "continue",
    memory = "transient", garbage_collection = TRUE, seed = 20261011L,
    controller = crew::crew_controller_local(
        workers = workers, seconds_idle = 600, seconds_launch = 300,
        seconds_timeout = 300, reset_globals = FALSE, crashes_max = 3L,
        options_local = crew::crew_options_local(log_directory = file.path(run, "worker-logs"))
    )
)

coverage_targets <- list(
    tar_target(manifest_file, file.path(run, "manifest.rds"), format = "file", deployment = "main"),
    tar_target(manifest, readRDS(manifest_file), deployment = "main"),
    tar_target(settings, list(signature = manifest$signature_hash,
                             config = manifest$signature$config), deployment = "main"),
    tar_target(task, lapply(seq_len(nrow(manifest$signature$plan)), function(i)
        as.list(manifest$signature$plan[i, ])), iteration = "list", deployment = "main"),
    tar_target(input_files, file.path(run, "inputs", paste0(manifest$signature$plan$id, ".rds")),
               format = "file", deployment = "main"),
    tar_target(inputs_ready, {
        stopifnot(all(unname(tools::md5sum(input_files)) == manifest$signature$plan$input_md5))
        TRUE
    }, deployment = "main"),
    tar_target(certified, {
        stopifnot(inputs_ready)
        compute_in_worker(task, file.path(run, "inputs", paste0(task$id, ".rds")), settings, run)
    }, pattern = map(task), iteration = "list", format = "file"),
    tar_target(coverage_ready, {
        stopifnot(length(certified) == manifest$signature$expected$tasks,
                  all(base::lengths(certified) > 0L), all(file.exists(unlist(certified))),
                  all(file.exists(file.path(run, "rows", paste0(manifest$signature$plan$id, ".csv")))))
        source("common/inference.R")
        save_atomic(list(phase = "SUMMARIZING",
                         updated = as.character(Sys.time()), pid = Sys.getpid(),
                         workers = workers), file.path(run, "execution.json"))
        TRUE
    }, deployment = "main", error = "trim")
)
c(coverage_targets, list(
    tar_target(report, {
        stopifnot(isTRUE(coverage_ready))
        write_report(unlist(certified, use.names = FALSE), run)
    }, format = "file", deployment = "main")
))
