#!/usr/bin/env Rscript
# Install exact source archives into a project library; no private path needed.
toolchain_output <- function(command, arguments) {
    value <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(value, "status")) || !length(value))
        stop("Cannot inspect recorded toolchain: ", command)
    value
}

compiler_contract <- function(project = ".", runtime = NULL) {
    if (is.null(runtime))
        runtime <- as.list(read.dcf(file.path(project, "environment", "runtime.dcf"))[1, ])
    environment <- Sys.getenv()
    overrides <- grepl(paste0("^(PKG_|R_MAKEVARS|R_XTRA_|USE_CXX|SHLIB_|",
        "CC$|CXX|CPP|CFLAGS$|CPICFLAGS$|FC$|F77$|FFLAGS$|FCFLAGS$|LDFLAGS$|",
        "MAKEFLAGS$|MFLAGS$|MAKEOVERRIDES$|MAKEFILES$|GNUMAKEFLAGS$|",
        "CPATH$|CPLUS_INCLUDE_PATH$|C_INCLUDE_PATH$|LIBRARY_PATH$|COMPILER_PATH$|",
        "GCC_EXEC_PREFIX$|R_PKG_CXX_STD$|R_CROSS_BUILD$)"), names(environment))
    if (any(overrides & nzchar(environment)))
        stop("Unsupported compiler environment override: ",
             paste(names(environment)[overrides & nzchar(environment)], collapse = ", "))
    if (!identical(Sys.getenv("MAKE", "make"), "make") || nzchar(Sys.getenv("R_ARCH")) ||
        !identical(Sys.getenv("R_PLATFORM", R.version$platform), R.version$platform))
        stop("The recorded compiler requires MAKE=make, an empty R_ARCH and the recorded R_PLATFORM.")
    makevars <- c(file.path(R.home("etc"), "Makevars.site"),
        path.expand(paste0("~/.R/Makevars-", R.version$platform)), path.expand("~/.R/Makevars"))
    if (any(file.exists(makevars)))
        stop("Custom Makevars are unsupported; use the recorded R compiler configuration.")
    fields <- c("CXX17", "CXX17FLAGS", "CXX17STD", "CXX17PICFLAGS", "CPPFLAGS",
                "LDFLAGS", "SHLIB_CXX17LD", "SHLIB_CXX17LDFLAGS")
    config <- stats::setNames(lapply(fields, function(name)
        toolchain_output(file.path(R.home("bin"), "R"), c("CMD", "config", name))), fields)
    for (name in fields)
        if (!identical(config[[name]], runtime[[name]]))
            stop("Recorded R compiler configuration differs: ", name)
    compiler <- paste(runtime$CXX17, toolchain_output(runtime$CXX17, "-dumpfullversion"))
    if (!identical(compiler, runtime$Compiler))
        stop("Recorded compiler is ", runtime$Compiler, "; found ", compiler)
    rules <- c(Makeconf = file.path(R.home("etc"), "Makeconf"),
               shlib.mk = file.path(R.home("share"), "make", "shlib.mk"))
    stopifnot(all(file.exists(rules)))
    rule_hashes <- stats::setNames(unname(tools::md5sum(rules)), names(rules))
    expected <- c(Makeconf = runtime$`Makeconf-MD5`, shlib.mk = runtime$`Shlib-MD5`)
    if (!identical(rule_hashes, expected))
        stop("Recorded R compiler build rules differ; see environment/runtime.dcf.")
    list(config = config, compiler = compiler, build_rules = rule_hashes)
}

check_toolchain_contract <- function(runtime) {
    python <- toolchain_output("python3", "--version")
    if (!identical(python, paste("Python", runtime$Python)))
        stop("Recorded Python is ", runtime$Python, "; found ", paste(python, collapse = " "))
    optimize <- toolchain_output("python3", c("-c", shQuote("import sys; print(sys.flags.optimize)")))
    if (!identical(optimize, "0"))
        stop("Python optimization disables scientific assertions; unset PYTHONOPTIMIZE and do not use -O.")
    compiler_contract(runtime = runtime)
    blas <- basename(unname(extSoftVersion()["BLAS"]))
    if (!identical(blas, runtime$`BLAS-library`) || !identical(La_version(), runtime$LAPACK))
        stop("Recorded BLAS/LAPACK differs; see environment/runtime.dcf.")
    invisible(TRUE)
}

check_runtime_contract <- function(root, library, check_packages = TRUE) {
    runtime <- as.list(read.dcf(file.path(root, "environment", "runtime.dcf"))[1, ])
    if (.Platform$OS.type != "unix" || !grepl("linux", R.version$platform))
        stop("Use Linux or WSL; Windows-native R is not the recorded environment.")
    if (as.character(getRversion()) != runtime$R || R.version$platform != runtime$Platform)
        stop("Recorded R/platform is ", runtime$R, " / ", runtime$Platform,
             "; found ", getRversion(), " / ", R.version$platform)
    commands <- strsplit(runtime$`Required-commands`, " ", fixed = TRUE)[[1]]
    missing <- commands[!nzchar(Sys.which(commands))]
    if (length(missing)) stop("Missing system commands: ", paste(missing, collapse = ", "))
    check_toolchain_contract(runtime)
    if (!check_packages) return(invisible(TRUE))
    if (!dir.exists(library)) stop("Project library is missing; run setup.R install first: ", library)
    library <- normalizePath(library, mustWork = TRUE)
    lock <- read.dcf(file.path(root, "environment", "packages.dcf"))
    stopifnot(!anyDuplicated(lock[, "Package"]), all(grepl("^[a-f0-9]{64}$", lock[, "SHA256"])))
    for (i in seq_len(nrow(lock))) {
        p <- lock[i, "Package"]; version <- lock[i, "Version"]
        allowed <- if ("Scope" %in% colnames(lock) && isTRUE(lock[i, "Scope"] == "R-recommended"))
            c(library, .Library) else library
        path <- find.package(p, lib.loc = allowed, quiet = TRUE)
        if (length(path) != 1L || read.dcf(file.path(path, "DESCRIPTION"))[1, "Version"] != version)
            stop("Locked package mismatch: ", p, " ", version, "; run setup.R install.")
        active <- find.package(p, quiet = TRUE)
        if (length(active) != 1L || !identical(normalizePath(active), normalizePath(path)))
            stop("A different library shadows locked package: ", p)
        if (isNamespaceLoaded(p) && as.character(getNamespaceVersion(p)) != version)
            stop("A different package version is already loaded: ", p, "; restart R with --vanilla.")
    }
    invisible(TRUE)
}

main <- function(args = commandArgs(TRUE)) {
    if (length(args) > 2L || (length(args) && !args[1] %in% c("install", "check")))
        stop("Usage: Rscript --vanilla setup.R install|check [LIBRARY]")
    script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
    if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
    root <- dirname(normalizePath(script, mustWork = TRUE))
    source(file.path(root, "common", "reproduce.R"), local = .GlobalEnv)
    library <- if (length(args) == 2L) args[2] else file.path(root, ".R-library")
    mode <- if (length(args)) args[1] else "check"
    runtime <- as.list(read.dcf(file.path(root, "environment", "runtime.dcf"))[1, ])
    check_runtime_contract(root, library, check_packages = FALSE)
    lock <- read.dcf(file.path(root, "environment", "packages.dcf"))
    stopifnot(!anyDuplicated(lock[, "Package"]), all(grepl("^[a-f0-9]{64}$", lock[, "SHA256"])))
    if (mode == "install") dir.create(library, recursive = TRUE, showWarnings = FALSE)
    library <- normalizePath(library, mustWork = TRUE)
    .libPaths(c(library, .Library))
    cache <- file.path(root, ".cache", "package-sources")
    if (mode == "install") dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    inventory_libraries <- inventories <- list()
    package_inventory <- function(allowed) {
        index <- which(vapply(inventory_libraries, identical, logical(1), allowed))
        if (!length(index)) {
            index <- length(inventories) + 1L
            inventory_libraries[[index]] <<- allowed
            inventories[[index]] <<- installed.packages(lib.loc = allowed, noCache = TRUE)
        }
        inventories[[index]]
    }
    for (i in seq_len(nrow(lock))) {
        p <- lock[i, "Package"]; version <- lock[i, "Version"]
        allowed <- if ("Scope" %in% colnames(lock) && isTRUE(lock[i,"Scope"]=="R-recommended"))
            c(library,.Library) else library
        installed <- package_inventory(allowed)
        if (p %in% rownames(installed) && installed[p, "Version"] == version) next
        if (mode == "check") stop("Missing locked package: ", p, " ", version)
        archive <- file.path(cache, paste0(p, "_", version, ".tar.gz"))
        if (!file.exists(archive)) {
            temporary <- tempfile(tmpdir = cache)
            tryCatch({
                urls <- unique(c(lock[i,"Source"],sprintf(
                    "https://cloud.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz",p,p,version)))
                downloaded <- FALSE
                for (url in urls) {
                    status <- tryCatch(suppressWarnings(download.file(url,temporary,mode="wb",quiet=TRUE)),
                                       error=function(e) 1L)
                    if (identical(status,0L)) {downloaded <- TRUE;break}
                }
                if (!downloaded) stop("Cannot fetch locked source: ",p," ",version)
                stopifnot(sha256_file(temporary) == lock[i, "SHA256"])
                stopifnot(file.rename(temporary, archive))
            }, finally = unlink(temporary))
        }
        if (sha256_file(archive) != lock[i, "SHA256"]) stop("Source checksum mismatch: ", archive)
        cat("Installing", p, version, "from its verified source archive\n")
        install.packages(archive, repos = NULL, type = "source", lib = library,
                         INSTALL_opts = c("--no-multiarch", "--no-test-load"))
        # A project installation can change both the project-only inventory and
        # the first matching version in the project-plus-system inventory.
        for (j in seq_along(inventories))
            if (library %in% inventory_libraries[[j]])
                inventories[[j]] <- installed.packages(lib.loc = inventory_libraries[[j]], noCache = TRUE)
        installed <- package_inventory(library)
        stopifnot(p %in% rownames(installed), installed[p, "Version"] == version)
    }
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    check_runtime_contract(root, library)
    for (p in strsplit(runtime$`Direct-packages`, " ", fixed = TRUE)[[1]])
        stopifnot(requireNamespace(p, quietly = TRUE))
    cat("PASS:", nrow(lock), "locked packages; R", runtime$R, "; library", library, "\n")
    cat("Recorded compiler/Python/BLAS are in environment/runtime.dcf.\n")
}

if (sys.nframe() == 0L) main()
