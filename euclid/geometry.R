# Euclid: observed geometry, inference and simulation law.

euclid_bounds <- function(Y, draws, witnesses = NULL) {
    Y <- as.matrix(Y)
    if (!is.numeric(Y) || nrow(Y) < 2L || ncol(Y) < 1L || any(!is.finite(Y)) || nrow(draws) != nrow(Y))
        stop("invalid Euclidean observations/draw rows")
    if (!is.null(witnesses))
        stop("Euclidean calibration uses the analytic direction space")
    raw <- NULL
    support_certificate <- "EXACT_RATIONAL_ELIMINATION"
    if (ncol(Y) <= 3L && ncol(Y) < nrow(Y)) {
        trial <- euclid_projection_cpp(Y, draws, seq_len(ncol(Y)))
        if (isTRUE(trial$inverse_verified)) {
            raw <- trial
            columns <- seq_len(ncol(Y))
            support_certificate <- "VERIFIED_POSITIVE_GRAM_PIVOTS"
        }
    }
    if (is.null(raw)) {
        columns <- euclid_exact(Y)$column
        raw <- euclid_projection_cpp(Y, draws, as.integer(columns))
    }
    initial_inverse_verified <- raw$inverse_verified
    basis_enclosure <- NULL
    if (!initial_inverse_verified && length(columns)) {
        exact_basis <- euclid_exact(Y, "basis")
        r <- length(columns)
        basis_enclosure <- list(lower = as.matrix(exact_basis[, 2L * seq_len(r) - 1L, drop = FALSE]), upper = as.matrix(exact_basis[,
            2L * seq_len(r), drop = FALSE]))
        raw <- euclid_projection_cpp(Y, draws, as.integer(columns), basis_enclosure$lower, basis_enclosure$upper)
    }
    bounds <- list(lower = cbind(raw$lower, raw$lower), upper = cbind(raw$upper, raw$upper))
    metadata <- list(exact_rank = length(columns), support_columns = columns, support_certificate = support_certificate,
        initial_inverse_verified = initial_inverse_verified, exact_basis_preconditioned = !is.null(basis_enclosure),
        basis_enclosure = basis_enclosure, inverse_verified = raw$inverse_verified, fallback = if (raw$inverse_verified) "NONE" else "CENTERED_GAUSSIAN_NORM_UPPER")
    domain <- "full_euclidean"
    identity <- list(Y = Y)
    trace <- list()
    list(Y = Y, bounds = bounds, metadata = metadata, domain = domain, identity = identity, trace = trace)
}

euclid_exact <- function(Y, operation = "support", q = NULL, cutoff = NULL) {
    rows <- as.data.frame(Y)
    names(rows) <- paste0("y", seq_len(ncol(Y)) - 1L)
    if (operation == "query") {
        if (!is.numeric(q) || length(q) != ncol(Y) || any(!is.finite(q)))
            stop("invalid Euclidean query")
        for (j in seq_len(ncol(Y))) rows[[paste0("q", j - 1L)]] <- q[j]
        rows$cutoff <- cutoff
    }
    input <- tempfile()
    output <- tempfile()
    errors <- tempfile()
    on.exit(unlink(c(input, output, errors)), add = TRUE)
    write_exact_input(rows, input)
    status <- run_exact_solver("python3", c("common/exact.py", "euclidean", shQuote(input), operation), output,
        errors, 60)
    if (status != 0L)
        stop("exact Euclidean algebra failed: ", paste(readLines(errors, warn = FALSE), collapse = " "))
    read.csv(output, stringsAsFactors = FALSE)
}

sample_euclid <- function(n, delta, Z = NULL) {
    check_sample_size(n)
    if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta < 0 || delta > 1)
        stop("delta must lie in [0,1]")
    if (is.null(Z))
        Z <- matrix(rnorm(3L * n), n, 3L)
    Z <- as.matrix(Z)
    if (!is.numeric(Z) || !identical(dim(Z), c(as.integer(n), 3L)) || any(!is.finite(Z)))
        stop("invalid Euclidean latent sample")
    # Use intrinsic orthonormal coordinates for the singular model. Forming
    # a binary64 ambient projection introduces a spurious residual direction.
    # Draw all three latent columns first to preserve the reserved RNG stream.
    if (delta == 0) return(Z[, 1:2, drop = FALSE])
    u <- rep(1/sqrt(3), 3L)
    projection <- diag(3L) - tcrossprod(u)
    Z %*% projection + sqrt(delta) * tcrossprod(drop(Z %*% u), u)
}
