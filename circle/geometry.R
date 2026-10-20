# Circle: observed geometry, inference and simulation law.

circle_backend <- function (y, counts = rep(1, length(y)))
{
    check_counts(counts)
    if (length(y) != length(counts) || any(!is.finite(y)))
        stop("invalid circular sample")
    list(raw_observations = y, counts = counts, n = sum(counts), domain = "full_circle")
}

circle_empirical_mean <- function(y, tolerance = 1e-10) {
    if (!is.numeric(y) || !length(y) || any(!is.finite(y)) || !is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance < 0)
        stop("invalid circle observations")
    y <- (y + pi)%%(2 * pi) - pi
    n <- length(y)
    lifts <- ifelse(y > 0, y - 2 * pi, y)
    idx <- which(y != 0)
    cuts <- ifelse(y[idx] > 0, y[idx] - pi, y[idx] + pi)
    ord <- order(cuts, method = "radix")
    cuts <- cuts[ord]
    idx <- idx[ord]
    if (length(cuts)) {
        groups <- rle(cuts)
        ends <- cumsum(groups$lengths)
        levels <- groups$values
        crossed <- c(0L, ends)
        crossed_lift <- c(0, cumsum(lifts[idx])[ends])
    }
    else {
        levels <- numeric()
        crossed <- 0L
        crossed_lift <- 0
    }
    s1 <- sum(lifts) + 2 * pi * crossed
    s2 <- sum(lifts * lifts) + 4 * pi * crossed_lift + 4 * pi * pi * crossed
    left <- c(-pi, levels)
    right <- c(levels, pi)
    candidates <- pmin(right, pmax(left, s1/n))
    risks <- candidates * candidates - 2 * candidates * s1/n + s2/n
    best <- min(risks)
    winners <- candidates[risks <= best + tolerance]
    winners <- unique((winners + pi)%%(2 * pi) - pi)
    winner <- winners[which.min(abs(winners))]
    list(point = winner, points = sort(winners), risk = best)
}

circle_witnesses <- function(backend, eps) {
    # Proposals use observations only; the interval solver verifies feasibility.
    grid <- -pi + (seq_len(16L) - .5) * (2 * pi / 16)
    mean <- circle_empirical_mean(backend$raw_observations, tolerance = 0)$point
    q <- unique(c(grid, wrap_angle(mean + c(-1, -.5, 0, .5, 1) * sqrt(eps))))
    as.matrix(expand.grid(q = q, r = grid))
}

circle_nonnegative_coefficients <- function(lambda) {
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0 || lambda > 1)
        stop("invalid fixed-risk lambda: require a finite scalar in [0,1]")
    c(a = (4 - lambda)/6, b = 2 * (1 - lambda)/3)
}

circle_kde_curvature <- function(y, base, bandwidth = length(y)^(-1/5)) {
    x <- circle_score_data(y, base)
    if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
        stop("invalid periodic bandwidth")
    kappa <- bandwidth^(-2)
    if (!is.finite(kappa))
        stop("periodic kernel concentration overflow")
    normalizer <- 2 * pi * base::besselI(kappa, 0, expon.scaled = TRUE)
    density <- mean(exp(kappa * (-cos(x$z) - 1)))/normalizer
    list(curvature = 2 - 4 * pi * density, antipodal_density = density, bandwidth = bandwidth, kappa = kappa, kernel = "PERIODIC_VON_MISES")
}

circle_wald <- function(y, q = 0, alpha = 0.05) {
    x <- circle_score_data(y, q)
    n <- length(y)
    method <- "WALD_S1_KDE_CENTERED"
    result <- function(statistic = NA_real_, status, ...) circle_method_result(method, statistic, status, alpha, n = n,
        ...)
    minimum <- circle_empirical_mean(x$y, tolerance = 0)
    if (length(minimum$points) != 1L || !is.finite(minimum$point))
        return(result(status = "NONUNIQUE_EMPIRICAL_MEAN"))
    base <- minimum$point
    z <- wrap_angle(x$y - base)
    delta <- wrap_angle(minimum$point - x$q)
    if (any(abs(z) == pi) || abs(delta) == pi)
        return(result(status = "WALD_CUT_LOCUS"))
    kde <- circle_kde_curvature(x$y, base)
    score_variance <- mean((-2 * z - mean(-2 * z))^2)
    details <- c(list(mean = minimum$point, evaluation_point = base, displacement = delta, score_variance = score_variance,
        regularity = "POSITIVE_CURVATURE_ASSUMPTION_NOT_CERTIFIED"), kde)
    if (!is.finite(score_variance) || score_variance <= 0)
        return(do.call(result, c(list(status = "WALD_ZERO_SCORE_VARIANCE"), details)))
    if (!is.finite(kde$curvature) || kde$curvature <= 0)
        return(do.call(result, c(list(status = "WALD_NONPOSITIVE_ESTIMATED_CURVATURE"), details)))
    statistic <- n * delta^2 * kde$curvature^2/score_variance
    do.call(result, c(list(statistic = statistic, status = if (is.finite(statistic)) "WALD_COMPUTED" else "WALD_NUMERICAL_UNRESOLVED"),
        details))
}

circle_distribution_coefficients <- function(lambda) {
    bound <- (12 - 8 * sqrt(13))/43
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= bound || lambda > 1)
        stop("require lambda in (lambda_min,1]")
    c(a = (4 - lambda)/6, b = 2 * (1 - lambda)/3)
}

circle_population_means <- function(lambda) {
    circle_distribution_coefficients(lambda)
    m <- if (lambda < 0)
        acos(1 + 3 * lambda/(4 * (1 - lambda)))
    else 0
    means <- if (lambda < 0)
        c(-m, m)
    else 0
    list(means = means)
}

sample_circle <- function(n, lambda) {
    ab <- circle_distribution_coefficients(lambda)
    if (length(n) != 1L || !is.finite(n) || n < 3L || n != floor(n))
        stop("invalid n")
    if (lambda >= 0)
        return(sample_circle_nonnegative(n, lambda))
    envelope <- 1 + sum(ab)
    answer <- numeric(n)
    filled <- 0L
    while (filled < n) {
        size <- max(64L, as.integer(ceiling(1.1 * envelope * (n - filled))))
        x <- runif(size, -pi, pi)
        u <- runif(size)
        numerator <- 1 + ab[["a"]] * cos(x) + ab[["b"]] * cos(2 * x)
        accepted <- x[u * envelope <= numerator]
        take <- min(length(accepted), n - filled)
        if (take) {
            answer[filled + seq_len(take)] <- accepted[seq_len(take)]
            filled <- filled + take
        }
    }
    answer
}

circle_score_data <- function(y, q = 0) {
    if (!is.numeric(y) || length(y) < 3L || any(!is.finite(y)) || !is.numeric(q) || length(q) != 1L || !is.finite(q))
        stop("finite Circle observations (n >= 3) and a scalar q are required")
    yy <- wrap_angle(y)
    qq <- wrap_angle(q)
    list(y = yy, q = qq, z = wrap_angle(yy - qq))
}

circle_method_result <- function(method, statistic = NA_real_, status, alpha = 0.05, ...) {
    if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
        stop("invalid alpha")
    cutoff <- stats::qchisq(1 - alpha, 1L)
    resolved <- status == "WALD_COMPUTED"
    list(method = method, statistic = statistic, cutoff = cutoff, statistic_scale = "CHI_SQUARE_1", decision = if (!resolved) "UNRESOLVED" else if (statistic <=
        cutoff) "IN" else "OUT", status = status, diagnostic = FALSE, ...)
}

wrap_angle <- function(x) (x + pi)%%(2 * pi) - pi

circle_membership <- function(backend, q, cutoff, timeout = 120) {
    if (!length(q) %in% 1:2 || !is.numeric(q) || any(!is.finite(q)) || any(abs(q) > pi) || is.unsorted(q) ||
        length(cutoff) != 1L || !is.finite(cutoff) ||
        cutoff < 0)
        stop("invalid circle query")
    {
        load_numerics()
        if (all(backend$counts == floor(backend$counts))) {
            ans <- tryCatch(circle_membership_cpp(rep(backend$raw_observations, backend$counts), q[1L], cutoff,
                q[length(q)]), error = function(e) NULL)
            if (!is.null(ans) && ans$decision != "UNRESOLVED")
                return(membership_result(ans$decision, certified = TRUE, detail = paste("True-pi interval all-chart certificate; cells",
                  ans$r_cells)))
        }
    }
    if (length(q) == 2L && q[1L] != q[2L])
        return(membership_result("UNRESOLVED", detail = "Query interval could not be classified uniformly"))
    q <- q[1L]
    input <- tempfile()
    output <- tempfile()
    errors <- tempfile()
    on.exit(unlink(c(input, output, errors)), add = TRUE)
    fmt <- function(x) sprintf("%.17g", x)
    write.csv(data.frame(y = fmt(backend$raw_observations), count = backend$counts, q = fmt(q), cutoff = fmt(cutoff)),
        input, row.names = FALSE)
    status <- system2("python3", c("common/exact.py", "circle", shQuote(input)), stdout = output, stderr = errors,
        timeout = timeout)
    if (status != 0)
        return(membership_result("UNRESOLVED", detail = paste(readLines(errors), collapse = " ")))
    answer <- read.csv(output)
    if (nrow(answer) != 1L || !answer$decision %in% c("IN", "OUT", "UNRESOLVED"))
        stop("invalid circle certificate")
    membership_result(answer$decision, certified = answer$decision != "UNRESOLVED", detail = paste("Rational interval Bernstein certificate with mathematical pi; cells",
        answer$cells))
}

sample_circle_nonnegative <- function(n, lambda) {
    stopifnot(length(n) == 1L, is.finite(n), n == as.integer(n), n >= 1L)
    ab <- circle_nonnegative_coefficients(lambda)
    envelope <- 1 + sum(ab)
    out <- numeric(n)
    filled <- 0L
    while (filled < n) {
        size <- max(64L, as.integer(ceiling(1.1 * envelope * (n - filled))))
        proposal <- runif(size, -pi, pi)
        uniform <- runif(size)
        numerator <- 1 + ab[["a"]] * cos(proposal) + ab[["b"]] * cos(2 * proposal)
        accepted <- proposal[uniform * envelope <= numerator]
        take <- min(length(accepted), n - filled)
        if (take) {
            out[filled + seq_len(take)] <- accepted[seq_len(take)]
            filled <- filled + take
        }
    }
    out
}
