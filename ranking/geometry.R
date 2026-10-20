# Ranking: observed geometry, inference and simulation law.

calibrate_ranking <- function(observations, metric = c("euclidean", "kendall", "footrule"), draws, ...) {
    b <- ranking_discrete_backend(observations, match.arg(metric))
    if (nrow(draws) != 24L)
        stop("original 24-row ranking draws required")
    result <- calibrate_gl(backend = b, draws = draws[b$seen, , drop = FALSE], ...)
    result$original_draw_hash <- object_hash(draws)
    result
}

ranking_discrete_backend <- function(observations, metric = c("euclidean", "kendall", "footrule")) {
    metric <- match.arg(metric)
    geometry <- ranking_discrete_geometry(metric)
    if (!is.numeric(observations) || length(observations) < 3L || anyNA(observations) || any(!observations %in%
        seq_len(24L)))
        stop("ranking IDs 1:24 required")
    counts <- tabulate(observations, 24L)
    seen <- which(counts > 0L)
    list(domain = "finite_ranking", counts = counts[seen], seen = seen, observations = observations, metric = metric,
        loss = geometry$loss[seen, , drop = FALSE], pairs = geometry$pairs)
}

ranking_finite_bounds <- function(backend, draws, eps, tau) {
    loss <- backend$loss
    counts <- backend$counts
    pairs <- backend$pairs
    n <- sum(counts)
    if (any(loss != floor(loss)))
        stop("integer ranking loss required")
    delta <- loss[, pairs[, 1L], drop = FALSE] - loss[, pairs[, 2L], drop = FALSE]
    if (n^2 * max(abs(delta))^2 >= 2^53)
        stop("integer variance envelope exceeded")
    sum_delta <- drop(crossprod(counts, delta))
    variance_num <- n * drop(crossprod(counts, delta^2)) - sum_delta^2
    stopifnot(all(variance_num >= 0), all(variance_num == floor(variance_num)))
    risk_sums <- drop(crossprod(counts, loss))
    masks <- if (all(loss >= 0) && n * max(loss) < 2^53)
        ranking_selector_masks_cpp(risk_sums, sum_delta, variance_num, pairs, n, eps, tau)
    else rep(NA_integer_, nrow(pairs))
    unresolved <- which(is.na(masks))
    if (length(unresolved)) {
        pieces <- lapply(seq_len(ncol(loss)), function(j) list(left = c(0, 0), right = c(0, 0), b = matrix(0, nrow(loss),
            2L), c = cbind(loss[, j], loss[, j])))
        points <- cbind(pairs[unresolved, 1L] - 1L, pairs[unresolved, 2L] - 1L, 0, 0)
        masks[unresolved] <- exact_selector_masks(pieces, counts, points, eps, tau)
    }
    positive <- variance_num > 0
    stopifnot(identical(bitwAnd(masks, 1L) != 0L, positive))
    z <- ranking_maxima_cpp(loss, counts, draws, pairs, masks, positive, sum_delta, variance_num)
    c(z, list(masks = masks, variance_numerator = variance_num, total_pairs = nrow(pairs), minimum = min(risk_sums)/n,
        selector_interval_certified = nrow(pairs) - length(unresolved), selector_exact_fallback = length(unresolved)))
}

permutations <- function(x) {
    if (length(x) == 1L)
        return(matrix(x, 1L))
    do.call(rbind, lapply(seq_along(x), function(j) cbind(x[j], permutations(x[-j]))))
}

ranking_finite_membership <- function(q_lower, q_upper, loss_lower, loss_upper, counts, cutoff) {
    if (!exists("ranking_membership_cpp", mode = "function"))
        stop("Numerical engine is not loaded")
    if (!is.finite(cutoff) || cutoff < 0)
        stop("finite fixed cutoff required")
    ans <- ranking_membership_cpp(as.numeric(q_lower), as.numeric(q_upper), as.matrix(loss_lower), as.matrix(loss_upper),
        as.numeric(counts), cutoff)
    class(ans) <- "membership_result"
    ans
}

# Enumerate the entire finite target, rather than only the population means used
# by the Coverage evaluation. G and L invert the SAME empirical statistic; only
# their multiplier cutoffs differ. The Local comparison selector is a subset of
# the Global selector, hence c_L <= c_G and C_L is a subset of C_G.
ranking_candidate_sets <- function(fit) {
    if (!inherits(fit, "gl_certificate") || !identical(fit$domain, "finite_ranking") ||
        !identical(fit$backend$domain, "finite_ranking"))
        stop("a certified finite-ranking fit is required")
    backend <- fit$backend
    load_numerics()
    geometry <- ranking_discrete_geometry(backend$metric)
    if (ncol(backend$loss) != 24L || !identical(unname(backend$loss),
        unname(geometry$loss[backend$seen, , drop = FALSE])))
        stop("the complete canonical 24-ranking loss table is required")
    check_counts(backend$counts)
    if (length(backend$counts) != nrow(backend$loss))
        stop("observed counts and ranking loss rows differ")
    lower <- fit$cutoff$lower
    upper <- fit$cutoff$upper
    if (!identical(names(lower), c("G", "L")) || !identical(names(upper), c("G", "L")) ||
        any(!is.finite(c(lower, upper))) || any(lower < 0) || any(lower > upper) ||
        lower["L"] > lower["G"] || upper["L"] > upper["G"])
        stop("finite ordered G/L cutoff enclosures are required")

    # The admissible binary states are (G,L)=(0,0),(1,0),(1,1). Filtering these
    # states separately at each candidate is sharper than subtracting aggregate
    # marginal counts. It does not assert that every combination across distinct
    # candidates is achievable by the common unknown pair of true cutoffs.
    states <- data.frame(G = c(0L, 1L, 1L), L = c(0L, 0L, 1L))
    states$G_only <- states$G - states$L
    states$L_only <- 0L
    states$intersection <- states$L
    states$symmetric_difference <- states$G_only
    set_names <- names(states)
    endpoint_cutoffs <- c(lower["G"], upper["G"], lower["L"], upper["L"])
    endpoint_names <- c("G_lower_cutoff", "G_upper_cutoff", "L_lower_cutoff", "L_upper_cutoff")
    rows <- lapply(seq_len(24L), function(q) {
        # This is the finite oracle used by confidence_set_contains, with the
        # numerical engine checked once per fit instead of hashing it 24 times.
        endpoints <- lapply(endpoint_cutoffs, function(cutoff) ranking_finite_membership(
            backend$loss[, q], backend$loss[, q], backend$loss, backend$loss, backend$counts, cutoff))
        raw <- vapply(endpoints, function(x) x$decision, character(1))
        certified <- vapply(endpoints, function(x) isTRUE(x$certified), logical(1))
        if (length(raw) != 4L || any(!raw %in% c("IN", "OUT", "UNRESOLVED")))
            stop("invalid finite membership oracle result")
        decided <- raw
        proof <- certified & raw != "UNRESOLVED"
        # Each comparator inequality is monotone in c. Carry a certified IN
        # upward and a certified OUT downward, including equal-cutoff queries.
        for (j in seq_len(4L)) {
            inside <- any(proof & raw == "IN" & endpoint_cutoffs <= endpoint_cutoffs[j])
            outside <- any(proof & raw == "OUT" & endpoint_cutoffs >= endpoint_cutoffs[j])
            if (inside && outside)
                stop("finite candidate cutoff-monotonicity conflict")
            decided[j] <- if (inside) "IN" else if (outside) "OUT" else "UNRESOLVED"
        }
        definite_in <- decided[c(1L, 3L)] == "IN"
        possible_in <- decided[c(2L, 4L)] != "OUT"
        possible <- states$G >= definite_in[1L] & states$G <= possible_in[1L] &
            states$L >= definite_in[2L] & states$L <= possible_in[2L]
        if (!any(possible))
            stop("finite candidate G/L nesting conflict")
        ans <- list(candidate_id = q, ranking_key = geometry$keys[q])
        for (j in seq_len(4L)) {
            ans[[paste0(endpoint_names[j], "_decision")]] <- decided[j]
            ans[[paste0(endpoint_names[j], "_certified")]] <- decided[j] != "UNRESOLVED"
            ans[[paste0(endpoint_names[j], "_raw_decision")]] <- raw[j]
        }
        for (name in set_names) {
            lo <- min(states[[name]][possible])
            hi <- max(states[[name]][possible])
            ans[[paste0(name, "_lower")]] <- lo
            ans[[paste0(name, "_upper")]] <- hi
            ans[[paste0(name, "_decision")]] <- if (lo == 1L) "IN" else if (hi == 0L) "OUT" else "UNRESOLVED"
        }
        as.data.frame(ans, stringsAsFactors = FALSE, row.names = NULL)
    })
    candidates <- do.call(rbind, rows)
    counts <- do.call(rbind, lapply(set_names, function(name) {
        lo <- candidates[[paste0(name, "_lower")]]
        hi <- candidates[[paste0(name, "_upper")]]
        data.frame(set = name, lower = sum(lo), upper = sum(hi), unresolved = sum(lo != hi),
            resolved = sum(lo == hi), total_candidates = 24L)
    }))
    calibration_candidates <- NULL
    masks <- fit$metadata$masks
    if (!is.null(masks)) {
        if (length(masks) != 552L || anyNA(masks) || any(masks != floor(masks)) ||
            any(!masks %in% c(0L, 1L, 3L)) || !identical(backend$pairs, geometry$pairs))
            stop("invalid exhaustive finite calibration selector metadata")
        global <- bitwAnd(as.integer(masks), 1L) != 0L
        local <- bitwAnd(as.integer(masks), 2L) != 0L
        if (any(local & !global))
            stop("finite Local selector is not nested in Global")
        calibration_candidates <- data.frame(method = c("G", "L"), total_ordered_pairs = 552L,
            feasible_ordered_pairs = c(sum(global), sum(local)),
            first_endpoint_candidates = c(length(unique(backend$pairs[global, 1L])),
                length(unique(backend$pairs[local, 1L]))),
            assurance = "SAVED_EXHAUSTIVE_SELECTOR_MASKS_NO_RECALIBRATION")
    }
    list(candidates = candidates, counts = counts,
        cutoffs = data.frame(method = c("G", "L"), lower = unname(lower), upper = unname(upper)),
        calibration_candidates = calibration_candidates,
        nesting = "L_SUBSET_G_FULL_FINITE_COMPARATOR_INVERSION",
        assurance = "ALL_24_CANDIDATES_CERTIFIED_ENDPOINTS_MONOTONICITY_AND_NESTING",
        scope = "FIXED_OBSERVATIONS_AND_SUPPLIED_DRAWS_TRUE_FINITE_B_CUTOFF_ENCLOSURES")
}


ranking_graph_backend <- function (graph, vertex_ids)
{
    N <- nrow(graph$rankings)
    if (!length(vertex_ids) || any(vertex_ids < 1L | vertex_ids > N | vertex_ids != floor(vertex_ids)))
        stop("invalid observed vertex IDs")
    count_all <- tabulate(vertex_ids, N)
    seen <- which(count_all > 0L)
    counts <- count_all[seen]
    p <- counts/sum(counts)
    pieces <- lapply(seq_len(nrow(graph$edges)), function(j) {
        u <- graph$edges[j, 1L]
        v <- graph$edges[j, 2L]
        a <- graph$distance[seen, u]
        s <- graph$distance[seen, v] - a
        if (any(abs(s) != 1))
            stop("bipartite adjacent-swap identity failed")
        list(left = 0, right = 1, b = 2 * a * s, c = a^2, id = j, kind = "ranking_edge")
    })
    list(pieces = pieces, p = p, n = length(vertex_ids), counts = counts, seen = seen, graph = graph, domain = "continuous_ranking_graph",
        assurance = "FLOAT64_REFERENCE")
}

ranking_discrete_distribution <- function(metric = c("euclidean", "kendall", "footrule")) {
    metric <- match.arg(metric)
    original <- ranking_distribution("contamination", epsilon = c(1, 2))
    stopifnot(sum(original$population$weights) == 48L, sum(original$population$weights == 9L) == 3L, sum(original$population$weights ==
        1L) == 21L)
    original$geometry <- list(space = "ranking_discrete", metric = metric, m = 4L)
    original
}

ranking_discrete_geometry <- function(metric = c("euclidean", "kendall", "footrule")) {
    metric <- match.arg(metric)
    graph <- ranking_graph(4L)
    orders <- graph$rankings
    orders <- matrix(as.integer(orders), nrow = 24L)
    ranks <- t(apply(orders, 1L, function(z) match(1:4, z)))
    loss <- if (metric == "kendall")
        graph$distance^2
    else outer(seq_len(24L), seq_len(24L), Vectorize(function(i, j) if (metric == "footrule")
        sum(abs(ranks[i, ] - ranks[j, ]))^2
    else sum((ranks[i, ] - ranks[j, ])^2)))
    stopifnot(all(loss == floor(loss)), identical(dim(loss), c(24L, 24L)), all(diag(loss) == 0), all(loss == t(loss)))
    pairs <- which(row(loss) != col(loss), arr.ind = TRUE)
    list(keys = graph$key, loss = loss, pairs = pairs)
}

ranking_discrete_population_means <- function(metric = c("euclidean", "kendall", "footrule")) {
    metric <- match.arg(metric)
    g <- ranking_discrete_geometry(metric)
    dgp <- ranking_discrete_distribution(metric)
    w <- dgp$population$weights
    numerator <- drop(crossprod(w, g$loss))
    means <- which(numerator == min(numerator))
    list(mean_ids = means)
}

ranking_exact <- function(operation, tables, timeout = 120) {
    python <- Sys.which("python3")
    if (!nzchar(python))
        stop("python3 required for rational oracle")
    paths <- vapply(tables, function(x) tempfile(fileext = ".csv"), character(1))
    output <- tempfile()
    errors <- tempfile()
    on.exit(unlink(c(paths, output, errors)), add = TRUE)
    for (i in seq_along(tables)) write_exact_input(tables[[i]], paths[i])
    status <- run_exact_solver(python, c("common/exact.py", "ranking", operation, shQuote(paths)), output, errors,
        timeout)
    if (status != 0)
        stop("ranking rational solver failed: ", paste(readLines(errors, warn = FALSE), collapse = " "))
    read.csv(output, stringsAsFactors = FALSE, colClasses = "character")
}

ranking_graph <- local({
    cache <- vector("list", 3L)
    function(m = 4L) {
        if (!m %in% 3:5)
            stop("reference graph restricted to 3--5 objects")
        if (!is.null(cache[[m - 2L]]))
            return(cache[[m - 2L]])
        load_numerics()
        ranks <- permutations(seq_len(m))
        N <- nrow(ranks)
        key <- apply(ranks, 1L, paste, collapse = ",")
        distance <- ranking_distances_cpp(ranks)
        edges <- matrix(integer(), 0L, 2L)
        for (i in seq_len(N)) for (k in seq_len(m - 1L)) {
            r <- ranks[i, ]
            r[c(k, k + 1L)] <- r[c(k + 1L, k)]
            j <- match(paste(r, collapse = ","), key)
            if (i < j)
                edges <- rbind(edges, c(i, j))
        }
        cache[[m - 2L]] <<- list(rankings = ranks, distance = distance, edges = edges, key = key)
        cache[[m - 2L]]
    }
})

ranking_distribution <- function(family = "contamination", epsilon = c(1, 2)) {
    stopifnot(identical(family, "contamination"), identical(epsilon, c(1, 2)))
    graph <- ranking_graph(4L)
    weights <- rep(1, 24)
    weights[match(c("1,2,3,4", "2,3,1,4", "3,1,2,4"), graph$key)] <- 9
    probabilities <- weights/sum(weights)
    sampler <- local({
        p <- probabilities
        function(n) {
            check_sample_size(n)
            sample.int(24L, n, TRUE, p)
        }
    })
    list(family = family, parameter = epsilon, geometry = list(space = "ranking_graph", m = 4L), sample = sampler,
        population = list(weights = weights, probabilities = probabilities))
}

evaluate_graph_mean_set <- function(backend, oracle, cutoffs) {
    stopifnot(identical(backend$domain, "continuous_ranking_graph"), length(backend$pieces) == 36L, nrow(oracle$points) >
        0L, identical(names(cutoffs), c("G", "L")), all(is.finite(cutoffs)), all(cutoffs >= 0))
    rows <- do.call(rbind, lapply(backend$pieces, function(x) data.frame(edge = x$id, count = backend$counts, b = x$b,
        c = x$c)))
    ans <- ranking_exact("evaluate", list(rows, oracle$points, data.frame(method = names(cutoffs), cutoff = unname(cutoffs))))
    stopifnot(identical(ans$method, names(cutoffs)), all(ans$decision %in% c("IN", "OUT")))
    ans
}

ranking_graph_population_means <- function(dgp) {
    graph <- ranking_graph(dgp$geometry$m)
    w <- dgp$population$weights
    W <- sum(w)
    if (length(w) != 24L || any(!is.finite(w)) || any(w < 0) || any(w != floor(w)) || W <= 0 || W * 36 > 2^52)
        stop("invalid exact population weights")
    rows <- do.call(rbind, lapply(seq_len(nrow(graph$edges)), function(j) {
        u <- graph$edges[j, 1]
        v <- graph$edges[j, 2]
        a <- graph$distance[, u]
        s <- graph$distance[, v] - a
        stopifnot(all(abs(s) == 1))
        data.frame(edge = j, u = u, v = v, W = W, A = sum(w * a * a), B = sum(w * 2 * a * s))
    }))
    all_edges <- ranking_exact("oracle", list(rows))
    selected <- all_edges[all_edges$minimum == "1", , drop = FALSE]
    keys <- ifelse(selected$numerator == "0", paste0("v", selected$u), ifelse(selected$numerator == selected$denominator,
        paste0("v", selected$v), paste0("e", selected$edge)))
    points <- selected[!duplicated(keys), c("edge", "numerator", "denominator"), drop = FALSE]
    list(points = points)
}
