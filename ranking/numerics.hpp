#pragma once
namespace finite_kernel {
IntegerVector ranking_selector_masks_cpp(NumericVector risk, NumericVector sums,
                                         NumericVector variances, IntegerMatrix pairs, double n,
                                         double eps, double tau) {
    fc::check_platform();
    using namespace gl;
    int P = pairs.nrow();
    if (pairs.ncol() != 2 || sums.size() != P || variances.size() != P || risk.size() < 1 ||
        !R_finite(n) || n < 1 || n != std::floor(n) || !R_finite(eps) || eps < 0 ||
        !R_finite(tau) || tau <= 0)
        stop("invalid finite selector dimensions/controls");
    auto integer = [](double x) {
        return R_finite(x) && std::abs(x) < 9007199254740992.0 && x == std::floor(x);
    };
    for (double x : risk)
        if (!integer(x) || x < 0)
            stop("exact nonnegative risk sums required");
    double minimum = *std::min_element(risk.begin(), risk.end());
    auto leq = [](I a, I b) { return a.hi <= b.lo ? 1 : (a.lo > b.hi ? 0 : -1); };
    IntegerVector masks(P);
    for (int h = 0; h < P; ++h) {
        int q = pairs(h, 0) - 1, r = pairs(h, 1) - 1;
        if (q < 0 || r < 0 || q >= risk.size() || r >= risk.size() || !integer(sums[h]) ||
            !integer(variances[h]) || variances[h] < 0)
            stop("invalid exact finite moments");
        if (variances[h] == 0)
            continue;
        int anchor = leq(pt(risk[q] - minimum), mul(pt(n), pt(eps)));
        int slack =
            sums[h] >= 0 ? 1 : leq(mul(pt(n), sq(pt(sums[h]))), mul(sq(pt(tau)), pt(variances[h])));
        // Only separated interval comparisons are decisions. Exact boundary or
        // overflow ambiguity goes to the rational selector through NA_INTEGER.
        masks[h] = anchor == 0 || slack == 0 ? 1 : (anchor == 1 && slack == 1 ? 3 : NA_INTEGER);
    }
    return masks;
}

List ranking_maxima_cpp(NumericMatrix loss, NumericVector counts, NumericMatrix g,
                        IntegerMatrix pairs, IntegerVector masks, LogicalVector positive,
                        NumericVector sums, NumericVector variances) {
    fc::check_platform();
    using namespace gl;
    int s = loss.nrow(), B = g.ncol(), P = pairs.nrow();
    double n = 0;
    if (g.nrow() != s || pairs.ncol() != 2 || masks.size() != P || positive.size() != P ||
        sums.size() != P || variances.size() != P)
        stop("finite dimensions");
    for (double x : counts) {
        if (!R_finite(x) || x <= 0 || x != std::floor(x))
            stop("invalid count");
        n += x;
    }
    if (counts.size() != s || n > 4503599627370496.0)
        stop("count bound");
    for (double x : loss)
        if (!R_finite(x))
            stop("nonfinite loss");
    for (double x : g)
        if (!R_finite(x))
            stop("nonfinite draw");
    struct Pair {
        std::vector<double> centered;
        double variance;
        int q, r, sign, forward = 0, reverse = 0;
    };
    std::vector<Pair> grouped;
    std::map<std::vector<int64_t>, size_t> ids;
    for (int h = 0; h < P; ++h) {
        if (!positive[h])
            continue;
        int q = pairs(h, 0) - 1, r = pairs(h, 1) - 1;
        if (q < 0 || r < 0 || q >= loss.ncol() || r >= loss.ncol())
            stop("invalid finite pair");
        if (!R_finite(sums[h]) || std::abs(sums[h]) >= 0x1p53 || sums[h] != std::floor(sums[h]) ||
            !R_finite(variances[h]) || variances[h] <= 0 || variances[h] >= 0x1p53 ||
            variances[h] != std::floor(variances[h]))
            stop("exact finite moments required");
        std::vector<int64_t> key(s);
        int64_t divisor = 0;
        int sign = 0;
        for (int i = 0; i < s; ++i) {
            double x = n * (loss(i, q) - loss(i, r)) - sums[h];
            if (!R_finite(x) || std::abs(x) >= 0x1p53 || x != std::floor(x))
                stop("exact centered integer envelope exceeded");
            key[i] = static_cast<int64_t>(x);
            divisor = std::gcd(divisor, std::abs(key[i]));
            if (!sign && x != 0)
                sign = x > 0 ? 1 : -1;
        }
        if (!divisor)
            stop("positive variance with zero centered contrast");
        std::vector<double> centered(s);
        for (int i = 0; i < s; ++i) {
            centered[i] = sign * key[i];
            key[i] = sign * (key[i] / divisor);
        }
        auto entry = ids.emplace(std::move(key), grouped.size());
        if (entry.second)
            grouped.push_back({std::move(centered), variances[h], q, r, sign, 0, 0});
        auto &pair = grouped[entry.first->second];
        (sign > 0 ? pair.forward : pair.reverse) |= masks[h];
    }
    NumericMatrix lower(B, 2), upper(B, 2);
    std::vector<I> weighted(static_cast<size_t>(s) * B);
    for (int i = 0; i < s; ++i) {
        I mass = root(fc::ratio_positive(counts[i], n));
        for (int b = 0; b < B; ++b)
            weighted[static_cast<size_t>(b) * s + i] = mul_point(mass, g(i, b));
    }
    // A difference of candidate scores encloses every pair cheaply. Only pairs
    // that could beat a certified incumbent need the tighter direct dot product.
    std::vector<I> candidate(static_cast<size_t>(loss.ncol()) * B, pt(0));
    for (int q = 1; q < loss.ncol(); ++q) {
        double sum = 0;
        for (int i = 0; i < s; ++i)
            sum += counts[i] * (loss(i, q) - loss(i, 0));
        for (int i = 0; i < s; ++i) {
            double centered = n * (loss(i, q) - loss(i, 0)) - sum;
            for (int b = 0; b < B; ++b)
                candidate[static_cast<size_t>(q) * B + b] =
                    add(candidate[static_cast<size_t>(q) * B + b],
                        mul_point(weighted[static_cast<size_t>(b) * s + i], centered));
        }
    }
    for (size_t h = 0; h < grouped.size(); ++h) {
        if ((h & 63) == 0)
            checkUserInterrupt();
        const auto &pair = grouped[h];
        I sd = root(pt(pair.variance));
        I scale = {fc::down(1 / sd.hi), fc::up(1 / sd.lo)};
        // Exact centered integers and V eliminate cancellation in mean/variance.
        // Reuse weighted draws and normalize the whole dot product only once.
        for (int b = 0; b < B; ++b) {
            I screen = mul(mul_point(sub(candidate[static_cast<size_t>(pair.q) * B + b],
                                         candidate[static_cast<size_t>(pair.r) * B + b]),
                                     pair.sign),
                           scale);
            bool improves = false;
            for (int j = 0; j < 2; ++j)
                if (((pair.forward & (1 << j)) && screen.hi > lower(b, j)) ||
                    ((pair.reverse & (1 << j)) && -screen.lo > lower(b, j)))
                    improves = true;
            if (!improves)
                continue;
            I score = pt(0);
            for (int i = 0; i < s; ++i)
                score = add(score,
                            mul_point(weighted[static_cast<size_t>(b) * s + i], pair.centered[i]));
            score = mul(score, scale);
            // Reversing q,r negates the exact centered direction and its enclosure.
            // Both selector masks remain distinct; only the dot product is shared.
            for (int j = 0; j < 2; ++j) {
                if (pair.forward & (1 << j)) {
                    lower(b, j) = std::max(lower(b, j), score.lo);
                    upper(b, j) = std::max(upper(b, j), score.hi);
                }
                if (pair.reverse & (1 << j)) {
                    lower(b, j) = std::max(lower(b, j), -score.hi);
                    upper(b, j) = std::max(upper(b, j), -score.lo);
                }
            }
        }
    }
    return List::create(_["lower"] = lower, _["upper"] = upper,
                        _["evaluated_unsigned_pairs"] = grouped.size());
}

NumericMatrix ranking_distances_cpp(IntegerMatrix ranks) {
    int N = ranks.nrow(), m = ranks.ncol();
    if (N < 1 || m < 3 || m > 5)
        stop("invalid ranking dimensions");
    std::vector<std::array<int, 5>> positions(N);
    for (int i = 0; i < N; ++i) {
        std::array<bool, 5> seen{};
        for (int k = 0; k < m; ++k) {
            int a = ranks(i, k);
            if (a < 1 || a > m)
                stop("invalid permutation");
            --a;
            if (seen[a])
                stop("invalid permutation");
            seen[a] = true;
            positions[i][a] = k;
        }
    }
    NumericMatrix distance(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j <= i; ++j) {
            int count = 0;
            for (int k = 0; k < m; ++k)
                for (int l = k + 1; l < m; ++l)
                    count += positions[j][ranks(i, k) - 1] > positions[j][ranks(i, l) - 1];
            distance(i, j) = distance(j, i) = count;
        }
    return distance;
}

// The caller supplies intervals enclosing the TRUE loss values on the ENTIRE finite target.

Rcpp::List ranking_membership_cpp(Rcpp::NumericVector qlo, Rcpp::NumericVector qhi,
                                  Rcpp::NumericMatrix rlo, Rcpp::NumericMatrix rhi,
                                  Rcpp::NumericVector counts, double cutoff) {
    if (qlo.size() < 1 || rlo.ncol() < 1 || qlo.size() != qhi.size() ||
        qlo.size() != counts.size() || rlo.nrow() != qlo.size() || rhi.nrow() != rlo.nrow() ||
        rhi.ncol() != rlo.ncol())
        Rcpp::stop("incompatible dimensions");
    std::vector<fc::Interval> q(qlo.size()), r(qlo.size());
    std::vector<double> k(counts.begin(), counts.end());
    for (R_xlen_t i = 0; i < qlo.size(); ++i)
        q[i] = {qlo[i], qhi[i]};
    bool unresolved = false;
    for (int j = 0; j < rlo.ncol(); ++j) {
        if ((j & 127) == 0)
            Rcpp::checkUserInterrupt();
        for (R_xlen_t i = 0; i < qlo.size(); ++i)
            r[i] = {rlo(i, j), rhi(i, j)};
        fc::PairDecision result = fc::pair_decision(q, r, k, cutoff);
        if (result.decision == fc::Decision::Out)
            return Rcpp::List::create(
                Rcpp::Named("decision") = "OUT", Rcpp::Named("certified") = true,
                Rcpp::Named("witness") = j + 1, Rcpp::Named("assurance") = "FINITE_INPUT_INTERVALS",
                Rcpp::Named("margin") =
                    Rcpp::NumericVector::create(result.margin.lo, result.margin.hi));
        if (result.decision == fc::Decision::Unresolved)
            unresolved = true;
    }
    return Rcpp::List::create(
        Rcpp::Named("decision") = unresolved ? "UNRESOLVED" : "IN",
        Rcpp::Named("certified") = !unresolved, Rcpp::Named("assurance") = "FINITE_INPUT_INTERVALS",
        Rcpp::Named("scope") =
            "conditional on supplied loss enclosures, full finite target, and build contract");
}
} // namespace finite_kernel
