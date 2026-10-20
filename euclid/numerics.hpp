#pragma once
namespace finite_kernel {
List euclid_projection_cpp(NumericMatrix y, NumericMatrix g, IntegerVector columns,
                           Nullable<NumericMatrix> basis_lower = R_NilValue,
                           Nullable<NumericMatrix> basis_upper = R_NilValue) {
    fc::check_platform();
    using namespace gl;
    int n = y.nrow(), r = columns.size(), B = g.ncol();
    if (n < 2 || g.nrow() != n || B < 1)
        stop("invalid Euclidean draw dimensions");
    for (double x : y)
        if (!R_finite(x))
            stop("nonfinite observation");
    for (double x : g)
        if (!R_finite(x))
            stop("nonfinite draw");
    std::vector<std::vector<I>> z(n, std::vector<I>(r));
    if (basis_lower.isNotNull() != basis_upper.isNotNull())
        stop("both basis enclosures required");
    if (basis_lower.isNotNull()) {
        // These enclosures must come from the exact-input rational helper. An
        // arbitrary floating QR basis is not a certificate of the original space.
        NumericMatrix lo(basis_lower), hi(basis_upper);
        if (lo.nrow() != n || hi.nrow() != n || lo.ncol() != r || hi.ncol() != r)
            stop("basis dimensions");
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < r; ++j) {
                if (!R_finite(lo(i, j)) || !R_finite(hi(i, j)) || lo(i, j) > hi(i, j))
                    stop("invalid basis enclosure");
                z[i][j] = {lo(i, j), hi(i, j)};
            }
    } else
        for (int j = 0; j < r; ++j) {
            int k = columns[j] - 1;
            if (k < 0 || k >= y.ncol())
                stop("invalid support column");
            // Difference first: exact translation invariance and less cancellation.
            I mean = pt(0);
            for (int i = 0; i < n; ++i) {
                z[i][j] = sub(pt(y(i, k)), pt(y(0, k)));
                mean = add(mean, z[i][j]);
            }
            mean = div(mean, pt(n));
            for (int i = 0; i < n; ++i)
                z[i][j] = sub(z[i][j], mean);
        }
    std::vector<std::vector<I>> gram(r, std::vector<I>(r, pt(0))), inv;
    for (int j = 0; j < r; ++j)
        for (int k = 0; k < r; ++k)
            for (int i = 0; i < n; ++i)
                gram[j][k] = add(gram[j][k], mul(z[i][j], z[i][k]));
    bool verified = inverse(gram, inv);
    NumericVector lower(B), upper(B);
    for (int b = 0; b < B; ++b) {
        if ((b & 63) == 0)
            checkUserInterrupt();
        if (r == 0)
            continue;
        I gm = pt(0), cap = pt(0);
        for (int i = 0; i < n; ++i)
            gm = add(gm, pt(g(i, b)));
        gm = div(gm, pt(n));
        for (int i = 0; i < n; ++i)
            cap = add(cap, sq(sub(pt(g(i, b)), gm)));
        upper[b] = root(cap).hi;
        std::vector<I> h(r, pt(0));
        for (int j = 0; j < r; ++j)
            for (int i = 0; i < n; ++i)
                h[j] = add(h[j], mul_point(z[i][j], g(i, b)));
        // Safe one-column witnesses survive a failed inversion.
        for (int j = 0; j < r; ++j)
            if (gram[j][j].lo > 0) {
                I score = div(h[j], root(gram[j][j]));
                lower[b] = std::max(lower[b], std::max(score.lo, -score.hi));
            }
        if (verified) {
            I norm2 = pt(0);
            for (int j = 0; j < r; ++j)
                for (int k = 0; k < r; ++k)
                    norm2 = add(norm2, mul(h[j], mul(inv[j][k], h[k])));
            if (norm2.hi >= 0) {
                I v = root(norm2);
                lower[b] = std::max(lower[b], v.lo);
                upper[b] = std::min(upper[b], v.hi);
            }
        }
        if (lower[b] > upper[b])
            stop("Euclidean projection enclosure conflict");
    }
    return List::create(_["lower"] = lower, _["upper"] = upper, _["inverse_verified"] = verified);
}

} // namespace finite_kernel
