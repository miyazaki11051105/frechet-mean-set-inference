// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <cerrno>
#include <cstring>
// [[Rcpp::export(rng = false)]]
SEXP address_space_limits_cpp(SEXP requested = R_NilValue) {
    struct rlimit limits;
    if (getrlimit(RLIMIT_AS, &limits) != 0) Rcpp::stop("getrlimit: %s", std::strerror(errno));
    if (requested == R_NilValue)
        return Rcpp::NumericVector::create(limits.rlim_cur == RLIM_INFINITY ? R_PosInf : double(limits.rlim_cur),
                                          limits.rlim_max == RLIM_INFINITY ? R_PosInf : double(limits.rlim_max));
    if (TYPEOF(requested) != REALSXP || Rf_xlength(requested) != 2) Rcpp::stop("Two address-space limits required");
    rlim_t values[2];
    for (int i = 0; i < 2; ++i) {
        const double x = REAL(requested)[i];
        if (x == R_PosInf) values[i] = RLIM_INFINITY;
        else {
            if (!R_finite(x) || x < 0 || x != std::floor(x) || x > 9007199254740991.0)
                Rcpp::stop("Invalid address-space limit");
            values[i] = static_cast<rlim_t>(x);
        }
    }
    limits.rlim_cur = values[0]; limits.rlim_max = values[1];
    if (setrlimit(RLIMIT_AS, &limits) != 0) Rcpp::stop("setrlimit: %s", std::strerror(errno));
    return R_NilValue;
}
#include "common/interval.hpp"
using namespace Rcpp;
namespace finite_kernel {

#define FC_INLINE_BINARY64_STEPS

using namespace Rcpp;
namespace gl {
using I = fc::Interval;
I pt(double x) { return {x, x}; }
// Same enclosure as mul(a,{x,x}), omitting repeated endpoint products.
I mul_point(I a, double x) {
    if (x == 0 || (a.lo == 0 && a.hi == 0))
        return pt(0);
    double lo = (x < 0 ? a.hi : a.lo) * x, hi = (x < 0 ? a.lo : a.hi) * x;
    if (!R_finite(lo) || !R_finite(hi))
        return fc::whole();
    return {fc::down(lo), fc::up(hi)};
}
using fc::add;
using fc::mul;
using fc::sub;
I sq(I a) { return fc::square(a); }
I root(I a) {
    return {a.lo <= 0 ? 0 : std::max(0.0, fc::down(std::sqrt(a.lo))),
            fc::up(std::sqrt(std::max(0.0, a.hi)))};
}
I div(I a, I b) {
    if (b.lo <= 0)
        return fc::whole();
    return mul(a, {fc::down(1 / b.hi), fc::up(1 / b.lo)});
}
bool inverse(const std::vector<std::vector<I>> &gram, std::vector<std::vector<I>> &out) {
    int r = gram.size();
    auto a = gram;
    for (int i = 0; i < r; ++i) {
        a[i].resize(2 * r, pt(0));
        a[i][r + i] = pt(1);
    }
    for (int j = 0; j < r; ++j) {
        I p = a[j][j];
        if (!(p.lo > 0 && R_finite(p.hi)))
            return false;
        for (int k = 0; k < 2 * r; ++k)
            if (k != j)
                a[j][k] = div(a[j][k], p);
        a[j][j] = pt(1);
        for (int i = 0; i < r; ++i)
            if (i != j) {
                I f = a[i][j];
                for (int k = 0; k < 2 * r; ++k)
                    if (k != j)
                        a[i][k] = sub(a[i][k], mul(f, a[j][k]));
                a[i][j] = pt(0);
            }
    }
    out.assign(r, std::vector<I>(r));
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < r; ++j)
            out[i][j] = a[i][r + j];
    return true;
}
} // namespace gl

// Shared exact rank and interval-width helpers.

double cutoff_index_cpp(int B, double alpha) {
    fc::check_platform();
    if (B < 1 || !R_finite(alpha) || alpha <= 0 || alpha >= 1)
        stop("invalid rank controls");
#if defined(__SIZEOF_INT128__)
    int exponent;
    double fraction = std::frexp(alpha, &exponent);
    uint64_t mantissa = static_cast<uint64_t>(std::ldexp(fraction, 53)),
             N = static_cast<uint64_t>(B) + 1;
    __uint128_t product = static_cast<__uint128_t>(mantissa) * N;
    int shift = 53 - exponent;
    return static_cast<double>(N - (shift >= 128 ? 0 : static_cast<uint64_t>(product >> shift)));
#else
    return NA_REAL;
#endif
}

NumericVector interval_width_cpp(NumericVector lower, NumericVector upper) {
    fc::check_platform();
    if (lower.size() != upper.size())
        stop("width dimensions");
    NumericVector out(lower.size());
    for (int i = 0; i < lower.size(); ++i) {
        if (ISNAN(lower[i]) || ISNAN(upper[i]) || lower[i] > upper[i])
            stop("invalid width interval");
        out[i] = lower[i] == upper[i] ? 0 : fc::up(upper[i] - lower[i]);
    }
    return out;
}
} // namespace finite_kernel
namespace piecewise_kernel {
List circle_charts_cpp(NumericVector y);

// Outward-rounded pathwise bounds. Geometry is built without consulting draws.

#define FC_INLINE_BINARY64_STEPS

using namespace Rcpp;

namespace cone {
using I = fc::Interval;
I pt(double x) { return {x, x}; }
I neg(I a) { return {-a.hi, -a.lo}; }
I plus(I a, I b) { return fc::add(a, b); }
I minus(I a, I b) { return fc::sub(a, b); }
I times(I a, I b) { return fc::mul(a, b); }
// Identical enclosure to mul(a,{x,x}), with two endpoint products instead of
// four repeated products. Preserve zero, overflow and outward-rounding rules.
I times_point(I a, double x) {
    if (x == 0 || (a.lo == 0 && a.hi == 0))
        return {0, 0};
    double lo = (x < 0 ? a.hi : a.lo) * x, hi = (x < 0 ? a.lo : a.hi) * x;
    if (!R_finite(lo) || !R_finite(hi))
        return fc::whole();
    return {fc::down(lo), fc::up(hi)};
}
I sq(I a) { return fc::square(a); }
I root(I a) {
    if (a.hi < 0)
        stop("negative squared norm enclosure");
    return {a.lo <= 0 ? 0 : std::max(0.0, fc::down(std::sqrt(a.lo))),
            fc::up(std::sqrt(std::max(0.0, a.hi)))};
}
I divide(I a, I b) {
    if (b.lo <= 0)
        return fc::whole();
    return times(a, {fc::down(1 / b.hi), fc::up(1 / b.lo)});
}
I input(NumericVector x) {
    if (x.size() != 2 || !R_finite(x[0]) || !R_finite(x[1]) || x[0] > x[1])
        stop("invalid interval input");
    return {x[0], x[1]};
}
struct Piece {
    I left, right, mb, mc, vb;
    std::vector<I> b, c, bc, cc;
    int left_id = -1, right_id = -1;
    int circle_index = -1, circle_support = 0;
};
struct Box {
    int j, k, depth;
    double u0, u1, v0, v1;
};
using Basis = std::vector<std::vector<I>>;
// Raw features a^power 1{i <= cut}, centered with empirical weights.
// Five terms suffice for every circle comparison, independently of n.
struct PrefixTerm {
    int cut, power;
    I coefficient;
};
using PrefixVector = std::vector<PrefixTerm>;
struct CirclePrefix {
    std::vector<int> order;
    std::vector<I> lift, axis;
    std::vector<std::array<I, 3>> moment;
    std::vector<std::array<I, 4>> covariance;
    bool constant_sample = true;
    std::vector<std::vector<I>> gaussian;
    std::vector<std::vector<std::array<I, 2>>> gaussian_prefix;
    std::vector<I> gaussian_norm;
    std::vector<double> gaussian_cap;
    CirclePrefix(NumericVector y, const std::vector<I> &p, const std::vector<I> &weights)
        : axis(weights) {
        I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1};
        for (int i = 0; i < y.size(); ++i) {
            constant_sample = constant_sample && y[i] == y[0];
            order.push_back(i);
            lift.push_back(y[i] >= 0 ? minus(pt(y[i]), times(pt(2), pi)) : pt(y[i]));
        }
        std::sort(order.begin(), order.end(), [&](int i, int j) {
            return (y[i] >= 0) != (y[j] >= 0) ? y[i] >= 0 : y[i] < y[j];
        });
        moment.resize(y.size() + 1);
        for (size_t k = 0; k < order.size(); ++k) {
            int i = order[k];
            moment[k + 1] = moment[k];
            I power = pt(1);
            for (int r = 0; r < 3; ++r) {
                moment[k + 1][r] = plus(moment[k + 1][r], times(p[i], power));
                power = times(power, lift[i]);
            }
        }
        // The empirical weights sum to exactly one, despite their individual
        // interval enclosures. Constant centered features vanish symbolically.
        moment.back()[0] = pt(1);
        size_t N = moment.size();
        covariance.resize(N * N);
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                for (int a = 0; a < 2; ++a)
                    for (int b = 0; b < 2; ++b)
                        covariance[i * N + j][2 * a + b] =
                            (i == 0 || j == 0 || (i + 1 == N && a == 0) ||
                             (j + 1 == N && b == 0))
                                ? pt(0)
                                : minus(moment[std::min(i, j)][a + b],
                                        times(moment[i][a], moment[j][b]));
    }
    PrefixVector contrast(int e, int f, I t, I u) const {
        I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1}, fourpi = times(pt(4), pi);
        if (e == f) {
            I gap = minus(t, u);
            return {{static_cast<int>(order.size()), 1, times(pt(-2), gap)},
                    {e, 0, neg(times(fourpi, gap))}};
        }
        return {{static_cast<int>(order.size()), 1, times(pt(-2), minus(t, u))},
                {e, 0, times(fourpi, minus(pi, t))},
                {e, 1, fourpi},
                {f, 0, times(fourpi, minus(u, pi))},
                {f, 1, neg(fourpi)}};
    }
    PrefixVector displacement(int e, int f, I t, I u) const {
        I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1}, fourpi = times(pt(4), pi);
        return {{static_cast<int>(order.size()), 1, times(pt(-2), minus(t, u))},
                {e, 0, neg(times(fourpi, t))}, {f, 0, times(fourpi, u)}};
    }
    I gram(const PrefixVector &a, const PrefixVector &b) const {
        I out = pt(0);
        for (const auto &x : a)
            for (const auto &y : b) {
                I h = covariance[x.cut * moment.size() + y.cut][2 * x.power + y.power];
                out = plus(out, times(times(x.coefficient, y.coefficient), h));
            }
        return out;
    }
    std::vector<std::array<I, 2>> draws(const std::vector<I> &g) const {
        std::vector<std::array<I, 2>> out(order.size() + 1);
        for (size_t k = 0; k < order.size(); ++k) {
            int i = order[k];
            I z = times_point(axis[i], g[i].lo);
            out[k + 1][0] = plus(out[k][0], z);
            out[k + 1][1] = plus(out[k][1], times(z, lift[i]));
        }
        return out;
    }
    I dot(const PrefixVector &a, const std::vector<std::array<I, 2>> &g) const {
        I out = pt(0);
        for (const auto &x : a)
            if (x.cut && !(x.cut == static_cast<int>(order.size()) && x.power == 0))
            out = plus(out,
                       times(x.coefficient,
                             minus(g[x.cut][x.power], times(moment[x.cut][x.power], g.back()[0]))));
        return out;
    }
    void prepare_draws(NumericMatrix g) {
        if (g.nrow() != static_cast<int>(axis.size()))
            stop("circle draw dimensions");
        for (int b = 0; b < g.ncol(); ++b) {
            std::vector<I> draw(axis.size());
            I projection = pt(0), raw2 = pt(0), centered2 = pt(0);
            for (size_t i = 0; i < axis.size(); ++i) {
                if (!R_finite(g(i, b))) stop("nonfinite circle draw");
                draw[i] = pt(g(i, b));
                projection = plus(projection, times_point(axis[i], g(i, b)));
                raw2 = plus(raw2, sq(draw[i]));
            }
            for (size_t i = 0; i < axis.size(); ++i)
                centered2 = plus(centered2, sq(minus(draw[i], times(axis[i], projection))));
            gaussian_prefix.push_back(draws(draw));
            gaussian.push_back(std::move(draw));
            gaussian_norm.push_back(root(raw2));
            gaussian_cap.push_back(root(centered2).hi);
        }
    }
};
// A cell image is the convex hull of its four affine corner images. Its
// nonzero directions equal those of their cone; corners are generators,
// not a discretization of the optimization domain.
struct Cone {
    std::function<int(const std::array<double, 4> &, int)> radial;
    Basis rays;
    std::vector<PrefixVector> prefix_rays;
    std::vector<std::vector<I>> gram;
    struct Face {
        std::vector<int> ids;
        std::vector<std::vector<I>> inverse;
    };
    std::vector<Face> faces;
    bool face_envelope = false, enumeration_complete = true;
    std::vector<double> shift_weights;
    std::vector<I> shift_dot;
    std::function<Basis()> dense_builder;
    const CirclePrefix *prefix_model = nullptr;
    int structural_rank = 3;
    explicit Cone(Basis input, std::vector<PrefixVector> compressed = {},
                  const CirclePrefix *prefix = nullptr, bool verified_face_envelope = false,
                  int rank_bound = 3)
        : rays(std::move(input)), prefix_rays(std::move(compressed)),
          face_envelope(verified_face_envelope), prefix_model(prefix), structural_rank(rank_bound) {
        int m = rays.size();
        gram.assign(m, std::vector<I>(m, pt(0)));
        for (int j = 0; j < m; ++j)
            for (int k = 0; k <= j; ++k) {
                I h = prefix ? prefix->gram(prefix_rays[j], prefix_rays[k]) : pt(0);
                if (!rays[j].empty()) {
                    I dense = pt(0);
                    for (size_t a = 0; a < rays[j].size(); ++a)
                        dense = plus(dense, times(rays[j][a], rays[k][a]));
                    if (!prefix)
                        h = dense;
                    else {
                        h = {std::max(dense.lo, h.lo), std::min(dense.hi, h.hi)};
                        if (h.lo > h.hi) stop("inconsistent circle Gram enclosures");
                    }
                }
                gram[j][k] = gram[k][j] = h;
            }
        for (int mask = 1; mask < (1 << m); ++mask) {
            Face f;
            for (int j = 0; j < m; ++j)
                if (mask & (1 << j))
                    f.ids.push_back(j);
            int d = f.ids.size();
            if (d > rank_bound)
                continue;
            std::vector<std::vector<I>> a(d, std::vector<I>(2 * d, pt(0)));
            for (int j = 0; j < d; ++j) {
                for (int k = 0; k < d; ++k)
                    a[j][k] = gram[f.ids[j]][f.ids[k]];
                a[j][d + j] = pt(1);
            }
            bool valid = true;
            // Positive pivots certify invertibility of this true Gram submatrix.
            // Interval elimination encloses its exact inverse; no numerical rank cut.
            for (int j = 0; j < d; ++j) {
                I pivot = a[j][j];
                if (!(pivot.lo > 0 && R_finite(pivot.hi))) {
                    valid = false;
                    break;
                }
                for (int k = 0; k < 2 * d; ++k)
                    if (k != j)
                        a[j][k] = divide(a[j][k], pivot);
                a[j][j] = pt(1);
                for (int r = 0; r < d; ++r)
                    if (r != j) {
                        I factor = a[r][j];
                        for (int k = 0; k < 2 * d; ++k)
                            if (k != j)
                                a[r][k] = minus(a[r][k], times(factor, a[j][k]));
                        a[r][j] = pt(0);
                    }
            }
            if (valid) {
                f.inverse.assign(d, std::vector<I>(d));
                for (int j = 0; j < d; ++j)
                    for (int k = 0; k < d; ++k)
                        f.inverse[j][k] = a[j][d + k];
                faces.push_back(std::move(f));
            } else
                enumeration_complete = false;
        }
        if (face_envelope) {
            // A strictly positive dual direction also repairs uncertain KKT signs
            // when every inverse was verified: interval inverse widths can still
            // leave a possible-face envelope loose. No rank is inferred here.
            // Candidates are arbitrary binary64 coefficients; every ray' e > 0
            // is verified below, including redundant or nearly parallel rays.
            double best_margin = 0;
            for (const auto &f : faces) {
                std::vector<double> weights(m, 0);
                bool finite = true;
                for (size_t j = 0; j < f.ids.size(); ++j) {
                    I value = pt(0);
                    for (size_t k = 0; k < f.ids.size(); ++k)
                        value = plus(value, f.inverse[j][k]);
                    weights[f.ids[j]] = value.lo / 2 + value.hi / 2;
                    if (!R_finite(weights[f.ids[j]]))
                        finite = false;
                }
                if (!finite)
                    continue;
                std::vector<I> dots(m, pt(0));
                double margin = R_PosInf;
                I norm2 = pt(0);
                for (int j = 0; j < m; ++j) {
                    for (int k = 0; k < m; ++k)
                        dots[j] = plus(dots[j], times_point(gram[j][k], weights[k]));
                    margin = std::min(margin, dots[j].lo);
                    norm2 = plus(norm2, times_point(dots[j], weights[j]));
                }
                if (!(margin > 0 && norm2.hi > 0 && R_finite(norm2.hi)))
                    continue;
                double quality = margin / std::sqrt(norm2.hi); // work selection only
                if (quality > best_margin) {
                    best_margin = quality;
                    shift_weights = std::move(weights);
                    shift_dot = std::move(dots);
                }
            }
        }
    }
    void refine_geometry() {
        if (!dense_builder)
            return;
        Cone dense(dense_builder(), prefix_rays, prefix_model, face_envelope, structural_rank);
        rays = std::move(dense.rays);
        gram = std::move(dense.gram);
        faces = std::move(dense.faces);
        enumeration_complete = dense.enumeration_complete;
        shift_weights = std::move(dense.shift_weights);
        shift_dot = std::move(dense.shift_dot);
        dense_builder = {};
    }
    I optimize(const std::vector<I> &g, double cap, const I *known_dot = nullptr,
               std::array<double, 2> *selected = nullptr, int requested = 2) const {
        int m = rays.size();
        std::array<I, 4> h{};
        if (known_dot)
            for (int j = 0; j < m; ++j)
                h[j] = known_dot[j];
        else
            for (int j = 0; j < m; ++j)
                for (size_t a = 0; a < g.size(); ++a)
                    h[j] = plus(h[j], times_point(rays[j][a], g[a].lo));
        double lo = 0, hi = cap, possible_upper = 0;
        bool possible_face = false;
        std::array<double, 4> best_weights{};
        bool polar = true;
        for (int j = 0; j < m; ++j)
            if (h[j].hi > 0)
                polar = false;
        if (polar)
            return pt(0);
        for (const auto &f : faces) {
            int d = f.ids.size();
            std::array<I, 3> lambda{};
            bool positive = true;
            for (int j = 0; j < d; ++j) {
                for (int k = 0; k < d; ++k)
                    lambda[j] = plus(lambda[j], times(f.inverse[j][k], h[f.ids[k]]));
                if (lambda[j].lo < 0)
                    positive = false;
            }
            // Any nonnegative binary64 coefficients give a valid primal witness.
            std::array<double, 3> candidate{};
            for (int j = 0; j < d; ++j)
                candidate[j] = std::max(0.0, lambda[j].lo / 2 + lambda[j].hi / 2);
            bool finite = true;
            for (int j = 0; j < d; ++j)
                if (!R_finite(candidate[j]))
                    finite = false;
            if (!finite)
                continue;
            I numerator = pt(0), denominator = pt(0);
            for (int j = 0; j < d; ++j) {
                numerator = plus(numerator, times(pt(candidate[j]), h[f.ids[j]]));
                for (int k = 0; k < d; ++k)
                    denominator =
                        plus(denominator, times(pt(candidate[j]),
                                                times(pt(candidate[k]), gram[f.ids[j]][f.ids[k]])));
            }
            if (denominator.lo > 0) {
                double value = divide(numerator, root(denominator)).lo;
                if (value > lo) {
                    lo = value;
                    if (!shift_weights.empty()) {
                        I scale = divide(numerator, denominator);
                        double point = std::max(0.0, scale.lo / 2 + scale.hi / 2);
                        if (R_finite(point)) {
                            best_weights.fill(0);
                            for (int j = 0; j < d; ++j)
                                best_weights[f.ids[j]] = candidate[j] * point;
                        }
                    }
                }
                if (selected && radial && (requested & 2) && value > (*selected)[1]) {
                    std::array<double, 4> weights{};
                    for (int j = 0; j < d; ++j)
                        weights[f.ids[j]] = candidate[j];
                    int mask = radial(weights, requested) & requested;
                    for (int h = 1; h < 2; ++h)
                        if (mask & (1 << h))
                            (*selected)[h] = std::max((*selected)[h], value);
                }
            }
            if (!positive)
                continue;
            bool kkt = true;
            for (int j = 0; j < m; ++j)
                if (std::find(f.ids.begin(), f.ids.end(), j) == f.ids.end()) {
                    I residual = neg(h[j]);
                    for (int k = 0; k < d; ++k)
                        residual = plus(residual, times(gram[j][f.ids[k]], lambda[k]));
                    if (residual.lo < 0)
                        kkt = false;
                }
            if (kkt) {
                I norm2 = pt(0);
                for (int j = 0; j < d; ++j)
                    norm2 = plus(norm2, times(h[f.ids[j]], lambda[j]));
                if (norm2.hi >= 0)
                    hi = std::min(hi, root(norm2).hi);
            }
            // Optional remaining primal witnesses can be omitted once verified
            // endpoints are tight. Their actual gap is retained, never set to zero.
            // Calls needing selected Local witnesses still enumerate every face.
            if (face_envelope && (!selected || requested == 0) &&
                hi - lo <= 1e-10 * std::max(1.0, hi))
                return {lo, hi};
        }
        // Defer the boundary-face envelope until the ordinary KKT check has failed
        // to give a tight interval; resolved cones pay no extra face-enumeration cost.
        if (face_envelope && enumeration_complete && hi - lo > 1e-12 * std::max(1.0, hi)) {
            // The true projection has a nonnegative representation on an independent
            // face. ALL such faces through the structural rank bound have verified
            // inverses here. Keep every face whose KKT conditions remain possible.
            // One unresolved inverse disables this entire improvement.
            for (const auto &f : faces) {
                int d = f.ids.size();
                std::array<I, 3> lambda{};
                bool possible = true;
                for (int j = 0; j < d; ++j) {
                    for (int k = 0; k < d; ++k)
                        lambda[j] = plus(lambda[j], times(f.inverse[j][k], h[f.ids[k]]));
                    if (lambda[j].hi < 0)
                        possible = false;
                }
                for (int j = 0; j < m && possible; ++j)
                    if (std::find(f.ids.begin(), f.ids.end(), j) == f.ids.end()) {
                        I residual = neg(h[j]);
                        for (int k = 0; k < d; ++k)
                            residual = plus(residual, times(gram[j][f.ids[k]], lambda[k]));
                        if (residual.hi < 0)
                            possible = false;
                    }
                if (possible) {
                    I norm2 = pt(0);
                    for (int j = 0; j < d; ++j)
                        norm2 = plus(norm2, times(h[f.ids[j]], lambda[j]));
                    if (norm2.hi >= 0) {
                        possible_upper = std::max(possible_upper, root(norm2).hi);
                        possible_face = true;
                    }
                }
            }
            if (possible_face)
                hi = std::min(hi, possible_upper);
        }
        if (!shift_weights.empty() && hi - lo > 1e-12 * std::max(1.0, hi)) {
            // x=A*best_weights, residual r=g-x. Choose eta so A'(r-eta*e)<=0.
            // Thus z=r-eta*e is in the polar cone and, for every unit v in the
            // original cone, g'v <= (g-z)'v <= ||x+eta*e||. All signs and the norm
            // round outward, including the previously ambiguous boundary faces.
            double eta = 0;
            bool finite = true;
            for (int j = 0; j < m; ++j) {
                I residual = h[j];
                for (int k = 0; k < m; ++k)
                    residual = minus(residual, times_point(gram[j][k], best_weights[k]));
                if (residual.hi > 0)
                    eta = std::max(eta, fc::up(residual.hi / shift_dot[j].lo));
            }
            std::array<I, 4> coefficient{};
            for (int j = 0; j < m; ++j) {
                coefficient[j] = plus(pt(best_weights[j]), times_point(pt(shift_weights[j]), eta));
                if (!R_finite(coefficient[j].lo) || !R_finite(coefficient[j].hi))
                    finite = false;
            }
            if (finite) {
                I norm2 = pt(0);
                for (int j = 0; j < m; ++j)
                    for (int k = 0; k < m; ++k)
                        norm2 =
                            plus(norm2, times(coefficient[j], times(gram[j][k], coefficient[k])));
                if (norm2.hi >= 0)
                    hi = std::min(hi, root(norm2).hi);
            }
        }
        return {lo, hi};
    }
};
// Model owns a fixed-size basis table throughout draw_bounds; no per-cell copies.
struct LocalLimit {
    std::vector<I> normal;
    I norm2;
    I dual_y, dual_z;
};
struct CircleProfile {
    PrefixVector intercept, slope;
    std::array<I, 3> variance, delta;
    I sample_size, tau_squared;
    std::vector<I> boundaries;
    // Includes every possible constraint/variance root for the entire q strip.
    double upper(const CirclePrefix &, const std::vector<std::array<I, 2>> &,
                 double cap) const;
};
struct Model;
struct Direction {
    std::vector<I> w;
    int mask;
    bool sphere;
    const Basis *basis = nullptr;
    Basis rays;
    std::shared_ptr<Cone> cone;
    // Graph leaf cones depend only on immutable observed geometry. Rebuild a
    // single cone at its original draw-evaluation position instead of retaining
    // every leaf's dense rays, Gram matrix and face inverses simultaneously.
    Model *lazy_cone_model = nullptr;
    Box lazy_cone_box{};
    int cone_certain = 0;
    PrefixVector compact;
    const std::vector<PrefixVector> *compact_basis = nullptr;
    std::shared_ptr<LocalLimit> local;
    double local_eta = 0;
    std::shared_ptr<CircleProfile> profile;
};
struct Evaluation {
    std::vector<Direction> directions;
    int possible, certain;
    double width;
};
struct CoverPlan {
    int start, count, constructed;
    std::vector<Box> leaves;
};
struct Model {
    bool radial_feasibility = true;
    bool sharp_rectangles = false;
    std::shared_ptr<CirclePrefix> prefix;
    std::vector<Piece> pieces;
    std::vector<I> p, axis;
    std::vector<double> counts;
    I fmin, eps, tau2, n;
    int size;
    std::map<int, CoverPlan> cover_plans;
    size_t plan_bytes = 0;
    // std::map keeps element addresses stable; allocate only visited pair entries.
    std::map<size_t, Basis> bases;
    std::map<size_t, std::shared_ptr<LocalLimit>> local_limits;
    std::map<size_t, std::vector<PrefixVector>> compact_bases;
    // Fixed observed geometry only, bounded to at most 201^2 canonical roots.
    // Child cells and mutable draw-specific cones are never stored here.
    std::vector<std::shared_ptr<CircleProfile>> circle_root_profiles;
    std::vector<std::shared_ptr<const Evaluation>> circle_rank_root_evaluations;
    const Evaluation &rank_root_evaluation(const Box &box) {
        if (circle_rank_root_evaluations.empty())
            circle_rank_root_evaluations.resize(pieces.size()*pieces.size());
        size_t key=box.j*pieces.size()+box.k;
        if (!circle_rank_root_evaluations[key]) {
            circle_rank_root_evaluations[key]=std::make_shared<const Evaluation>(evaluate(box));
        }
        return *circle_rank_root_evaluations[key];
    }
    const std::vector<PrefixVector> &compressed_basis(int j, int k) {
        size_t key = j * pieces.size() + k;
        auto found = compact_bases.find(key);
        if (found != compact_bases.end())
            return found->second;
        auto &out = compact_bases[key];
        int e = pieces[j].circle_index, f = pieces[k].circle_index;
        int crossed = std::abs(e - f), s = size,
            d = (crossed == 0 || crossed == s) ? 1 : ((crossed == 1 || crossed == s - 1) ? 2 : 3);
        I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1}, fourpi = times(pt(4), pi);
        std::vector<PrefixVector> candidates = {{{s, 1, pt(-2)}, {e, 0, neg(fourpi)}}};
        if (d >= 2)
            candidates.push_back({{s, 1, pt(-2)}, {f, 0, neg(fourpi)}});
        if (d == 3)
            candidates.push_back({{e, 1, fourpi},
                                  {e, 0, times(fourpi, pi)},
                                  {f, 1, neg(fourpi)},
                                  {f, 0, neg(times(fourpi, pi))}});
        for (auto x : candidates) {
            for (const auto &q : out) {
                I h = prefix->gram(x, q);
                for (auto term : q) {
                    term.coefficient = neg(times(h, term.coefficient));
                    x.push_back(term);
                }
            }
            // Consolidate equal raw features; the interval coefficients still enclose
            // the exact Gram-Schmidt expression, without any rank tolerance.
            std::map<std::pair<int, int>, I> coefficients;
            for (auto term : x) {
                auto key = std::make_pair(term.cut, term.power);
                auto it = coefficients.find(key);
                if (it == coefficients.end())
                    coefficients[key] = term.coefficient;
                else
                    it->second = plus(it->second, term.coefficient);
            }
            x.clear();
            for (auto term : coefficients)
                x.push_back({term.first.first, term.first.second, term.second});
            I norm2 = prefix->gram(x, x);
            if (norm2.lo <= 0) {
                out.clear();
                return out;
            }
            I norm = root(norm2);
            for (auto &term : x)
                term.coefficient = divide(term.coefficient, norm);
            out.push_back(std::move(x));
        }
        return out;
    }
    std::map<size_t, int> constant;
    std::vector<int> exact_columns(const Basis &raw) const {
        std::vector<int> all(raw.size());
        std::iota(all.begin(), all.end(), 0);
        std::vector<std::vector<int64_t>> a(raw.size(), std::vector<int64_t>(size));
        for (size_t j = 0; j < raw.size(); ++j) {
            double sum = 0;
            for (int i = 0; i < size; ++i) {
                I x = raw[j][i];
                if (x.lo != x.hi || !R_finite(x.lo) || x.lo != std::floor(x.lo) ||
                    n.lo * std::abs(x.lo) > 0x1p50)
                    return all;
                sum += counts[i] * x.lo;
            }
            for (int i = 0; i < size; ++i) {
                double x = n.lo * raw[j][i].lo - sum;
                // Three-by-three minors have at most six products; this bound makes
                // every determinant exact in int64. Outside it, retain all columns.
                if (std::abs(x) > 0x1p19)
                    return all;
                a[j][i] = static_cast<int64_t>(x);
            }
        }
        std::vector<int> columns, rows;
        for (size_t j = 0; j < a.size(); ++j)
            for (int i = 0; i < size; ++i) {
                int64_t det = a[j][i];
                if (columns.size() == 1) {
                    const auto &u = a[columns[0]];
                    det = u[rows[0]] * a[j][i] - u[i] * a[j][rows[0]];
                }
                if (columns.size() == 2) {
                    const auto &u = a[columns[0]];
                    const auto &v = a[columns[1]];
                    int r = rows[0], s = rows[1];
                    det = u[r] * (v[s] * a[j][i] - v[i] * a[j][s]) -
                          v[r] * (u[s] * a[j][i] - u[i] * a[j][s]) +
                          a[j][r] * (u[s] * v[i] - u[i] * v[s]);
                }
                if (det != 0) {
                    columns.push_back(j);
                    rows.push_back(i);
                    break;
                }
            }
        return columns;
    }
    bool constant_contrast(int j, int k) {
        size_t key = j * pieces.size() + k;
        auto found = constant.find(key);
        if (found != constant.end())
            return found->second;
        constant[key] = 0;
        const auto &a = pieces[j], &b = pieces[k];
        std::pair<double, double> first;
        for (int i = 0; i < size; ++i) {
            if (a.b[i].lo != a.b[i].hi || b.b[i].lo != b.b[i].hi || a.b[i].lo != a.b[0].lo ||
                b.b[i].lo != b.b[0].lo || a.c[i].lo != a.c[i].hi || b.c[i].lo != b.c[i].hi)
                return false;
            double ac = a.c[i].lo, bc = b.c[i].lo;
            if (std::abs(ac) > 0x1p+500 || std::abs(bc) > 0x1p+500)
                return false;
            // Error-free TwoDiff: x+tail is the EXACT difference of the input doubles.
            // No rounded equality test of a single difference is used to drop variance.
            double x = ac - bc, bvirt = ac - x, avirt = x + bvirt, bround = bvirt - bc,
                   around = ac - avirt;
            std::pair<double, double> expansion = {x, around + bround};
            if (i == 0)
                first = expansion;
            else if (expansion != first)
                return false;
        }
        constant[key] = 1;
        return true;
    }
    const Basis &basis(int j, int k) {
        size_t key = j * pieces.size() + k;
        auto found = bases.find(key);
        if (found != bases.end())
            return found->second;
        bases.emplace(key, Basis{});
        auto &a = pieces[j], &b = pieces[k];
        bool shared = (a.left_id >= 0 && (a.left_id == b.left_id || a.left_id == b.right_id)) ||
                      (a.right_id >= 0 && (a.right_id == b.left_id || a.right_id == b.right_id));
        int dimension = shared ? 2 : 3;
        if (a.circle_index >= 0 && b.circle_index >= 0 && a.circle_support == b.circle_support) {
            int crossed = std::abs(a.circle_index - b.circle_index), total = a.circle_support;
            // Crossing all cuts shifts EVERY lift by 2*pi: one centered feature.
            // A singleton crossed set or complement adds only one indicator feature.
            if (crossed == 0 || crossed == total)
                dimension = 1;
            else if (crossed == 1 || crossed == total - 1)
                dimension = 2;
        }
        std::vector<std::vector<I>> candidates = {a.bc};
        if (dimension >= 2)
            candidates.push_back(b.bc);
        if (dimension == 3) {
            std::vector<I> cc(size);
            for (int i = 0; i < size; ++i)
                cc[i] = minus(a.cc[i], b.cc[i]);
            candidates.push_back(cc);
        }
        if (sharp_rectangles && !prefix) {
            Basis raw = {a.b};
            if (dimension >= 2)
                raw.push_back(b.b);
            if (dimension == 3) {
                std::vector<I> cc(size);
                for (int i = 0; i < size; ++i) {
                    double x = a.c[i].lo - b.c[i].lo;
                    cc[i] =
                        (a.c[i].lo == a.c[i].hi && b.c[i].lo == b.c[i].hi &&
                         a.c[i].lo == std::floor(a.c[i].lo) && b.c[i].lo == std::floor(b.c[i].lo) &&
                         std::abs(a.c[i].lo) <= 0x1p40 && std::abs(b.c[i].lo) <= 0x1p40)
                            ? pt(x)
                            : minus(a.c[i], b.c[i]);
                }
                raw.push_back(std::move(cc));
            }
            Basis independent;
            for (int j : exact_columns(raw))
                independent.push_back(candidates[j]);
            candidates = std::move(independent);
        }
        std::vector<std::vector<I>> orthogonal;
        for (auto x : candidates) {
            for (int i = 0; i < size; ++i)
                x[i] = times(axis[i], x[i]);
            for (const auto &q : orthogonal) {
                I dot = pt(0);
                for (int i = 0; i < size; ++i)
                    dot = plus(dot, times(q[i], x[i]));
                for (int i = 0; i < size; ++i)
                    x[i] = minus(x[i], times(dot, q[i]));
            }
            I norm2 = pt(0);
            for (I z : x)
                norm2 = plus(norm2, sq(z));
            if (norm2.lo <= 0)
                return bases[key]; // Never infer exact rank from a small numerical residual.
            I norm = root(norm2);
            for (I &z : x)
                z = divide(z, norm);
            orthogonal.push_back(x);
        }
        bases[key] = std::move(orthogonal);
        return bases[key];
    }
    I coordinate(const Piece &a, double lo, double hi) {
        // The SAME parameter appears twice. Use monotonicity of interpolation,
        // rather than treating (1-u) and u as independent intervals.
        I first = plus(times(a.left, pt(1 - lo)), times(a.right, pt(lo)));
        I last = plus(times(a.left, pt(1 - hi)), times(a.right, pt(hi)));
        return {first.lo, last.hi};
    }
    bool circle_shared_endpoint(int j, int k, int &ae, int &be) {
        if (!prefix || j == k) return false;
        const auto &a = pieces[j], &b = pieces[k];
        for (int x = 0; x < 2; ++x) for (int y = 0; y < 2; ++y) {
            int ai = x ? a.right_id : a.left_id, bi = y ? b.right_id : b.left_id;
            if (ai >= 0 && ai == bi) { ae = x; be = y; return true; }
        }
        return false;
    }
    I circle_endpoint_offset(const Piece &a, double lo, double hi, int endpoint) {
        // q-e = width*u at a left endpoint and -width*(1-u) at a right
        // endpoint. Avoid subtracting independently enclosed absolute angles.
        I width = minus(a.right, a.left);
        width = {std::max(0.0, width.lo), std::max(0.0, width.hi)};
        I fraction = endpoint ? minus(pt(1), {lo, hi}) : I{lo, hi};
        I distance = times(width, fraction);
        return endpoint ? neg(distance) : distance;
    }
    bool circle_endpoint_offsets(const Box &box, I &dq, I &dr, int &ae, int &be) {
        if (!circle_shared_endpoint(box.j, box.k, ae, be)) return false;
        dq = circle_endpoint_offset(pieces[box.j], box.u0, box.u1, ae);
        dr = circle_endpoint_offset(pieces[box.k], box.v0, box.v1, be);
        return true;
    }
    PrefixVector circle_contrast(const Box &box) {
        const auto &a = pieces[box.j], &b = pieces[box.k];
        I dq, dr; int ae, be;
        // Equal endpoint IDs certify identical squared losses, including the
        // two real-coordinate representatives -pi and +pi. Centering removes
        // each chart's scalar quadratic term, leaving this exact linear image.
        if (circle_endpoint_offsets(box, dq, dr, ae, be))
            return prefix->displacement(a.circle_index, b.circle_index, dq, dr);
        return prefix->contrast(a.circle_index, b.circle_index,
            coordinate(a, box.u0, box.u1), coordinate(b, box.v0, box.v1));
    }
    I circle_offset_risk_difference(int j, int k, I dq, I dr, int ae, int be) {
        const auto &a = pieces[j], &b = pieces[k];
        I ea = ae ? a.right : a.left, eb = be ? b.right : b.left;
        return plus(minus(times(dq, plus(times_point(ea, 2), a.mb)),
                          times(dr, plus(times_point(eb, 2), b.mb))),
                    minus(sq(dq), sq(dr)));
    }
    I circle_risk_difference(const Box &box) {
        I dq, dr; int ae, be;
        if (circle_endpoint_offsets(box, dq, dr, ae, be))
            return circle_offset_risk_difference(box.j, box.k, dq, dr, ae, be);
        return minus(risk(pieces[box.j], coordinate(pieces[box.j], box.u0, box.u1)),
                     risk(pieces[box.k], coordinate(pieces[box.k], box.v0, box.v1)));
    }
    std::vector<I> circle_dense_contrast(const Box &box) {
        const auto &a = pieces[box.j], &b = pieces[box.k];
        I dq, dr; int ae, be;
        bool shared = circle_endpoint_offsets(box, dq, dr, ae, be);
        I q = coordinate(a, box.u0, box.u1), r = coordinate(b, box.v0, box.v1);
        std::vector<I> w(size);
        for (int i = 0; i < size; ++i)
            w[i] = times(axis[i], shared ? minus(times(a.bc[i], dq), times(b.bc[i], dr))
                : plus(minus(times(a.bc[i], q), times(b.bc[i], r)), minus(a.cc[i], b.cc[i])));
        return w;
    }
    std::shared_ptr<LocalLimit> local_limit(int j, int k, int ae, int be) {
        size_t key = 4 * (j * pieces.size() + k) + 2 * ae + be;
        auto found = local_limits.find(key);
        if (found != local_limits.end())
            return found->second;
        auto &out = local_limits[key];
        const auto &a = pieces[j];
        const auto &b = pieces[k];
        if (!sharp_rectangles) return out;
        if (prefix) {
            int ai = ae ? a.right_id : a.left_id, bi = be ? b.right_id : b.left_id;
            // Prefix metadata was checked against regenerated true-pi charts.
            // Matching endpoint IDs therefore prove identical squared losses.
            if (ai < 0 || ai != bi) return out;
        } else {
            if (a.left.lo != 0 || a.left.hi != 0 || a.right.lo != 1 || a.right.hi != 1 ||
                b.left.lo != 0 || b.left.hi != 0 || b.right.lo != 1 || b.right.hi != 1)
                return out;
            // Prove the common graph endpoint losses equal with bounded integers.
            for (int i = 0; i < size; ++i) {
                for (I x : {a.b[i], b.b[i], a.c[i], b.c[i]})
                    if (x.lo != x.hi || x.lo != std::floor(x.lo) || std::abs(x.lo) > 0x1p20)
                        return out;
                if (ae + ae * a.b[i].lo + a.c[i].lo != be + be * b.b[i].lo + b.c[i].lo)
                    return out;
            }
        }
        std::vector<I> y(size), z(size);
        I yy = pt(0), yz = pt(0), zz = pt(0);
        for (int i = 0; i < size; ++i) {
            y[i] = times_point(times(axis[i], a.bc[i]), ae ? -1 : 1);
            z[i] = times_point(times(axis[i], b.bc[i]), be ? 1 : -1);
            yy = plus(yy, sq(y[i]));
            yz = plus(yz, times(y[i], z[i]));
            zz = plus(zz, sq(z[i]));
        }
        I det = minus(times(yy, zz), sq(yz));
        if (det.lo <= 0)
            return out;
        I iy = divide(zz, det), iz = divide(yy, det), iyz = neg(divide(yz, det));
        I qa = prefix ? times_point(ae ? a.right : a.left, 2) : pt(2 * ae);
        I rb = prefix ? times_point(be ? b.right : b.left, 2) : pt(2 * be);
        I my = times_point(plus(qa, a.mb), ae ? -1 : 1);
        I mz = times_point(plus(rb, b.mb), be ? 1 : -1);
        I cy = plus(times(iy, my), times(iyz, mz)), cz = plus(times(iyz, my), times(iz, mz));
        auto limit = std::make_shared<LocalLimit>();
        limit->norm2 = pt(0);
        for (int i = 0; i < size; ++i) {
            I x = plus(times(cy, y[i]), times(cz, z[i]));
            limit->normal.push_back(x);
            limit->norm2 = plus(limit->norm2, sq(x));
        }
        limit->dual_y = root(iy);
        limit->dual_z = root(iz);
        out = limit;
        return out;
    }
    I risk(const Piece &a, I t) {
        I natural = plus(plus(sq(t), times(a.mb, t)), a.mc);
        if ((a.circle_index < 0 && !sharp_rectangles) || !R_finite(t.lo) || !R_finite(t.hi) ||
            !R_finite(a.mb.lo) || !R_finite(a.mb.hi))
            return natural;
        // At a point t there is no repeated-variable dependency to remove.
        if (t.lo == t.hi)
            return natural;
        // Circle piece risk is t^2+b*t+c. Convexity in t and affinity in b
        // give the exact extrema over the enclosing independent (t,b) box.
        // This is a relaxation of the true dependent coefficients, hence safe.
        auto value = [&](double x, double b) { return plus(sq(pt(x)), times(pt(b), pt(x))); };
        double high = -R_PosInf;
        for (double x : {t.lo, t.hi})
            high = std::max(high, value(x, x >= 0 ? a.mb.hi : a.mb.lo).hi);
        double low = R_PosInf;
        auto minimum_on = [&](double left, double right, double b) {
            if (left > right)
                return;
            low = std::min(low, value(left, b).lo);
            low = std::min(low, value(right, b).lo);
            // Even a subnormal b/2 may be unrepresentable. The outward vertex
            // interval guards the inclusion test; -b^2/4 is then enclosed directly.
            I vertex = times_point(pt(b), -0.5);
            if (vertex.hi >= left && vertex.lo <= right)
                low = std::min(low, neg(times(sq(pt(b)), pt(.25))).lo);
        };
        if (t.lo <= 0)
            minimum_on(t.lo, std::min(0.0, t.hi), a.mb.hi);
        if (t.hi >= 0)
            minimum_on(std::max(0.0, t.lo), t.hi, a.mb.lo);
        I sharp = plus({low, high}, a.mc);
        I out = {std::max(natural.lo, sharp.lo), std::min(natural.hi, sharp.hi)};
        if (!fc::valid(out))
            stop("inconsistent circle quadratic risk enclosure");
        return out;
    }
    std::shared_ptr<Cone> cell_cone(const Box &box) {
        const auto &a = pieces[box.j];
        const auto &b = pieces[box.k];
        Basis corners;
        std::vector<PrefixVector> compressed;
        std::vector<std::pair<I, I>> coordinates;
        std::vector<std::pair<double, double>> unit_coordinates;
        std::vector<std::pair<I, I>> endpoint_offsets;
        int shared_ae = 0, shared_be = 0;
        bool shared_circle = circle_shared_endpoint(box.j, box.k, shared_ae, shared_be);
        for (double u : {box.u0, box.u1})
            for (double v : {box.v0, box.v1}) {
                int ai = u == 0 ? a.left_id : (u == 1 ? a.right_id : -1),
                    bi = v == 0 ? b.left_id : (v == 1 ? b.right_id : -1);
                if (ai >= 0 && ai == bi)
                    continue; // Identical graph point: exactly zero.
                I t = coordinate(a, u, u), s = coordinate(b, v, v);
                std::vector<I> w(prefix ? 0 : size);
                if (!prefix)
                    for (int i = 0; i < size; ++i)
                        w[i] = times(axis[i], plus(minus(times(a.bc[i], t), times(b.bc[i], s)),
                                                   minus(a.cc[i], b.cc[i])));
                corners.push_back(std::move(w));
                coordinates.push_back({t, s});
                unit_coordinates.push_back({u, v});
                if (shared_circle)
                    endpoint_offsets.push_back({circle_endpoint_offset(a, u, u, shared_ae),
                                                 circle_endpoint_offset(b, v, v, shared_be)});
                if (prefix)
                    compressed.push_back(circle_contrast({box.j, box.k, box.depth, u, u, v, v}));
            }
        bool shared = (a.left_id >= 0 && (a.left_id == b.left_id || a.left_id == b.right_id)) ||
                      (a.right_id >= 0 && (a.right_id == b.left_id || a.right_id == b.right_id));
        int rank_bound = shared ? 2 : 3;
        if (a.circle_index >= 0 && b.circle_index >= 0 && a.circle_support == b.circle_support) {
            int crossed = std::abs(a.circle_index - b.circle_index), total = a.circle_support;
            if (crossed == 0 || crossed == total)
                rank_bound = 1;
            else if (crossed == 1 || crossed == total - 1)
                rank_bound = 2;
        }
        auto cone = std::make_shared<Cone>(std::move(corners), std::move(compressed), prefix.get(),
                                           sharp_rectangles, rank_bound);
        if (prefix)
            cone->dense_builder = [this, box, unit_coordinates]() {
                Basis out;
                for (auto uv : unit_coordinates)
                    out.push_back(circle_dense_contrast({box.j, box.k, box.depth,
                        uv.first, uv.first, uv.second, uv.second}));
                return out;
            };
        // For any observed support size, nonnegative cone coefficients give a
        // genuine point in the affine rectangle after division by their sum.
        // Certify Local at that point; never assume that an unrestricted G optimum
        // is selected. This also covers circle's rank-three affine contrasts.
        if (radial_feasibility && size >= 3) {
            cone->radial = [this, box, coordinates, endpoint_offsets, shared_circle, shared_ae, shared_be](const std::array<double, 4> &weights,
                                                    int requested) {
                I sum = pt(0), t = pt(0), u = pt(0), dq = pt(0), dr = pt(0);
                for (size_t j = 0; j < coordinates.size(); ++j) {
                    sum = plus(sum, pt(weights[j]));
                    t = plus(t, times_point(coordinates[j].first, weights[j]));
                    u = plus(u, times_point(coordinates[j].second, weights[j]));
                    if (shared_circle) {
                        dq = plus(dq, times_point(endpoint_offsets[j].first, weights[j]));
                        dr = plus(dr, times_point(endpoint_offsets[j].second, weights[j]));
                    }
                }
                if (sum.lo <= 0)
                    return 0;
                t = divide(t, sum);
                u = divide(u, sum);
                if (shared_circle) { dq = divide(dq, sum); dr = divide(dr, sum); }
                const auto &a = pieces[box.j];
                const auto &b = pieces[box.k];
                I variance = pt(0);
                if (prefix) {
                    auto w = shared_circle
                        ? prefix->displacement(a.circle_index, b.circle_index, dq, dr)
                        : prefix->contrast(a.circle_index, b.circle_index, t, u);
                    variance = prefix->gram(w, w);
                } else
                    for (int i = 0; i < size; ++i) {
                        I z = plus(minus(times(a.bc[i], t), times(b.bc[i], u)),
                                   minus(a.cc[i], b.cc[i]));
                        variance = plus(variance, times(p[i], sq(z)));
                    }
                I delta = shared_circle
                    ? circle_offset_risk_difference(box.j, box.k, dq, dr, shared_ae, shared_be)
                    : minus(risk(a, t), risk(b, u));
                return select(delta, variance, minus(risk(a, t), fmin),
                              true) &
                       requested;
            };
        }
        // With three positive support weights, two loss contrasts determine the
        // entire centered vector. Certify invertibility, never threshold its rank.
        if (radial_feasibility && size == 3 && !prefix) {
            I aa = minus(a.b[0], a.b[2]), ab = neg(minus(b.b[0], b.b[2]));
            I ba = minus(a.b[1], a.b[2]), bb = neg(minus(b.b[1], b.b[2]));
            I ca = minus(minus(a.c[0], a.c[2]), minus(b.c[0], b.c[2]));
            I cb = minus(minus(a.c[1], a.c[2]), minus(b.c[1], b.c[2]));
            I determinant = minus(times(aa, bb), times(ab, ba));
            if (determinant.hi < 0) {
                aa = neg(aa);
                ab = neg(ab);
                ca = neg(ca);
                determinant = neg(determinant);
            }
            if (determinant.lo > 0) {
                I t0 = divide(minus(times(ab, cb), times(bb, ca)), determinant);
                I u0 = divide(minus(times(ba, ca), times(aa, cb)), determinant);
                cone->radial = [this, box, coordinates, t0,
                                u0](const std::array<double, 4> &weights, int requested) {
                    I sum = pt(0), t = pt(0), u = pt(0);
                    for (size_t j = 0; j < coordinates.size(); ++j) {
                        sum = plus(sum, pt(weights[j]));
                        t = plus(t, times_point(coordinates[j].first, weights[j]));
                        u = plus(u, times_point(coordinates[j].second, weights[j]));
                    }
                    if (sum.lo <= 0)
                        return 0;
                    t = divide(t, sum);
                    u = divide(u, sum);
                    const auto &a = pieces[box.j];
                    const auto &b = pieces[box.k];
                    I tl = coordinate(a, box.u0, box.u0), tr = coordinate(a, box.u1, box.u1);
                    I ul = coordinate(b, box.v0, box.v0), ur = coordinate(b, box.v1, box.v1);
                    // Bounded witness search only: a failure is UNKNOWN, never infeasible.
                    // No upper bound or domain restriction depends on this search budget.
                    double maximum = R_PosInf;
                    auto restrict_radius = [&](I origin, I direction, I left, I right) {
                        if (direction.lo > 0)
                            maximum = std::min(maximum, divide(minus(right, origin), direction).hi);
                        else if (direction.hi < 0)
                            maximum =
                                std::min(maximum, divide(minus(origin, left), neg(direction)).hi);
                    };
                    restrict_radius(t0, minus(t, t0), tl, tr);
                    restrict_radius(u0, minus(u, u0), ul, ur);
                    if (!R_finite(maximum) || maximum <= 0)
                        return 0;
                    int mask = 0;
                    std::deque<std::pair<I, int>> work;
                    work.push_back({{0, maximum}, 0});
                    for (int tries = 0; !work.empty() && tries < 63; ++tries) {
                        auto node = work.front();
                        work.pop_front();
                        I r = node.first;
                        auto check = [&](I radius, bool certain) {
                            I x = plus(t0, times(radius, minus(t, t0))),
                              y = plus(u0, times(radius, minus(u, u0)));
                            if (certain
                                    ? (x.lo < tl.hi || x.hi > tr.lo || y.lo < ul.hi || y.hi > ur.lo)
                                    : (x.hi < tl.lo || x.lo > tr.hi || y.hi < ul.lo ||
                                       y.lo > ur.hi))
                                return 0;
                            I variance = pt(0);
                            for (int i = 0; i < size; ++i) {
                                I z = plus(minus(times(a.bc[i], x), times(b.bc[i], y)),
                                           minus(a.cc[i], b.cc[i]));
                                variance = plus(variance, times(p[i], sq(z)));
                            }
                            return select(minus(risk(a, x), risk(b, y)), variance,
                                          minus(risk(a, x), fmin), certain);
                        };
                        if (!(check(r, false) & requested))
                            continue;
                        double mid = r.lo / 2 + r.hi / 2;
                        if (mid > 0)
                            mask |= check(pt(mid), true);
                        if ((mask & requested) == requested)
                            break;
                        if (node.second < 8) {
                            work.push_back({{r.lo, mid}, node.second + 1});
                            work.push_back({{mid, r.hi}, node.second + 1});
                        }
                    }
                    return mask;
                };
            }
        }
        return cone;
    }
    int select(I delta, I variance, I anchor, bool certain) {
        if (certain ? variance.lo <= 0 : variance.hi <= 0)
            return 0;
        int mask = 1;
        I slack = minus(times(tau2, variance), times(n, sq(delta)));
        bool pass = certain ? (anchor.hi <= eps.lo && (delta.lo >= 0 || slack.lo >= 0))
                            : (anchor.lo <= eps.hi && (delta.hi >= 0 || slack.hi >= 0));
        if (pass)
            mask |= 2;
        return mask;
    }
    Evaluation circle_evaluation(const Box &box, bool upper_geometry);
    std::shared_ptr<CircleProfile> circle_profile(const Box &box);
    std::shared_ptr<CircleProfile> circle_shared_profile(const Box &box);
    Evaluation evaluate(const Box &box, bool upper_geometry = true, bool compact_geometry = true) {
        if (prefix && compact_geometry)
            return circle_evaluation(box, upper_geometry);
        const auto &a = pieces[box.j], &b = pieces[box.k];
        if (constant_contrast(box.j, box.k))
            return {{}, 0, 0, 0};
        if (box.j == box.k && a.left.lo == a.left.hi && a.right.lo == a.right.hi &&
            a.left.lo == a.right.lo)
            return {{}, 0, 0, 0};
        if (box.j == box.k && box.u0 == box.u1 && box.v0 == box.v1 && box.u0 == box.v0)
            return {{}, 0, 0, 0};
        I t = coordinate(a, box.u0, box.u1), s = coordinate(b, box.v0, box.v1);
        I delta, variance = pt(0);
        std::vector<I> z(size);
        I gap = minus(t, s);
        I endpoint_dq = pt(0), endpoint_dr = pt(0);
        int endpoint_ae = 0, endpoint_be = 0;
        bool shared_circle = circle_endpoint_offsets(box, endpoint_dq, endpoint_dr, endpoint_ae, endpoint_be);
        if (box.j == box.k) {
            // Exact symbolic cancellation of the diagonal factor; no variance ridge.
            delta = times(gap, plus(plus(t, s), a.mb));
            variance = times(sq(gap), a.vb);
        } else {
            delta = shared_circle
                ? circle_offset_risk_difference(box.j, box.k, endpoint_dq, endpoint_dr, endpoint_ae, endpoint_be)
                : minus(risk(a, t), risk(b, s));
            for (int i = 0; i < size; ++i) {
                z[i] = shared_circle ? minus(times(a.bc[i], endpoint_dq), times(b.bc[i], endpoint_dr))
                    : plus(minus(times(a.bc[i], t), times(b.bc[i], s)), minus(a.cc[i], b.cc[i]));
                variance = plus(variance, times(p[i], sq(z[i])));
            }
            if (sharp_rectangles) {
                // Norm of an affine vector is convex: its maximum on a rectangle
                // occurs at a corner. This encloses the true maximum, including
                // coefficient/pi rounding, and strengthens the Local slack exclusion.
                double corner_upper = 0;
                for (double u : {box.u0, box.u1})
                    for (double v : {box.v0, box.v1}) {
                        I x = coordinate(a, u, u), y = coordinate(b, v, v), at = pt(0);
                        I dq = pt(0), dr = pt(0); int ae = 0, be = 0;
                        bool shared = circle_endpoint_offsets({box.j, box.k, box.depth, u, u, v, v}, dq, dr, ae, be);
                        for (int i = 0; i < size; ++i) {
                            I w = shared ? minus(times(a.bc[i], dq), times(b.bc[i], dr))
                                : plus(minus(times(a.bc[i], x), times(b.bc[i], y)), minus(a.cc[i], b.cc[i]));
                            at = plus(at, times(p[i], sq(w)));
                        }
                        corner_upper = std::max(corner_upper, at.hi);
                    }
                variance.hi = std::min(variance.hi, corner_upper);
                if (variance.lo > variance.hi)
                    stop("corner variance enclosure conflict");
            }
        }
        variance.lo = std::max(0.0, variance.lo);
        I anchor = minus(risk(a, t), fmin);
        int possible = select(delta, variance, anchor, false),
            certain = select(delta, variance, anchor, true);
        if (!possible)
            return {{}, 0, 0, 0};
        std::vector<Direction> dirs;
        double width = 0;
        if (box.j == box.k && a.vb.lo > 0) {
            I sd = root(a.vb);
            int sign_possible = 0, sign_certain = 3;
            for (int sign : {-1, 1}) {
                if ((sign < 0 && box.u0 >= box.v1) || (sign > 0 && box.u1 <= box.v0))
                    continue;
                // Cancel the NONZERO diagonal factor in BOTH direction and selector.
                // Its sign remains essential: T=sqrt(n)*sign*(t+s+mb)/sqrt(vb).
                I scaled = times_point(plus(plus(t, s), a.mb), sign);
                bool factored = sharp_rectangles && !prefix;
                int pm = factored ? select(scaled, a.vb, anchor, false) : possible;
                int cm = factored ? select(scaled, a.vb, anchor, true) : certain;
                sign_possible |= pm;
                sign_certain &= cm;
                Direction d;
                d.mask = pm;
                d.sphere = false;
                for (int i = 0; i < size; ++i)
                    d.w.push_back(times(pt(sign), divide(times(axis[i], a.bc[i]), sd)));
                if (prefix) {
                    I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1};
                    d.compact = {{size, 1, divide(pt(-2 * sign), sd)},
                                 {a.circle_index, 0, divide(times(pt(-4 * sign), pi), sd)}};
                }
                dirs.push_back(d);
            }
            possible = sign_possible;
            certain = sign_certain & possible;
        } else if (variance.lo > 0) {
            Direction d;
            d.mask = possible;
            d.sphere = false;
            I sd = root(variance);
            for (int i = 0; i < size; ++i)
                d.w.push_back(divide(times(axis[i], z[i]), sd));
            if (prefix) {
                d.compact = circle_contrast(box);
                for (auto &term : d.compact)
                    term.coefficient = divide(term.coefficient, sd);
            }
            dirs.push_back(d);
        } else
            dirs.push_back({{}, possible, true});
        if (upper_geometry && box.j != box.k)
            for (auto &d : dirs) {
                d.basis = &basis(box.j, box.k);
                if (prefix) {
                    d.compact_basis = &compressed_basis(box.j, box.k);
                    if (d.compact_basis->empty())
                        d.compact_basis = nullptr;
                }
                for (int ae = 0; ae < 2; ++ae)
                    for (int be = 0; be < 2; ++be) {
                        int ai = ae ? a.right_id : a.left_id, bi = be ? b.right_id : b.left_id;
                        bool touches =
                            (ae ? box.u1 == 1 : box.u0 == 0) && (be ? box.v1 == 1 : box.v0 == 0);
                        if (ai >= 0 && ai == bi) {
                            if (possible & 2) {
                                d.local = local_limit(box.j, box.k, ae, be);
                                if (d.local) {
                                    I dq = prefix ? circle_endpoint_offset(a, box.u0, box.u1, ae) : pt(0);
                                    I dr = prefix ? circle_endpoint_offset(b, box.v0, box.v1, be) : pt(0);
                                    I remainder = plus(
                                        times_point(d.local->dual_y, prefix ? std::max(0.0, (ae ? neg(dq) : dq).hi) : (ae ? 1 - box.u0 : box.u1)),
                                        times_point(d.local->dual_z, prefix ? std::max(0.0, (be ? neg(dr) : dr).hi) : (be ? 1 - box.v0 : box.v1)));
                                    d.local_eta = plus(divide(root(tau2), root(n)), remainder).hi;
                                }
                            }
                            if (!touches)
                                continue;
                            std::vector<I> ar(size), br(size);
                            for (int i = 0; i < size; ++i) {
                                ar[i] = times(pt(ae ? -1 : 1), times(axis[i], a.bc[i]));
                                br[i] = times(pt(be ? 1 : -1), times(axis[i], b.bc[i]));
                            }
                            d.rays = {ar, br};
                        }
                    }
            }
        for (auto &d : dirs) {
            if (d.sphere) {
                width = 2;
                continue;
            }
            double sum = 0;
            for (auto &w : d.w) {
                w.lo = std::max(-1.0, w.lo);
                w.hi = std::min(1.0, w.hi);
                if (w.lo > w.hi)
                    stop("empty direction enclosure");
                sum += (w.hi - w.lo) * (w.hi - w.lo);
            }
            width = std::max(width, std::sqrt(sum));
        }
        return {dirs, possible, certain, width};
    }
};
Model model(List source, NumericVector counts, double eps, double tau) {
    Model m;
    m.size = counts.size();
    double n = 0;
    for (double x : counts) {
        if (!R_finite(x) || x <= 0 || std::floor(x) != x)
            stop("invalid counts");
        n += x;
    }
    if (m.size < 1 || n > 4503599627370496.0)
        stop("invalid count total");
    m.n = pt(n);
    m.eps = pt(eps);
    m.tau2 = sq(pt(tau));
    m.counts = as<std::vector<double>>(counts);
    for (double x : counts) {
        m.p.push_back(fc::ratio_positive(x, n));
        m.axis.push_back(root(m.p.back()));
    }
    for (SEXP el : source) {
        List a(el);
        Piece x;
        x.left = input(a["left"]);
        x.right = input(a["right"]);
        if (a.containsElementNamed("left_id"))
            x.left_id = as<int>(a["left_id"]);
        if (a.containsElementNamed("right_id"))
            x.right_id = as<int>(a["right_id"]);
        if (a.containsElementNamed("circle_index"))
            x.circle_index = as<int>(a["circle_index"]);
        if (a.containsElementNamed("circle_support"))
            x.circle_support = as<int>(a["circle_support"]);
        if (x.left.lo > x.right.hi)
            stop("invalid piece domain");
        NumericMatrix bb(a["b"]), cc(a["c"]);
        if (bb.nrow() != m.size || cc.nrow() != m.size || bb.ncol() != 2 || cc.ncol() != 2)
            stop("coefficient dimensions");
        x.mb = x.mc = pt(0);
        for (int i = 0; i < m.size; ++i) {
            x.b.push_back(input(NumericVector::create(bb(i, 0), bb(i, 1))));
            x.c.push_back(input(NumericVector::create(cc(i, 0), cc(i, 1))));
            x.mb = plus(x.mb, times(m.p[i], x.b.back()));
            x.mc = plus(x.mc, times(m.p[i], x.c.back()));
        }
        x.vb = pt(0);
        for (int i = 0; i < m.size; ++i) {
            x.bc.push_back(minus(x.b[i], x.mb));
            x.cc.push_back(minus(x.c[i], x.mc));
            x.vb = plus(x.vb, times(m.p[i], sq(x.bc.back())));
        }
        x.vb.lo = std::max(0.0, x.vb.lo);
        m.pieces.push_back(x);
    }
    if (m.pieces.empty())
        stop("empty domain");
    m.fmin = {R_PosInf, R_PosInf};
    for (auto &a : m.pieces) {
        I vertex = times(pt(-.5), a.mb);
        I minimizer = {std::min(a.right.lo, std::max(a.left.lo, vertex.lo)),
                       std::min(a.right.hi, std::max(a.left.hi, vertex.hi))};
        I val = m.risk(a, minimizer);
        m.fmin.lo = std::min(m.fmin.lo, val.lo);
        m.fmin.hi = std::min(m.fmin.hi, val.hi);
    }
    return m;
}
// A floating candidate polar vector is only a proposal. Every dual inequality
// is then certified by interval dot products; failure returns the prior bound.
double cone_bound(const std::vector<std::vector<I>> &rays, const std::vector<I> &draw,
                  double prior) {
    if (rays.size() != 2)
        return prior;
    int s = draw.size();
    std::vector<double> a(s), b(s), h(s), p(s), best(s, 0);
    double aa = 0, bb = 0, ab = 0, ag = 0, bg = 0, bestnorm = 0;
    for (int i = 0; i < s; ++i) {
        a[i] = (rays[0][i].lo + rays[0][i].hi) / 2;
        b[i] = (rays[1][i].lo + rays[1][i].hi) / 2;
        aa += a[i] * a[i];
        bb += b[i] * b[i];
        ab += a[i] * b[i];
        ag += a[i] * draw[i].lo;
        bg += b[i] * draw[i].lo;
        bestnorm += draw[i].lo * draw[i].lo;
    }
    if (!(aa > 0 && bb > 0) || !R_finite(aa + bb + bestnorm))
        return prior;
    auto propose = [&](double x, double y) {
        if (!R_finite(x) || !R_finite(y) || x < 0 || y < 0)
            return;
        double residual = 0;
        std::vector<double> k(s);
        for (int i = 0; i < s; ++i) {
            k[i] = x * a[i] + y * b[i];
            double r = draw[i].lo - k[i];
            residual += r * r;
        }
        if (residual < bestnorm) {
            bestnorm = residual;
            best = k;
        }
    };
    propose(std::max(0.0, ag / aa), 0);
    propose(0, std::max(0.0, bg / bb));
    double det = aa * bb - ab * ab;
    if (det > 0)
        propose((ag * bb - bg * ab) / det, (bg * aa - ag * ab) / det);
    for (int i = 0; i < s; ++i) {
        p[i] = draw[i].lo - best[i];
        h[i] = a[i] / std::sqrt(aa) + b[i] / std::sqrt(bb);
    }
    double correction = 0;
    for (const auto &ray : rays) {
        I hp = pt(0), pp = pt(0);
        for (int i = 0; i < s; ++i) {
            hp = plus(hp, times(ray[i], pt(h[i])));
            pp = plus(pp, times(ray[i], pt(p[i])));
        }
        if (hp.lo <= 0)
            return prior;
        correction = std::max(correction, fc::up(pp.hi / hp.lo));
    }
    double shift = correction > 0 ? fc::up(2 * correction) : 0;
    for (int attempt = 0; attempt < 5; ++attempt) {
        std::vector<double> candidate(s);
        for (int i = 0; i < s; ++i)
            candidate[i] = p[i] - shift * h[i];
        bool feasible = true;
        for (const auto &ray : rays) {
            I value = pt(0);
            for (int i = 0; i < s; ++i)
                value = plus(value, times(ray[i], pt(candidate[i])));
            if (value.hi > 0 || !R_finite(value.hi))
                feasible = false;
        }
        if (feasible) {
            I residual = pt(0);
            for (int i = 0; i < s; ++i)
                residual = plus(residual, sq(minus(draw[i], pt(candidate[i]))));
            return std::min(prior, root(residual).hi);
        }
        shift = std::max(1e-14, shift * 4);
    }
    return prior;
}
std::vector<Direction> deduplicate_inner(const std::vector<Direction> &inner, int s) {
    // Coincident corners and parent/child representatives often yield the exact
    // same interval vector. Its dot-product lower bound is identical; for each
    // selected method an original feasible representative still witnesses it.
    // No distance tolerance, rounding, or change of the direction class is used.
    std::map<std::vector<double>, Direction> unique_inner;
    for (const auto &d : inner)
        if (!d.sphere) {
            std::vector<double> key;
            key.reserve(2 * s);
            for (I w : d.w) {
                key.push_back(w.lo);
                key.push_back(w.hi);
            }
            if (d.w.empty())
                for (const auto &term : d.compact) {
                    key.push_back(term.cut);
                    key.push_back(term.power);
                    key.push_back(term.coefficient.lo);
                    key.push_back(term.coefficient.hi);
                }
            auto inserted = unique_inner.emplace(std::move(key), d);
            if (!inserted.second)
                inserted.first->second.mask |= d.mask;
        }
    std::vector<Direction> result;
    result.reserve(unique_inner.size());
    for (auto &entry : unique_inner)
        result.push_back(std::move(entry.second));
    return result;
}
// The graph construction owns these representatives outright. Move their
// dense vectors into the same exact-key map instead of retaining both the
// original and a copied bank during deduplication. Map order, first witnesses
// and mask unions are identical to deduplicate_inner above.
std::vector<Direction> deduplicate_owned_inner(std::vector<Direction> inner, int s) {
    std::map<std::vector<double>, Direction> unique_inner;
    for (auto &d : inner) if (!d.sphere) {
        std::vector<double> key;
        key.reserve(2 * s);
        for (I w : d.w) { key.push_back(w.lo); key.push_back(w.hi); }
        if (d.w.empty()) for (const auto &term : d.compact) {
            key.push_back(term.cut); key.push_back(term.power);
            key.push_back(term.coefficient.lo); key.push_back(term.coefficient.hi);
        }
        auto found = unique_inner.find(key);
        if (found == unique_inner.end()) unique_inner.emplace(std::move(key), std::move(d));
        else found->second.mask |= d.mask;
    }
    std::vector<Direction>().swap(inner);
    std::vector<Direction> result;
    result.reserve(unique_inner.size());
    for (auto &entry : unique_inner) result.push_back(std::move(entry.second));
    return result;
}
// A group is used ONLY to avoid redundant lower-witness dot products. Every
// member remains feasible under its original mask; no outer cell is removed.
struct WitnessGroup {
    size_t first, last, seed;
    double radius;
    int mask;
};
std::vector<WitnessGroup> witness_groups(const std::vector<Direction> &inner,
                                        const CirclePrefix *prefix = nullptr) {
    std::vector<WitnessGroup> groups;
    if (inner.size() < 128 ||
        (!prefix && std::any_of(inner.begin(), inner.end(),
                                [](const Direction &d) { return d.w.empty(); })))
        return groups; // Small banks are cheaper to scan directly.
    for (size_t first = 0; first < inner.size(); first += 32) {
        size_t last = std::min(first + 32, inner.size()), seed = first + (last - first) / 2;
        double radius = 0;
        int mask = 0;
        for (size_t j = first; j < last; ++j) {
            mask |= inner[j].mask;
            if (j == seed)
                continue;
            I distance2 = pt(0);
            if (prefix && !inner[j].compact.empty() && !inner[seed].compact.empty()) {
                auto difference = inner[j].compact;
                for (auto term : inner[seed].compact) {
                    term.coefficient = neg(term.coefficient);
                    difference.push_back(term);
                }
                distance2 = prefix->gram(difference, difference);
            } else
                for (size_t i = 0; i < inner[j].w.size(); ++i)
                    distance2 = plus(distance2, sq(minus(inner[j].w[i], inner[seed].w[i])));
            radius = std::max(radius, std::min(2.0, root(distance2).hi));
        }
        groups.push_back({first, last, seed, radius, mask});
    }
    return groups;
}
List draw_bounds(const std::vector<Direction> &inner, const std::vector<Direction> &outer,
                 NumericMatrix g, const std::vector<I> &axis, bool already_unique = false,
                 const CirclePrefix *prefix = nullptr,
                 const std::vector<WitnessGroup> *prepared_groups = nullptr,
                 bool shared_lower = true, bool lower_only = false,
                 const NumericMatrix *prior_lower = nullptr,
                 const NumericMatrix *prior_upper = nullptr,
                 const IntegerVector *draw_ids = nullptr,
                 const NumericMatrix *cutoff_targets = nullptr, int covered_methods = 3) {
    int B = draw_ids ? draw_ids->size() : g.ncol();
    int s = draw_ids ? prefix->axis.size() : g.nrow();
    NumericMatrix lower(B, 2), upper(B, 2);
    // Internal complete-cover accumulator. Prior lower witnesses are valid on
    // the ORIGINAL full domain; prior uppers cover every previously visited
    // block. A cell may be skipped only if its certified upper is no larger
    // than this accumulated upper. Both independent G covers receive the same
    // accumulator, so intersecting them still covers the union of all blocks.
    if (prior_lower || prior_upper) {
        if (!prior_lower || !prior_upper || prior_lower->nrow() != B || prior_upper->nrow() != B ||
            prior_lower->ncol() != 2 || prior_upper->ncol() != 2)
            stop("invalid cover accumulator dimensions");
        for (int j = 0; j < 2; ++j)
            for (int b = 0; b < B; ++b) {
                double lo = (*prior_lower)(b, j), hi = (*prior_upper)(b, j);
                if (!R_finite(lo) || ISNAN(hi) || lo < 0 || lo > hi)
                    stop("invalid cover accumulator bounds");
                lower(b, j) = lo;
                upper(b, j) = hi;
            }
    }
    std::vector<Direction> local_unique;
    if (!already_unique)
        local_unique = deduplicate_inner(inner, s);
    const auto &unique_inner = already_unique ? inner : local_unique;
    std::vector<WitnessGroup> local_groups;
    if (shared_lower && !prepared_groups)
        local_groups = witness_groups(unique_inner, prefix);
    const auto &groups = prepared_groups ? *prepared_groups : local_groups;
    std::vector<I> seed_scores(groups.size());
    std::vector<int> seed_stamps(groups.size(), -1);
    double witness_dots = 0, witness_skipped = 0;
    // Give each shared subspace an index once, instead of allocating a map and
    // searching it for every outer cell of every draw. Evaluation order is kept.
    std::map<const void *, int> basis_ids;
    std::vector<int> outer_basis(outer.size(), -1);
    for (size_t j = 0; j < outer.size(); ++j)
        if ((outer[j].basis && !outer[j].basis->empty()) ||
            (outer[j].compact_basis && !outer[j].compact_basis->empty())) {
            const void *ptr = outer[j].compact_basis
                                  ? static_cast<const void *>(outer[j].compact_basis)
                                  : static_cast<const void *>(outer[j].basis);
            auto entry = basis_ids.emplace(ptr, static_cast<int>(basis_ids.size()));
            outer_basis[j] = entry.first->second;
        }
    std::vector<double> projected_bounds(basis_ids.size(), -1);
    std::map<const LocalLimit *, int> local_ids;
    std::vector<int> outer_local(outer.size(), -1);
    for (size_t j = 0; j < outer.size(); ++j)
        if (outer[j].local)
            outer_local[j] =
                local_ids.emplace(outer[j].local.get(), local_ids.size()).first->second;
    std::vector<I> local_dots(local_ids.size());
    std::vector<int> local_stamps(local_ids.size(), -1);
    std::vector<I> local_draw(s);
    // Adjacent continuous cells share corner vectors. Cache their EXACT interval
    // dot products lazily, with a fresh stamp for each draw. No geometric merging
    // or tolerance comparison is performed, and cone generators stay separate.
    std::map<std::vector<double>, int> ray_ids;
    std::vector<const std::vector<I> *> unique_rays;
    std::vector<std::array<int, 4>> cone_ray_ids(outer.size());
    for (size_t j = 0; j < outer.size(); ++j)
        if (outer[j].cone) {
            const auto &rays = outer[j].cone->rays;
            for (size_t k = 0; k < rays.size(); ++k) {
                if (prefix && !outer[j].cone->prefix_rays.empty()) {
                    cone_ray_ids[j][k] = -1;
                    continue;
                }
                std::vector<double> key;
                key.reserve(2 * s);
                for (I x : rays[k]) {
                    key.push_back(x.lo);
                    key.push_back(x.hi);
                }
                auto entry = ray_ids.emplace(std::move(key), static_cast<int>(unique_rays.size()));
                if (entry.second)
                    unique_rays.push_back(&rays[k]);
                cone_ray_ids[j][k] = entry.first->second;
            }
        }
    // Map keys are only needed to construct IDs, not in the Gaussian loop.
    ray_ids.clear();
    std::vector<I> ray_dots(unique_rays.size());
    std::vector<int> ray_stamp(unique_rays.size(), -1);
    NumericVector global_cap(B);
    for (int b = 0; b < B; ++b) {
        if ((b & 31) == 0)
            checkUserInterrupt();
        int id = draw_ids ? (*draw_ids)[b] - 1 : -1;
        const auto &draw = draw_ids ? prefix->gaussian[id] : local_draw;
        I raw_norm = pt(0);
        double cap;
        std::vector<std::array<I, 2>> local_prefix;
        if (draw_ids) {
            raw_norm = prefix->gaussian_norm[id];
            cap = prefix->gaussian_cap[id];
        } else {
            I projection = pt(0), norm2 = pt(0), raw2 = pt(0);
            for (int i = 0; i < s; ++i) {
                if (!R_finite(g(i, b))) stop("nonfinite Gaussian draw");
                local_draw[i] = pt(g(i, b));
                raw2 = plus(raw2, sq(local_draw[i]));
                if (!axis.empty()) projection = plus(projection, times_point(axis[i], g(i, b)));
            }
            for (int i = 0; i < s; ++i) {
                I z = axis.empty() ? draw[i] : minus(draw[i], times(axis[i], projection));
                norm2 = plus(norm2, sq(z));
            }
            cap = root(norm2).hi;
            raw_norm = root(raw2);
            if (prefix) local_prefix = prefix->draws(draw);
        }
        global_cap[b] = cap;
        const auto &prefix_draw = draw_ids ? prefix->gaussian_prefix[id] : local_prefix;
        // The same pair subspace is shared by many cells. Evaluate it once per draw.
        std::fill(projected_bounds.begin(), projected_bounds.end(), -1);
        auto witness = [&](size_t index) {
            const auto &d = unique_inner[index];
            I score = pt(0);
            if (prefix && !d.compact.empty())
                score = prefix->dot(d.compact, prefix_draw);
            else
                for (int i = 0; i < s; ++i)
                    score = plus(score, times_point(d.w[i], draw[i].lo));
            for (int h = 0; h < 2; ++h)
                if (d.mask & (1 << h))
                    lower(b, h) = std::max(lower(b, h), score.lo);
            ++witness_dots;
            return score;
        };
        if (!shared_lower || groups.empty()) {
            for (size_t j = 0; j < unique_inner.size(); ++j)
                witness(j);
        } else {
            auto seed_score = [&](size_t j) {
                if (seed_stamps[j] != b) {
                    seed_scores[j] = witness(groups[j].seed);
                    seed_stamps[j] = b;
                }
                return seed_scores[j];
            };
            // At most 24 initial comparisons, shared by G/L. The choice depends
            // only on the prepared bank, not on oracle truth or a different RNG.
            size_t initial = std::min(size_t(24), groups.size());
            for (size_t j = 0; j < initial; ++j)
                seed_score(j * groups.size() / initial);
            I draw_norm = raw_norm;
            for (size_t j = 0; j < groups.size(); ++j) {
                const auto &group = groups[j];
                // For every exact member v and seed v0, ||v-v0|| <= radius, hence
                // g'v <= seed_score.hi + ||g|| radius. All arithmetic rounds outward.
                double bound = plus(seed_score(j), times(draw_norm, pt(group.radius))).hi;
                bool improves = false;
                for (int h = 0; h < 2; ++h)
                    if ((group.mask & (1 << h)) && bound > lower(b, h))
                        improves = true;
                if (!improves) {
                    witness_skipped += group.last - group.first - 1;
                    continue;
                }
                for (size_t k = group.first; k < group.last; ++k)
                    if (k != group.seed)
                        witness(k);
            }
        }
        for (size_t outer_index = 0; outer_index < outer.size(); ++outer_index) {
            const auto &d = outer[outer_index];
            double bound = cap;
            int id = outer_basis[outer_index];
            if (id >= 0) {
                if (projected_bounds[id] < 0) {
                    I projected2 = pt(0);
                    bool uncertain_projection = false;
                    if (prefix && d.compact_basis)
                        for (const auto &q : *d.compact_basis) {
                            I score = prefix->dot(q, prefix_draw);
                            projected2 = plus(projected2, sq(score));
                            if (score.hi - score.lo >
                                1e-7 *
                                    std::max(1.0, std::max(std::abs(score.lo), std::abs(score.hi))))
                                uncertain_projection = true;
                        }
                    else
                        for (const auto &q : *d.basis) {
                            I score = pt(0);
                            for (int i = 0; i < s; ++i)
                                score = plus(score, times_point(q[i], draw[i].lo));
                            projected2 = plus(projected2, sq(score));
                        }
                    projected_bounds[id] = root(projected2).hi;
                    if (uncertain_projection && d.basis && !d.basis->empty()) {
                        I dense2 = pt(0);
                        for (const auto &q : *d.basis) {
                            I score = pt(0);
                            for (int i = 0; i < s; ++i)
                                score = plus(score, times_point(q[i], draw[i].lo));
                            dense2 = plus(dense2, sq(score));
                        }
                        projected_bounds[id] = std::min(projected_bounds[id], root(dense2).hi);
                    }
                }
                bound = std::min(bound, projected_bounds[id]);
            }
            double span_cap = bound;
            double profile_cap = d.profile && (d.mask & 2) ? d.profile->upper(*prefix, prefix_draw, bound) : bound;
            // Rank-directed work changes no upper certificate: skipped cells
            // retain this valid coarse bound; the complete cover is still evaluated.
            if (cutoff_targets && prefix) {
                bool relevant = false;
                for (int h = 0; h < 2; ++h)
                    if ((d.mask & (1 << h)) &&
                        (h == 1 ? profile_cap : bound) > (*cutoff_targets)(b, h)) relevant = true;
                if (!relevant) {
                    for (int h = 0; h < 2; ++h)
                        if (d.mask & (1 << h))
                            upper(b, h) = std::max(upper(b, h), h == 1 ? profile_cap : bound);
                    continue;
                }
            }
            // Remaining operations can only reduce bound. Skip if no selected maximum
            // can change; this retains the identical floating-point output, not a tolerance.
            bool improves = false;
            for (int h = 0; h < 2; ++h)
                if ((d.mask & (1 << h)) && bound > upper(b, h))
                    improves = true;
            bool improves_lower = false;
            // Do not solve another primal problem just to improve a few rounding
            // ulps. This only skips optional LOWER witnesses: the remaining width
            // is still reported, and no selection/sign/upper bound uses this rule.
            double primal_slack =
                64 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(bound));
            bool has_cone = d.cone || d.lazy_cone_model;
            bool radial_cone = d.cone ? static_cast<bool>(d.cone->radial) :
                d.lazy_cone_model && d.lazy_cone_model->radial_feasibility && d.lazy_cone_model->size >= 3;
            if (has_cone)
                for (int h = 0; h < 2; ++h)
                    if (((d.cone_certain | (radial_cone ? d.mask : 0)) & (1 << h)) &&
                        bound - lower(b, h) > primal_slack)
                        improves_lower = true;
            if (!improves && !improves_lower)
                continue;
            if (has_cone) {
                auto active_cone = d.cone ? d.cone : d.lazy_cone_model->cell_cone(d.lazy_cone_box);
                std::array<I, 4> h{};
                for (size_t k = 0; k < active_cone->rays.size(); ++k) {
                    if (d.lazy_cone_model) {
                        I dot = pt(0);
                        const auto &w = active_cone->rays[k];
                        for (int i = 0; i < s; ++i)
                            dot = plus(dot, times_point(w[i], draw[i].lo));
                        h[k] = dot;
                        continue;
                    }
                    int ray = cone_ray_ids[outer_index][k];
                    if (prefix && !active_cone->prefix_rays.empty()) {
                        h[k] = prefix->dot(active_cone->prefix_rays[k], prefix_draw);
                        continue;
                    }
                    if (ray_stamp[ray] != b) {
                        I dot = pt(0);
                        const auto &w = *unique_rays[ray];
                        for (int i = 0; i < s; ++i)
                            dot = plus(dot, times_point(w[i], draw[i].lo));
                        ray_dots[ray] = dot;
                        ray_stamp[ray] = b;
                    }
                    h[k] = ray_dots[ray];
                }
                std::array<double, 2> selected = {lower(b, 0), lower(b, 1)};
                I optimum = active_cone->optimize(draw, bound, h.data(), &selected,
                                             d.mask & ~d.cone_certain & 2);
                // Compressed sums can lose precision through cancellation. Only if an
                // unresolved cone could raise a reported upper bound, intersect its
                // dot products with the independent dense enclosures and retry.
                // This threshold selects extra work; it never certifies a sign/rank.
                bool relevant = false;
                for (int j = 0; j < 2; ++j)
                    if ((d.mask & (1 << j)) && optimum.hi > upper(b, j))
                        relevant = true;
                if (prefix && relevant && optimum.hi - optimum.lo > 1e-8 * std::max(1.0, bound)) {
                    active_cone->refine_geometry();
                    for (size_t k = 0; k < active_cone->rays.size(); ++k) {
                        int ray = cone_ray_ids[outer_index][k];
                        if (ray < 0) {
                            I dot = pt(0);
                            const auto &w = active_cone->rays[k];
                            for (int i = 0; i < s; ++i)
                                dot = plus(dot, times_point(w[i], draw[i].lo));
                            h[k] = {std::max(h[k].lo, dot.lo), std::min(h[k].hi, dot.hi)};
                            if (h[k].lo > h[k].hi)
                                stop("inconsistent lazy circle dot enclosure");
                            continue;
                        }
                        if (ray_stamp[ray] != b) {
                            I dot = pt(0);
                            const auto &w = *unique_rays[ray];
                            for (int i = 0; i < s; ++i)
                                dot = plus(dot, times_point(w[i], draw[i].lo));
                            ray_dots[ray] = dot;
                            ray_stamp[ray] = b;
                        }
                        h[k] = {std::max(h[k].lo, ray_dots[ray].lo),
                                std::min(h[k].hi, ray_dots[ray].hi)};
                        if (h[k].lo > h[k].hi)
                            stop("inconsistent circle dot enclosures");
                    }
                    I refined = active_cone->optimize(draw, bound, h.data());
                    optimum = {std::max(optimum.lo, refined.lo), std::min(optimum.hi, refined.hi)};
                }
                bound = std::min(bound, optimum.hi);
                for (int h = 1; h < 2; ++h)
                    lower(b, h) = std::max(lower(b, h), selected[h]);
                for (int h = 0; h < 2; ++h)
                    if (d.cone_certain & (1 << h))
                        lower(b, h) = std::max(lower(b, h), optimum.lo);
            }
            // Retain the independently certified polar fallback near shared endpoints.
            // It is used only after the cheap compressed cone has passed pruning.
            if (!d.rays.empty())
                bound = cone_bound(d.rays, draw, bound);
            if (!d.sphere) {
                I score = pt(0);
                if (prefix && !d.compact.empty())
                    score = prefix->dot(d.compact, prefix_draw);
                else
                    for (int i = 0; i < s; ++i)
                        score = plus(score, times_point(d.w[i], draw[i].lo));
                bound = std::min(bound, score.hi);
            }
            double local_bound = std::min(bound, profile_cap);
            if (d.local && (d.mask & 2) && bound > upper(b, 1)) {
                // For v=w/||w||, Local implies a'v >= -eta. For EVERY lambda>=0,
                // g'v <= ||P_span(g)+lambda*a|| + lambda*eta. Lambda only proposes
                // work; the resulting bound is independently rounded outward.
                int local = outer_local[outer_index];
                if (local_stamps[local] != b) {
                    I dot = pt(0);
                    for (int i = 0; i < s; ++i)
                        dot = plus(dot, times_point(d.local->normal[i], draw[i].lo));
                    local_dots[local] = dot;
                    local_stamps[local] = b;
                }
                I dot = local_dots[local];
                double aa = d.local->norm2.lo / 2 + d.local->norm2.hi / 2;
                double h = dot.lo / 2 + dot.hi / 2, eta = d.local_eta, lambda = 0;
                if (aa > eta * eta && h < -eta * span_cap)
                    lambda = std::max(
                        0.0,
                        (-h - eta * std::sqrt(std::max(0.0, (aa * span_cap * span_cap - h * h) /
                                                                (aa - eta * eta)))) /
                            aa);
                if (lambda > 0 && R_finite(2 * lambda)) {
                    I norm2 = plus(plus(sq(pt(span_cap)), times_point(dot, 2 * lambda)),
                                   times(sq(pt(lambda)), d.local->norm2));
                    if (norm2.hi >= 0)
                        local_bound =
                            std::min(local_bound, plus(root(norm2), times(pt(lambda), pt(eta))).hi);
                }
            }
            for (int h = 0; h < 2; ++h)
                if (d.mask & (1 << h))
                    upper(b, h) = std::max(upper(b, h), h ? local_bound : bound);
        }
        for (int h = 0; h < 2; ++h)
            if (!lower_only && (covered_methods & (1 << h)) && lower(b, h) > upper(b, h))
                stop("pathwise bounds inconsistent: draw=%d method=%d lower=%.17g upper=%.17g", b,
                     h, lower(b, h), upper(b, h));
    }
    IntegerVector active(2);
    for (const auto &d : outer)
        for (int h = 0; h < 2; ++h)
            if (d.mask & (1 << h))
                ++active[h];
    return List::create(_["lower"] = lower, _["upper"] = upper, _["global_cap"] = global_cap,
                        _["representatives"] = inner.size(), _["outer_directions"] = outer.size(),
                        _["active_methods"] = active, _["witness_dot_products"] = witness_dots,
                        _["witness_dots_skipped"] = witness_skipped);
}
} // namespace cone

namespace cone {
// Blocks share the immutable observed model and its lazy coefficient caches.
// Shared ownership keeps direction pointers valid after releasing the R handle.
struct Geometry {
    std::shared_ptr<Model> owner;
    Model &model;
    std::vector<Direction> inner, outer, global_outer;
    std::vector<WitnessGroup> lower_groups;
    List metadata;
    bool witness_only = false;
    bool global_only_G = false;
    explicit Geometry(std::shared_ptr<Model> source)
        : owner(std::move(source)), model(*owner), metadata(List::create()) {}
};
} // namespace cone

SEXP build_piecewise_model_cpp(List pieces, NumericVector counts, double eps, double tau,
                               bool circle_prefix = true, bool radial_feasibility = true,
                               bool sharp_rectangles = false, Nullable<NumericMatrix> draws = R_NilValue) {
    fc::check_platform();
    if (!R_finite(eps) || eps < 0 || !R_finite(tau) || tau <= 0)
        stop("invalid controls");
    auto m = std::make_shared<cone::Model>(cone::model(pieces, counts, eps, tau));
    XPtr<std::shared_ptr<cone::Model>> out(new std::shared_ptr<cone::Model>(m), true);
    m->radial_feasibility = radial_feasibility;
    m->sharp_rectangles = sharp_rectangles;
    if (circle_prefix && pieces.hasAttribute("circle_observations")) {
        NumericVector y = pieces.attr("circle_observations");
        if (y.size() != counts.size() || pieces.size() != counts.size() + 1)
            stop("invalid circle prefix metadata");
        Function identical = Environment::base_env()["identical"];
        if (!as<bool>(identical(pieces, circle_charts_cpp(y))))
            stop("circle prefix metadata does not match coefficients");
        m->prefix = std::make_shared<cone::CirclePrefix>(y, m->p, m->axis);
    }
    if (draws.isNotNull() && m->prefix) m->prefix->prepare_draws(as<NumericMatrix>(draws));
    return out;
}
void release_piecewise_model_cpp(SEXP model_ptr) {
    XPtr<std::shared_ptr<cone::Model>> x(model_ptr);
    x.release();
}

SEXP build_piecewise_cover_cpp(SEXP model_ptr, List pieces, NumericVector counts, double eps,
                               double tau, int max_cells = 10000, int max_depth = 8,
                               double direction_tolerance = .1, int pair_start = 0,
                               int pair_count = -1, bool continuous_cones = false,
                               int requested_methods = 3, bool witness_only = false,
                               bool independent_G_cover = true) {
    fc::check_platform();
    if (!R_finite(eps) || eps < 0 || !R_finite(tau) || tau <= 0 ||
        max_cells < 1 || max_depth < 0 || max_depth > 30 || !R_finite(direction_tolerance) ||
        direction_tolerance <= 0)
        stop("invalid controls");
    if (requested_methods < 1 || requested_methods > 3)
        stop("invalid requested methods");
    // In G-only root mode the retained leaf cover IS the independent root
    // cover. Avoid computing the identical full-domain cones twice. Existing
    // callers keep their old independent-cover path by default.
    if (!independent_G_cover && (max_depth != 0 || requested_methods != 1 || witness_only))
        stop("shared G root cover requires depth zero, G only, and upper geometry");
    // At depth zero with G alone, the leaf and independent root covers have
    // exactly the same boxes, cones, feasible witnesses and certain G mask.
    // Evaluate that complete cover once. Local subdivisions still retain their
    // independent roots because their small cones can be ill conditioned.
    if (max_depth == 0 && requested_methods == 1 && !witness_only && continuous_cones)
        independent_G_cover = false;
    double total_pairs = static_cast<double>(pieces.size()) * pieces.size();
    if (total_pairs > 2147483647.0)
        stop("pair index limit exceeded");
    if (pair_count < 0)
        pair_count = static_cast<int>(total_pairs);
    if (pair_start < 0 || pair_count < 1 ||
        pair_start + static_cast<double>(pair_count) > total_pairs)
        stop("invalid pair block");
    if (pair_count > max_cells)
        stop("initial full-domain cell count exceeds budget");
    XPtr<std::shared_ptr<cone::Model>> prepared(model_ptr);
    XPtr<cone::Geometry> geometry(new cone::Geometry(*prepared), true);
    geometry->witness_only = witness_only;
    auto &m = geometry->model;
    geometry->global_only_G = m.prefix && independent_G_cover && continuous_cones && (requested_methods & 1);
    if (witness_only && !m.prefix)
        stop("witness-only geometry requires the circle prefix model");
    auto &inner = geometry->inner;
    auto &outer = geometry->outer;
    struct Work {
        cone::Box box;
        bool adaptive;
        int serial;
        cone::Evaluation evaluation;
    };
    struct Priority {
        bool operator()(const Work &a, const Work &b) const {
            if (a.adaptive != b.adaptive)
                return a.adaptive < b.adaptive;
            if (a.box.depth != b.box.depth)
                return a.box.depth > b.box.depth;
            return a.serial > b.serial;
        }
    };
    std::priority_queue<Work, std::vector<Work>, Priority> queue;
    int serial = 0;
    std::vector<cone::Box> retained;
    auto enqueue = [&](cone::Box box) {
        auto eval = m.evaluate(box);
        eval.possible &= requested_methods;
        eval.certain &= requested_methods;
        for (auto &d : eval.directions)
            d.mask &= requested_methods;
        if (!eval.possible) {
            if (m.prefix) retained.push_back(box);
            return;
        }
        // An independent root cover for G avoids losing a well-conditioned cone
        // when Local boundary refinement produces almost dependent tiny generators.
        // Both covers include the same full pair domain; their global bounds may
        // be intersected. Never use an unselected root cone as a Local lower bound.
        if (independent_G_cover && !witness_only && continuous_cones && (requested_methods & 1) &&
            box.depth == 0 && eval.possible) {
            auto cone = box.j != box.k ? m.cell_cone(box) : nullptr;
            for (auto d : eval.directions) {
                d.mask = 1;
                d.cone = cone;
                d.cone_certain = 1;
                geometry->global_outer.push_back(std::move(d));
            }
        }
        queue.push({box, static_cast<bool>(eval.possible & 2), serial++, std::move(eval)});
    };
    auto old = m.cover_plans.find(pair_start);
    bool reuse = m.prefix && old != m.cover_plans.end() && old->second.count == pair_count &&
                 old->second.constructed <= max_cells;
    if (reuse)
        for (const auto &box : old->second.leaves) if (box.depth > max_depth) reuse = false;
    // Root G cover remains independent of the retained Local subdivision.
    if (reuse && independent_G_cover && continuous_cones && (requested_methods & 1))
        for (int index = pair_start; index < pair_start + pair_count; ++index) {
            cone::Box box{index / static_cast<int>(pieces.size()), index % static_cast<int>(pieces.size()), 0, 0, 1, 0, 1};
            auto eval = m.evaluate(box);
            auto cone = box.j != box.k ? m.cell_cone(box) : nullptr;
            for (auto d : eval.directions) {
                d.mask = 1; d.cone = cone; d.cone_certain = 1;
                geometry->global_outer.push_back(std::move(d));
            }
        }
    if (reuse) {
        for (auto box : old->second.leaves) enqueue(box);
    } else
        for (int index = pair_start; index < pair_start + pair_count; ++index)
            enqueue({index / static_cast<int>(pieces.size()), index % static_cast<int>(pieces.size()), 0, 0, 1, 0, 1});
    int constructed = reuse ? old->second.constructed : pair_count;
    int reused_cells = reuse ? old->second.leaves.size() : 0;

    bool rational = true;
    for (const auto &a : m.pieces) {
        rational = rational && a.left.lo == a.left.hi && a.right.lo == a.right.hi;
        for (int i = 0; i < m.size; ++i)
            rational = rational && a.b[i].lo == a.b[i].hi && a.c[i].lo == a.c[i].hi;
    }
    std::vector<std::array<double, 4>> pending_points;
    std::vector<std::vector<cone::Direction>> pending_directions;
    std::map<std::array<double, 4>, cone::Evaluation> point_cache;
    // Graph representatives are evaluated once, in the original visit order.
    // A repeated point needs only its certain selector mask: its directions
    // were already copied to inner/pending on the first visit. Retaining that
    // dense Evaluation again grows with every dyadic grid point unnecessarily.
    std::map<std::array<double, 4>, int> graph_point_masks;
    int visited = 0, leaves = 0, budget_leaves = 0;
    IntegerVector possible_counts(2), certain_counts(2);
    while (!queue.empty()) {
        if ((visited & 127) == 0)
            checkUserInterrupt();
        cone::Box box = queue.top().box;
        auto eval = queue.top().evaluation;
        queue.pop();
        ++visited;
        if (!eval.possible) {
            if (m.prefix) retained.push_back(box);
            continue;
        }
        int representative_mask = 0;
        double um = (box.u0 + box.u1) / 2, vm = (box.v0 + box.v1) / 2;
        for (auto uv : std::vector<std::pair<double, double>>{{um, vm},
                                                              {box.u0, box.v0},
                                                              {box.u0, box.v1},
                                                              {box.u1, box.v0},
                                                              {box.u1, box.v1}}) {
            std::array<double, 4> key = {static_cast<double>(box.j), static_cast<double>(box.k),
                                         uv.first, uv.second};
            if (!m.prefix) {
                auto found = graph_point_masks.find(key);
                if (found != graph_point_masks.end()) {
                    representative_mask |= found->second & requested_methods;
                    continue;
                }
                auto at = m.evaluate({box.j, box.k, box.depth, uv.first, uv.first,
                                      uv.second, uv.second}, false);
                graph_point_masks.emplace(key, at.certain);
                representative_mask |= at.certain & requested_methods;
                for (auto d : at.directions) {
                    d.mask = at.certain & requested_methods;
                    if (d.mask && !d.sphere)
                        inner.push_back(d);
                }
                if (rational && (at.possible & requested_methods) !=
                    (at.certain & requested_methods) && !at.directions.empty()) {
                    pending_points.push_back(
                        {static_cast<double>(box.j), static_cast<double>(box.k), uv.first, uv.second});
                    pending_directions.push_back(at.directions);
                }
                continue;
            }
            auto found = point_cache.find(key);
            bool fresh = found == point_cache.end();
            if (fresh)
                found = point_cache
                            .emplace(key, m.evaluate({box.j, box.k, box.depth, uv.first, uv.first,
                                                      uv.second, uv.second},
                                                     false))
                            .first;
            const auto &at = found->second;
            representative_mask |= at.certain & requested_methods;
            if (fresh)
                for (auto d : at.directions) {
                    d.mask = at.certain & requested_methods;
                    if (d.mask && !d.sphere)
                        inner.push_back(d);
                }
            if (fresh && rational &&
                (at.possible & requested_methods) != (at.certain & requested_methods) &&
                !at.directions.empty()) {
                pending_points.push_back(
                    {static_cast<double>(box.j), static_cast<double>(box.k), uv.first, uv.second});
                pending_directions.push_back(at.directions);
            }
        }
        int adaptive_mask = geometry->global_only_G ? (eval.possible & 2) : eval.possible;
        bool needs = adaptive_mask && (eval.width > direction_tolerance || (adaptive_mask & ~representative_mask));
        int children = 4;
        if (needs && box.depth < max_depth &&
            (m.prefix ? constructed + children <= max_cells :
             visited + queue.size() + children <= static_cast<size_t>(max_cells))) {
            constructed += children;
            for (auto u : std::vector<std::pair<double, double>>{{box.u0, um}, {um, box.u1}})
                for (auto w : std::vector<std::pair<double, double>>{{box.v0, vm}, {vm, box.v1}})
                    enqueue({box.j, box.k, box.depth + 1, u.first, u.second, w.first, w.second});
        } else {
            if (m.prefix) retained.push_back(box);
            ++leaves;
            if (needs)
                ++budget_leaves;
            // Screening needs only certified feasible point witnesses.  Retain the
            // exact same subdivision, representative evaluations, and rational
            // feasibility checks, but do not build any upper-cover cones or leaves.
            if (!witness_only) {
                int leaf_mask = geometry->global_only_G ? (eval.possible & 2) : eval.possible;
                bool lazy_graph_cone = !m.prefix && max_cells > 8192 && leaf_mask &&
                    continuous_cones && box.j != box.k;
                auto cone = !lazy_graph_cone && leaf_mask && continuous_cones && box.j != box.k ? m.cell_cone(box) : nullptr;
                auto profile = m.prefix && box.j != box.k && (leaf_mask & 2) ? m.circle_profile(box) : nullptr;
                for (auto d : eval.directions) {
                    d.mask &= leaf_mask;
                    if (!d.mask) continue;
                    d.cone = cone;
                    if (lazy_graph_cone) { d.lazy_cone_model = &m; d.lazy_cone_box = box; }
                    d.profile = profile;
                    d.cone_certain = (1 | eval.certain) & requested_methods;
                    outer.push_back(d);
                }
            }
            for (int h = 0; h < 2; ++h) {
                if (eval.possible & (1 << h))
                    ++possible_counts[h];
                if (representative_mask & (1 << h))
                    ++certain_counts[h];
            }
        }
    }
    if (!m.prefix) {
        // Neither cache is read again. Release graph point masks and the empty
        // priority queue's retained capacity before the exact-key dedup peak.
        std::map<std::array<double, 4>, int>().swap(graph_point_masks);
        decltype(queue)().swap(queue);
    }
    // Exact rational decisions only strengthen feasible representative lower bounds.
    // Geometry and all upper masks remain conservative and independent of draws.
    if (!pending_points.empty()) {
        NumericMatrix points(pending_points.size(), 4);
        for (size_t j = 0; j < pending_points.size(); ++j)
            for (int k = 0; k < 4; ++k)
                points(j, k) = pending_points[j][k];
        Function resolve = Environment::global_env()["exact_selector_masks"];
        IntegerVector masks = resolve(pieces, counts, points, eps, tau);
        if (masks.size() != static_cast<int>(pending_points.size()))
            stop("exact mask dimensions");
        for (size_t j = 0; j < pending_points.size(); ++j)
            for (auto &d : pending_directions[j]) {
                if ((masks[j] & requested_methods) && !d.sphere) {
                    d.mask = masks[j] & requested_methods;
                    inner.push_back(std::move(d));
                }
            }
    }
    if (m.prefix) {
        if (old != m.cover_plans.end()) {
            m.plan_bytes -= old->second.leaves.size() * sizeof(cone::Box);
            m.cover_plans.erase(old);
        }
        m.plan_bytes += retained.size() * sizeof(cone::Box);
        m.cover_plans.emplace(pair_start, cone::CoverPlan{pair_start, pair_count, constructed, std::move(retained)});
        while (m.plan_bytes > 64 * 1024 * 1024 && !m.cover_plans.empty()) {
            auto first = m.cover_plans.begin();
            m.plan_bytes -= first->second.leaves.size() * sizeof(cone::Box);
            m.cover_plans.erase(first);
        }
    }
    List out = List::create();
    out["exact_selection_points"] = pending_points.size();
    out["pair_start"] = pair_start;
    out["pair_count"] = pair_count;
    out["total_pairs"] = total_pairs;
    out["visited_cells"] = visited;
    out["reused_cells"] = reused_cells;
    out["leaves"] = leaves;
    out["budget_leaves"] = budget_leaves;
    out["possible_cells"] = possible_counts;
    out["represented_cells"] = certain_counts;
    out["fmin"] = NumericVector::create(m.fmin.lo, m.fmin.hi);
    out["representatives"] = inner.size(); // Preserve the construction count.
    if (!m.prefix)
        inner = cone::deduplicate_owned_inner(std::move(inner), m.size);
    geometry->lower_groups = cone::witness_groups(inner, m.prefix.get());
    geometry->metadata = out;
    return geometry;
}

List evaluate_piecewise_draws_cpp(SEXP geometry_ptr, NumericMatrix g, bool shared_lower = true,
                                  Nullable<List> prior = R_NilValue,
                                  Nullable<IntegerVector> draw_ids = R_NilValue,
                                  Nullable<NumericMatrix> cutoff_targets = R_NilValue) {
    XPtr<cone::Geometry> x(geometry_ptr);
    if (x->witness_only)
        stop("witness-only geometry cannot supply a full-domain upper bound");
    IntegerVector ids;
    const IntegerVector *prepared_ids = nullptr;
    if (draw_ids.isNotNull()) {
        ids = as<IntegerVector>(draw_ids);
        if (!x->model.prefix || !ids.size()) stop("missing prepared circle draws");
        for (int id : ids)
            if (id < 1 || id > static_cast<int>(x->model.prefix->gaussian.size()))
                stop("invalid prepared draw index");
        prepared_ids = &ids;
    } else if (g.nrow() != x->model.size || g.ncol() < 1)
        stop("Gaussian dimensions");
    NumericMatrix pl, pu;
    NumericMatrix targets;
    const NumericMatrix *target_ptr = nullptr;
    if (cutoff_targets.isNotNull()) {
        targets = as<NumericMatrix>(cutoff_targets);
        if (targets.nrow() != (prepared_ids ? ids.size() : g.ncol()) || targets.ncol() != 2 ||
            is_true(any(is_na(targets)))) stop("invalid cutoff targets");
        target_ptr = &targets;
    }
    const NumericMatrix *lower = nullptr;
    const NumericMatrix *upper = nullptr;
    if (prior.isNotNull()) {
        List p(prior);
        pl = as<NumericMatrix>(p["lower"]);
        pu = as<NumericMatrix>(p["upper"]);
        lower = &pl;
        upper = &pu;
    }
    List out = cone::draw_bounds(x->inner, x->outer, g, x->model.axis, true, x->model.prefix.get(),
                                 &x->lower_groups, shared_lower, false, lower, upper, prepared_ids, target_ptr,
                                 x->global_only_G ? 2 : 3);
    if (!x->global_outer.empty()) {
        List root = cone::draw_bounds({}, x->global_outer, g, x->model.axis, true,
                                      x->model.prefix.get(), nullptr, true, false, lower, upper, prepared_ids, target_ptr);
        NumericMatrix lo = out["lower"], hi = out["upper"], rlo = root["lower"],
                      rhi = root["upper"];
        for (int b = 0; b < lo.nrow(); ++b) {
            lo(b, 0) = std::max(lo(b, 0), rlo(b, 0));
            hi(b, 0) = x->global_only_G ? rhi(b, 0) : std::min(hi(b, 0), rhi(b, 0));
            if (lo(b, 0) > hi(b, 0))
                stop("inconsistent independent G covers");
        }
    }
    CharacterVector keys = x->metadata.names();
    for (int j = 0; j < keys.size(); ++j)
        out[as<std::string>(keys[j])] = x->metadata[j];
    return out;
}
void release_piecewise_cover_cpp(SEXP geometry_ptr) {
    XPtr<cone::Geometry> x(geometry_ptr);
    x.release();
}
// Rank-threshold branch and bound; unvisited roots keep the previous full-domain
// upper, and every discarded cell retains its certified upper. No truth is used.
List refine_rank_draws_cpp(SEXP model_ptr, NumericMatrix g, NumericMatrix lower,
                           NumericMatrix upper, IntegerVector ids, double target,
                           int max_cells = 262144, int max_depth = 60) {
    XPtr<std::shared_ptr<cone::Model>> ptr(model_ptr);
    auto &m = **ptr;
    if (g.nrow() != m.size || lower.nrow() != g.ncol() || upper.nrow() != g.ncol() ||
        lower.ncol() != 2 || upper.ncol() != 2 || !R_finite(target) || target < 0 ||
        max_cells < static_cast<int>(m.pieces.size()*m.pieces.size()) || max_depth < 0 || max_depth > 60)
        stop("invalid rank refinement controls");
    lower = clone(lower); upper = clone(upper);
    for (int b = 0; b < g.ncol(); ++b) for (int j = 0; j < 2; ++j)
        if (!R_finite(lower(b, j)) || !R_finite(upper(b, j)) || lower(b, j) < 0 || lower(b, j) > upper(b, j))
            stop("invalid previous draw enclosure");
    struct Node { cone::Box box; double upper; };
    struct Priority { bool operator()(const Node &a, const Node &b) const { return a.upper < b.upper; } };
    IntegerVector visited(ids.size());
    NumericMatrix peaks(ids.size(), 7);
    for (int index = 0; index < ids.size(); ++index) {
        int b = ids[index]-1;
        if (b < 0 || b >= g.ncol()) stop("invalid draw index");
        double incumbent = lower(b, 1), discarded = incumbent;
        if (incumbent > target || upper(b, 1) <= target) continue;
        NumericMatrix draw(m.size, 1);
        for (int i = 0; i < m.size; ++i) draw(i, 0) = g(i, b);
        IntegerVector prepared = IntegerVector::create(b+1);
        const IntegerVector *draw_ids = m.prefix ? &prepared : nullptr;
        std::priority_queue<Node, std::vector<Node>, Priority> queue;
        Node terminal{{0, 0, 0, 0, 0, 0, 0}, 0};
        auto evaluate = [&](cone::Box box, double parent_upper) {
            if ((visited[index]++ & 127) == 0) checkUserInterrupt();
            cone::Evaluation uncached;
            bool root=m.prefix && m.pieces.size()<=201 && box.depth==0 &&
                box.u0==0 && box.u1==1 && box.v0==0 && box.v1==1;
            const auto &at=root ? m.rank_root_evaluation(box) : (uncached=m.evaluate(box),uncached);
            if (!(at.possible & 2)) return;
            std::vector<cone::Direction> inner, outer;
            double um = (box.u0+box.u1)/2, vm = (box.v0+box.v1)/2;
            for (auto uv : std::vector<std::pair<double, double>>{{um, vm}, {box.u0, box.v0},
                         {box.u0, box.v1}, {box.u1, box.v0}, {box.u1, box.v1}}) {
                auto point = m.evaluate({box.j, box.k, box.depth, uv.first, uv.first, uv.second, uv.second}, false);
                for (auto d : point.directions) if ((point.certain & 2) && !d.sphere) {
                    d.mask = 2; inner.push_back(std::move(d));
                }
            }
            auto cell = box.j != box.k ? m.cell_cone(box) : nullptr;
            for (auto d : at.directions) if (d.mask & 2) {
                d.mask = 2; d.cone = cell; d.cone_certain = at.certain & 2;
                outer.push_back(std::move(d));
            }
            NumericMatrix prior(1, 2); prior(0, 1) = incumbent;
            auto bound = cone::draw_bounds(inner, outer, draw, m.axis, true, m.prefix.get(),
                nullptr, false, false, &prior, &prior, draw_ids, nullptr, 2);
            NumericMatrix lo = bound["lower"], hi = bound["upper"];
            incumbent = std::max(incumbent, lo(0, 1));
            double cap = std::min(parent_upper, hi(0, 1));
            if (cap > target && m.prefix && box.j != box.k) {
                auto profile = m.circle_profile(box);
                cap = std::min(cap, profile->upper(*m.prefix, m.prefix->gaussian_prefix[b], cap));
                if (cap > target) {
                    auto shared_profile = m.circle_shared_profile(box);
                    if (shared_profile)
                        cap = std::min(cap,shared_profile->upper(*m.prefix,m.prefix->gaussian_prefix[b],cap));
                }
            }
            if (incumbent > upper(b, 1)) stop("rank witness exceeds previous certified upper");
            cap = std::max(cap, incumbent);
            if (cap <= target) discarded = std::max(discarded, cap);
            else queue.push({box, cap});
        };
        bool complete = true;
        for (int j = 0; j < static_cast<int>(m.pieces.size()) && complete; ++j)
            for (int k = 0; k < static_cast<int>(m.pieces.size()); ++k) {
                evaluate({j, k, 0, 0, 1, 0, 1}, upper(b, 1));
                if (incumbent > target) { complete = false; break; }
            }
        if (complete) {
            while (!queue.empty() && incumbent <= target && visited[index]+2 <= max_cells) {
                Node top = queue.top(); queue.pop();
                if (top.upper <= target || top.box.depth == max_depth) {
                    discarded = std::max(discarded, top.upper);
                    if (top.upper > terminal.upper) terminal = top;
                    continue;
                }
                const auto &box = top.box;
                double um = (box.u0+box.u1)/2, vm = (box.v0+box.v1)/2;
                const auto &a = m.pieces[box.j], &c = m.pieces[box.k];
                double qwidth = (a.right.hi-a.left.lo)*(box.u1-box.u0)*std::sqrt(a.vb.hi);
                double rwidth = (c.right.hi-c.left.lo)*(box.v1-box.v0)*std::sqrt(c.vb.hi);
                if ((qwidth >= rwidth && um > box.u0 && um < box.u1) || vm == box.v0 || vm == box.v1) {
                    if (um == box.u0 || um == box.u1) { discarded = std::max(discarded, top.upper); continue; }
                    evaluate({box.j, box.k, box.depth+1, box.u0, um, box.v0, box.v1}, top.upper);
                    evaluate({box.j, box.k, box.depth+1, um, box.u1, box.v0, box.v1}, top.upper);
                } else {
                    evaluate({box.j, box.k, box.depth+1, box.u0, box.u1, box.v0, vm}, top.upper);
                    evaluate({box.j, box.k, box.depth+1, box.u0, box.u1, vm, box.v1}, top.upper);
                }
            }
            double cap = std::max(discarded, incumbent);
            if (!queue.empty()) cap = std::max(cap, queue.top().upper);
            upper(b, 1) = std::min(upper(b, 1), cap);
        }
        lower(b, 1) = incumbent;
        if (!queue.empty() || terminal.upper > target) {
            const auto &box = queue.empty() || terminal.upper > queue.top().upper ? terminal.box : queue.top().box;
            double values[] = {double(box.j+1), double(box.k+1), double(box.depth), box.u0, box.u1, box.v0, box.v1};
            for (int j = 0; j < 7; ++j) peaks(index, j) = values[j];
        }
    }
    return List::create(_["lower"] = lower, _["upper"] = upper, _["visited_cells"] = visited, _["remaining_peaks"] = peaks);
}

// Re-evaluate a point bank with the SAME enclosing model as the full cover.
// Coordinates are (piece_q, u, piece_r, v), one-based IDs, u/v in [0,1].
// A point optimizer's raw objective is never admitted as a certified lower.

List evaluate_witnesses_cpp(SEXP model_ptr, List pieces, NumericVector counts, NumericMatrix points,
                            NumericMatrix g, double eps, double tau) {
    XPtr<std::shared_ptr<cone::Model>> ptr(model_ptr);
    auto &m = **ptr;
    if (points.ncol() != 4 || g.nrow() != m.size || g.ncol() < 1)
        stop("invalid witness dimensions");
    std::vector<cone::Direction> inner;
    std::vector<std::array<double, 4>> pending;
    std::vector<std::vector<cone::Direction>> uncertain;
    bool rational = !m.prefix;
    for (int h = 0; h < points.nrow(); ++h) {
        for (int k = 0; k < 4; ++k)
            if (!R_finite(points(h, k)))
                stop("nonfinite witness");
        int j = points(h, 0) - 1, k = points(h, 2) - 1;
        double u = points(h, 1), w = points(h, 3);
        if (points(h, 0) != j + 1 || points(h, 2) != k + 1 || j < 0 || k < 0 ||
            j >= static_cast<int>(m.pieces.size()) || k >= static_cast<int>(m.pieces.size()) ||
            u < 0 || u > 1 || w < 0 || w > 1)
            stop("invalid witness coordinates");
        auto at = m.evaluate({j, k, 0, u, u, w, w}, false);
        for (auto d : at.directions)
            if (!d.sphere && at.certain) {
                d.mask = at.certain;
                inner.push_back(d);
            }
        if (rational && at.possible != at.certain && !at.directions.empty()) {
            pending.push_back({static_cast<double>(j), static_cast<double>(k), u, w});
            uncertain.push_back(at.directions);
        }
    }
    if (!pending.empty()) {
        NumericMatrix p(pending.size(), 4);
        for (size_t j = 0; j < pending.size(); ++j)
            for (int k = 0; k < 4; ++k)
                p(j, k) = pending[j][k];
        Function resolve = Environment::global_env()["exact_selector_masks"];
        IntegerVector masks = resolve(pieces, counts, p, eps, tau);
        for (size_t j = 0; j < pending.size(); ++j)
            for (auto d : uncertain[j])
                if (!d.sphere && masks[j]) {
                    d.mask = masks[j];
                    inner.push_back(d);
                }
    }
    auto out = cone::draw_bounds(inner, {}, g, m.axis, false, m.prefix.get(), nullptr, true, true);
    return List::create(_["lower"] = out["lower"], _["certified_points"] = inner.size(),
                        _["exact_selection_points"] = pending.size(), _["global_cap"] = out["global_cap"],
                        _["fmin"] = NumericVector::create(m.fmin.lo, m.fmin.hi));
}

// Build circle coefficients with an enclosure of mathematical pi, not rounded pi.

List circle_charts_cpp(NumericVector y) {
    fc::check_platform();
    using namespace cone;
    I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1};
    I twopi = times(pt(2), pi);
    int s = y.size();
    if (s < 1)
        stop("empty circle data");
    struct Cut {
        I value;
        int observation;
    };
    std::vector<Cut> cuts;
    std::vector<I> lift;
    for (int i = 0; i < s; ++i) {
        if (!R_finite(y[i]) || y[i] < -pi.lo || y[i] > pi.lo)
            stop("circle input must be in [-pi_double,pi_double]");
        lift.push_back(y[i] >= 0 ? minus(pt(y[i]), twopi) : pt(y[i]));
        cuts.push_back({y[i] >= 0 ? minus(pt(y[i]), pi) : plus(pt(y[i]), pi), i});
    }
    std::sort(cuts.begin(), cuts.end(), [&](const Cut &a, const Cut &b) {
        bool ap = y[a.observation] >= 0, bp = y[b.observation] >= 0;
        return ap != bp ? ap : y[a.observation] < y[b.observation];
    });
    List pieces;
    I left = neg(pi);
    for (int j = 0; j <= s; ++j) {
        I right = j < s ? cuts[j].value : pi;
        NumericMatrix b(s, 2), c(s, 2);
        for (int i = 0; i < s; ++i) {
            I bb = times(pt(-2), lift[i]), cc = sq(lift[i]);
            b(i, 0) = bb.lo;
            b(i, 1) = bb.hi;
            c(i, 0) = cc.lo;
            c(i, 1) = cc.hi;
        }
        pieces.push_back(List::create(_["left"] = NumericVector::create(left.lo, left.hi),
                                      _["right"] = NumericVector::create(right.lo, right.hi),
                                      _["b"] = b, _["c"] = c, _["left_id"] = j,
                                      _["right_id"] = (j == s ? 0 : j + 1), _["circle_index"] = j,
                                      _["circle_support"] = s));
        if (j < s) {
            lift[cuts[j].observation] = plus(lift[cuts[j].observation], twopi);
            left = right;
        }
    }
    pieces.attr("circle_observations") = clone(y);
    return pieces;
}

} // namespace piecewise_kernel

#include "circle/numerics.hpp"
#include "circle/length.hpp"
#include "euclid/numerics.hpp"
#include "ranking/numerics.hpp"

// [[Rcpp::export]]
double cutoff_index_cpp(int B, double alpha) { return finite_kernel::cutoff_index_cpp(B, alpha); }
// [[Rcpp::export]]
NumericVector interval_width_cpp(NumericVector lower, NumericVector upper) {
    return finite_kernel::interval_width_cpp(lower, upper);
}
// [[Rcpp::export]]
List euclid_projection_cpp(NumericMatrix y, NumericMatrix g, IntegerVector columns,
                           Nullable<NumericMatrix> basis_lower = R_NilValue,
                           Nullable<NumericMatrix> basis_upper = R_NilValue) {
    return finite_kernel::euclid_projection_cpp(y, g, columns, basis_lower, basis_upper);
}
// [[Rcpp::export]]
IntegerVector ranking_selector_masks_cpp(NumericVector risk, NumericVector sums,
                                         NumericVector variances, IntegerMatrix pairs, double n,
                                         double eps, double tau) {
    return finite_kernel::ranking_selector_masks_cpp(risk, sums, variances, pairs, n, eps, tau);
}
// [[Rcpp::export]]
List ranking_maxima_cpp(NumericMatrix loss, NumericVector counts, NumericMatrix g,
                        IntegerMatrix pairs, IntegerVector masks, LogicalVector positive,
                        NumericVector sums, NumericVector variances) {
    return finite_kernel::ranking_maxima_cpp(loss, counts, g, pairs, masks, positive, sums,
                                             variances);
}
// [[Rcpp::export]]
NumericMatrix ranking_distances_cpp(IntegerMatrix ranks) {
    return finite_kernel::ranking_distances_cpp(ranks);
}
// [[Rcpp::export]]
void merge_piecewise_bounds_cpp(NumericMatrix lower, NumericMatrix upper, NumericMatrix chunk_lower,
                                NumericMatrix chunk_upper, int offset) {
    piecewise_kernel::merge_piecewise_bounds_cpp(lower, upper, chunk_lower, chunk_upper, offset);
}
// [[Rcpp::export]]
List piecewise_state_cpp(SEXP model_ptr, Nullable<List> restore = R_NilValue) {
    XPtr<std::shared_ptr<piecewise_kernel::cone::Model>> handle(model_ptr);
    auto &m = **handle;
    using piecewise_kernel::cone::Box;
    using piecewise_kernel::cone::CoverPlan;
    if (restore.isNotNull()) {
        if (!m.prefix) stop("subdivision state requires circle geometry");
        List state(restore);
        std::map<int, CoverPlan> plans;
        size_t bytes = 0;
        for (SEXP entry : state) {
            List record(entry);
            int start = as<int>(record["start"]), count = as<int>(record["count"]);
            NumericMatrix cells = record["cells"];
            if (start < 0 || count < 1 || start + static_cast<double>(count) >
                static_cast<double>(m.pieces.size()) * m.pieces.size() || cells.ncol() != 6)
                stop("invalid subdivision state");
            CoverPlan plan{start, count, 0, {}};
            std::vector<std::vector<std::pair<uint64_t, uint64_t>>> partitions(count);
            for (int row = 0; row < cells.nrow(); ++row) {
                for (int col = 0; col < 6; ++col)
                    if (!R_finite(cells(row, col))) stop("nonfinite subdivision state");
                if (cells(row, 0) < 0 || cells(row, 0) > 2147483647.0 ||
                    cells(row, 1) < 0 || cells(row, 1) > 30) stop("invalid subdivision index");
                int pair = cells(row, 0), depth = cells(row, 1);
                double lo = cells(row, 2), hi = cells(row, 3), width = std::ldexp(1.0, -depth);
                double vlo = cells(row, 4), vhi = cells(row, 5);
                if (pair != cells(row, 0) || depth != cells(row, 1) || pair < start ||
                    pair >= start + count || depth < 0 || depth > 30 || lo < 0 || hi > 1 ||
                    hi - lo != width || std::floor(lo / width) != lo / width || vlo < 0 || vhi > 1 ||
                    vhi - vlo != width || std::floor(vlo / width) != vlo / width)
                    stop("invalid dyadic subdivision state");
                // Morton intervals encode entire dyadic squares without sampling.
                uint64_t u = lo / width, v = vlo / width, path = 0;
                for (int bit = depth - 1; bit >= 0; --bit)
                    path = (path << 2) | (((u >> bit) & 1) << 1) | ((v >> bit) & 1);
                int shift = 2 * (30 - depth);
                partitions[pair - start].push_back({path << shift, (path + 1) << shift});
                plan.leaves.push_back({pair / static_cast<int>(m.pieces.size()),
                    pair % static_cast<int>(m.pieces.size()), depth, lo, hi, vlo, vhi});
            }
            for (auto &partition : partitions) {
                std::sort(partition.begin(), partition.end());
                uint64_t end = 0;
                for (auto cell : partition) {
                    if (cell.first != end) stop("incomplete or overlapping subdivision state");
                    end = cell.second;
                }
                if (end != (uint64_t(1) << 60)) stop("incomplete subdivision state");
            }
            plan.constructed = (4 * plan.leaves.size() - count) / 3;
            bytes += plan.leaves.size() * sizeof(Box);
            if (bytes > 64 * 1024 * 1024 || !plans.emplace(start, std::move(plan)).second)
                stop("excessive or duplicated subdivision state");
        }
        m.cover_plans = std::move(plans);
        m.plan_bytes = bytes;
    }
    List state(m.cover_plans.size());
    int index = 0;
    for (const auto &entry : m.cover_plans) {
        const auto &plan = entry.second;
        NumericMatrix cells(plan.leaves.size(), 6);
        for (size_t j = 0; j < plan.leaves.size(); ++j) {
            const auto &box = plan.leaves[j];
            cells(j, 0) = box.j * m.pieces.size() + box.k;
            cells(j, 1) = box.depth; cells(j, 2) = box.u0; cells(j, 3) = box.u1;
            cells(j, 4) = box.v0; cells(j, 5) = box.v1;
        }
        state[index++] = List::create(_["start"] = plan.start, _["count"] = plan.count, _["cells"] = cells);
    }
    return state;
}
// [[Rcpp::export]]
SEXP build_piecewise_model_cpp(List pieces, NumericVector counts, double eps, double tau,
                               bool circle_prefix = true, bool radial_feasibility = true,
                               bool sharp_rectangles = false, Nullable<NumericMatrix> draws = R_NilValue) {
    return piecewise_kernel::build_piecewise_model_cpp(pieces, counts, eps, tau, circle_prefix,
                                                       radial_feasibility, sharp_rectangles, draws);
}
// [[Rcpp::export]]
void release_piecewise_model_cpp(SEXP model_ptr) {
    piecewise_kernel::release_piecewise_model_cpp(model_ptr);
}
// [[Rcpp::export]]
SEXP build_piecewise_cover_cpp(SEXP model_ptr, List pieces, NumericVector counts, double eps,
                               double tau, int max_cells = 10000, int max_depth = 8,
                               double direction_tolerance = .1, int pair_start = 0,
                               int pair_count = -1, bool continuous_cones = false,
                               int requested_methods = 3, bool witness_only = false,
                               bool independent_G_cover = true) {
    return piecewise_kernel::build_piecewise_cover_cpp(
        model_ptr, pieces, counts, eps, tau, max_cells, max_depth, direction_tolerance,
        pair_start, pair_count, continuous_cones, requested_methods, witness_only,
        independent_G_cover);
}
// [[Rcpp::export]]
List evaluate_piecewise_draws_cpp(SEXP geometry_ptr, NumericMatrix g, bool shared_lower = true,
                                  Nullable<List> prior = R_NilValue,
                                  Nullable<IntegerVector> draw_ids = R_NilValue,
                                  Nullable<NumericMatrix> cutoff_targets = R_NilValue) {
    return piecewise_kernel::evaluate_piecewise_draws_cpp(geometry_ptr, g, shared_lower, prior, draw_ids, cutoff_targets);
}
// [[Rcpp::export]]
void release_piecewise_cover_cpp(SEXP geometry_ptr) {
    piecewise_kernel::release_piecewise_cover_cpp(geometry_ptr);
}
// [[Rcpp::export]]
List evaluate_witnesses_cpp(SEXP model_ptr, List pieces, NumericVector counts, NumericMatrix points,
                            NumericMatrix g, double eps, double tau) {
    return piecewise_kernel::evaluate_witnesses_cpp(model_ptr, pieces, counts, points, g, eps, tau);
}
// [[Rcpp::export]]
List refine_rank_draws_cpp(SEXP model_ptr, NumericMatrix g, NumericMatrix lower,
                           NumericMatrix upper, IntegerVector ids, double target,
                           int max_cells = 262144, int max_depth = 60) {
    return piecewise_kernel::refine_rank_draws_cpp(model_ptr, g, lower, upper, ids, target, max_cells, max_depth);
}
// [[Rcpp::export]]
List circle_charts_cpp(NumericVector y) { return piecewise_kernel::circle_charts_cpp(y); }
// [[Rcpp::export]]
Rcpp::List ranking_membership_cpp(Rcpp::NumericVector qlo, Rcpp::NumericVector qhi,
                                  Rcpp::NumericMatrix rlo, Rcpp::NumericMatrix rhi,
                                  Rcpp::NumericVector counts, double cutoff) {
    return finite_kernel::ranking_membership_cpp(qlo, qhi, rlo, rhi, counts, cutoff);
}
// [[Rcpp::export]]
NumericMatrix circle_population_mean_bounds_cpp(double lambda) {
    return circle_kernel::circle_population_mean_bounds_cpp(lambda);
}
// [[Rcpp::export]]
List circle_el_bounds_cpp(NumericVector y, NumericVector queries, double cutoff,
                          Nullable<NumericVector> query_upper = R_NilValue) {
    return circle_kernel::circle_el_bounds_cpp(y, queries, cutoff, query_upper);
}
// [[Rcpp::export]]
List circle_membership_cpp(NumericVector y, double q, double cutoff, double q_upper = NA_REAL) {
    return circle_kernel::circle_membership_cpp(y, q, cutoff, q_upper);
}

// [[Rcpp::export]]
List circle_length_cpp(NumericVector y, NumericVector lower, NumericVector upper, double el,
                       int budget=8192, double tolerance=.57) {
    return circle_length::measure(y, lower, upper, el, budget, tolerance);
}
