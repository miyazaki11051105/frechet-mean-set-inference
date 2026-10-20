#pragma once
namespace piecewise_kernel {
namespace cone {
namespace {
I polynomial(const std::vector<I> &c, I x) {
    I value = pt(0);
    for (auto it = c.rbegin(); it != c.rend(); ++it) value = plus(times(value, x), *it);
    return value;
}
// Interval exclusion, rather than approximate polynomial roots. All unexcluded
// bands remain candidates, including multiple roots and identically zero cases.
void root_bands(const std::vector<I> &c, std::vector<I> &out) {
    std::vector<std::pair<I, int>> pending = {{{0, 1}, 0}};
    int visited = 0;
    while (!pending.empty()) {
        auto [x, depth] = pending.back(); pending.pop_back();
        I value = polynomial(c, x);
        if (value.lo > 0 || value.hi < 0) continue;
        if (depth == 12 || ++visited >= 128) { out.push_back(x); continue; }
        double mid = (x.lo + x.hi) / 2;
        pending.push_back({{x.lo, mid}, depth + 1});
        pending.push_back({{mid, x.hi}, depth + 1});
    }
}
}
std::shared_ptr<CircleProfile> Model::circle_profile(const Box &box) {
    bool cache_root = prefix && pieces.size() <= 201 && box.depth == 0 &&
        box.u0 == 0 && box.u1 == 1 && box.v0 == 0 && box.v1 == 1;
    size_t root_key = box.j * pieces.size() + box.k;
    if (cache_root) {
        if (circle_root_profiles.empty()) circle_root_profiles.resize(pieces.size()*pieces.size());
        if (circle_root_profiles[root_key]) {
            return circle_root_profiles[root_key];
        }
    }
    const auto &a = pieces[box.j], &b = pieces[box.k];
    auto out = std::make_shared<CircleProfile>();
    I r = coordinate(b, box.v0, box.v0);
    int ae = 0, be = 0;
    bool shared = circle_shared_endpoint(box.j, box.k, ae, be);
    I length = shared ? times(minus(b.right, b.left), minus(pt(box.v1), pt(box.v0)))
                      : minus(coordinate(b, box.v1, box.v1), r);
    out->intercept = circle_contrast({box.j, box.k, box.depth, box.u0, box.u1, box.v0, box.v0});
    out->slope = prefix->displacement(a.circle_index, b.circle_index, pt(0), length);
    out->variance = {prefix->gram(out->intercept, out->intercept),
                     times_point(prefix->gram(out->intercept, out->slope), 2),
                     prefix->gram(out->slope, out->slope)};
    out->delta = {circle_risk_difference({box.j, box.k, box.depth, box.u0, box.u1, box.v0, box.v0}),
                  neg(times(length, plus(times_point(r, 2), b.mb))), neg(sq(length))};
    out->sample_size = n;
    out->tau_squared = tau2;
    std::vector<I> h(5, pt(0));
    for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
            h[j + k] = plus(h[j + k], times(n, times(out->delta[j], out->delta[k])));
    for (int j = 0; j < 3; ++j) h[j] = minus(h[j], times(tau2, out->variance[j]));
    out->boundaries = {pt(0), pt(1)};
    root_bands(h, out->boundaries);
    root_bands({out->variance.begin(), out->variance.end()}, out->boundaries);
    // delta=0 is harmless to include, and covers all branches of the selector.
    root_bands({out->delta.begin(), out->delta.end()}, out->boundaries);
    std::sort(out->boundaries.begin(), out->boundaries.end(), [](I x, I y) { return x.lo < y.lo; });
    std::vector<I> merged;
    for (I x : out->boundaries) {
        if (!merged.empty() && x.lo <= merged.back().hi) merged.back().hi = std::max(merged.back().hi, x.hi);
        else merged.push_back(x);
    }
    out->boundaries = std::move(merged);
    if (cache_root) {
        circle_root_profiles[root_key] = out;
    }
    return out;
}
std::shared_ptr<CircleProfile> Model::circle_shared_profile(const Box &box) {
    I dq, dr; int ae = 0, be = 0;
    if (!circle_endpoint_offsets(box, dq, dr, ae, be)) return nullptr;
    const auto &a = pieces[box.j], &b = pieces[box.k];
    I s = ae ? neg(dq) : dq, t = be ? neg(dr) : dr;
    s = {std::max(0.0,s.lo),std::max(0.0,s.hi)};
    t = {std::max(0.0,t.lo),std::max(0.0,t.hi)};
    if (s.hi == 0 && t.hi == 0) return nullptr; // Existing zero-variance path.
    I lower_den = plus(pt(s.hi),pt(t.lo)), upper_den = plus(pt(s.lo),pt(t.hi));
    double rho_lo = lower_den.lo > 0 ? std::max(0.0,divide(pt(t.lo),lower_den).lo) : 0;
    double rho_hi = upper_den.lo > 0 ? std::min(1.0,divide(pt(t.hi),upper_den).hi) : 1;
    if (rho_lo > rho_hi) return nullptr; // Retain the previous full-cover cap.
    I span = minus(pt(rho_hi),pt(rho_lo));
    auto Y = prefix->displacement(a.circle_index,b.circle_index,pt(ae ? -1 : 1),pt(0));
    auto Z = prefix->displacement(a.circle_index,b.circle_index,pt(0),pt(be ? -1 : 1));
    auto combine = [](const PrefixVector &x, I cx, const PrefixVector &y, I cy) {
        std::map<std::pair<int,int>,I> coefficients;
        for (const auto &part : std::vector<std::pair<const PrefixVector *,I>>{{&x,cx},{&y,cy}})
            for (const auto &term : *part.first) {
                auto key = std::make_pair(term.cut,term.power);
                I value = times(term.coefficient,part.second);
                auto found = coefficients.find(key);
                if (found == coefficients.end()) coefficients[key] = value;
                else found->second = plus(found->second,value);
            }
        PrefixVector answer;
        for (const auto &term : coefficients)
            answer.push_back({term.first.first,term.first.second,term.second});
        return answer;
    };
    auto out = std::make_shared<CircleProfile>();
    // For s+t>0, w/(s+t)=V(rho)=(1-rho)Y+rho Z.  This preserves
    // the positive-direction cone while removing the vanishing radial scale.
    out->intercept = combine(Y,minus(pt(1),pt(rho_lo)),Z,pt(rho_lo));
    out->slope = combine(Y,neg(span),Z,span);
    out->variance = {prefix->gram(out->intercept,out->intercept),
        times_point(prefix->gram(out->intercept,out->slope),2),prefix->gram(out->slope,out->slope)};
    I ea = ae ? a.right : a.left, eb = be ? b.right : b.left;
    I my = times_point(plus(times_point(ea,2),a.mb),ae ? -1 : 1);
    I mz = times_point(plus(times_point(eb,2),b.mb),be ? 1 : -1);
    I difference = minus(mz,my);
    double D = minus(pt(s.hi),pt(t.lo)).hi;
    // Local requires m(rho)+s-t >= -tau/sqrt(n)*||V||.  Since s-t<=D,
    // m(rho)+D is a necessary-condition relaxation.  Dropping the anchor
    // restriction only enlarges this UPPER domain; it admits no lower witness.
    out->delta = {plus(plus(my,times(difference,pt(rho_lo))),pt(D)),times(difference,span),pt(0)};
    out->sample_size = n;
    out->tau_squared = tau2;
    std::vector<I> h(5,pt(0));
    for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k)
        h[j+k] = plus(h[j+k],times(n,times(out->delta[j],out->delta[k])));
    for (int j = 0; j < 3; ++j) h[j] = minus(h[j],times(tau2,out->variance[j]));
    out->boundaries = {pt(0),pt(1)};
    root_bands(h,out->boundaries);
    root_bands({out->variance.begin(),out->variance.end()},out->boundaries);
    root_bands({out->delta.begin(),out->delta.end()},out->boundaries);
    std::sort(out->boundaries.begin(),out->boundaries.end(),[](I x,I y){return x.lo<y.lo;});
    std::vector<I> merged;
    for (I x : out->boundaries) {
        if (!merged.empty() && x.lo<=merged.back().hi) merged.back().hi=std::max(merged.back().hi,x.hi);
        else merged.push_back(x);
    }
    out->boundaries=std::move(merged);
    return out;
}
double CircleProfile::upper(const CirclePrefix &prefix,
                           const std::vector<std::array<I, 2>> &g,
                           double cap) const {
    I c = prefix.dot(intercept, g), a = prefix.dot(slope, g);
    // 2*N'*Q-N*Q' is linear: the quadratic terms cancel exactly.
    I s0 = minus(times_point(times(a, variance[0]), 2), times(c, variance[1]));
    I s1 = minus(times(a, variance[1]), times_point(times(c, variance[2]), 2));
    double answer = 0;
    auto candidate = [&](I x) {
        I q = plus(variance[0], times(x, plus(variance[1], times(x, variance[2]))));
        I d = plus(delta[0], times(x, plus(delta[1], times(x, delta[2]))));
        if (d.hi < 0 && times(sample_size, sq(d)).lo > times(tau_squared, q).hi) return;
        if (q.hi <= 0) { answer = cap; return; }
        I numerator = plus(c, times(a, x));
        double upper = q.lo > 0 ? std::min(cap, divide(numerator, root(q)).hi) : cap;
        answer = std::max(answer, upper);
    };
    for (I x : boundaries) candidate(x);
    if (s1.lo > 0 || s1.hi < 0) {
        I stationary = times(neg(s0), {fc::down(1 / s1.hi), fc::up(1 / s1.lo)});
        stationary = {std::max(0.0, stationary.lo), std::min(1.0, stationary.hi)};
        if (stationary.lo <= stationary.hi) candidate(stationary);
    } else {
        I derivative = plus(s0, times(s1, {0, 1}));
        if (derivative.lo <= 0 && derivative.hi >= 0) candidate({0, 1});
    }
    return std::min(cap, answer);
}
Evaluation Model::circle_evaluation(const Box &box, bool upper_geometry) {
    if (prefix->constant_sample ||
        (box.j == box.k && box.u0 == box.u1 && box.v0 == box.v1 && box.u0 == box.v0))
        return {{}, 0, 0, 0};
    const auto &a = pieces[box.j], &b = pieces[box.k];
    I t = coordinate(a, box.u0, box.u1), r = coordinate(b, box.v0, box.v1);
    I anchor = minus(risk(a, t), fmin);
    if (box.j == box.k && a.vb.lo > 0) {
        std::vector<Direction> out;
        int possible = 0, certain = 3;
        for (int sign : {-1, 1}) {
            if ((sign < 0 && box.u0 >= box.v1) || (sign > 0 && box.u1 <= box.v0))
                continue;
            I scaled = times_point(plus(plus(t, r), a.mb), sign);
            int pm = select(scaled, a.vb, anchor, false);
            int cm = select(scaled, a.vb, anchor, true);
            possible |= pm;
            certain &= cm;
            Direction d;
            d.mask = pm;
            d.sphere = false;
            I pi = {0x1.921fb54442d18p+1, 0x1.921fb54442d19p+1}, sd = root(a.vb);
            d.compact = {{size, 1, divide(pt(-2 * sign), sd)},
                         {a.circle_index, 0, divide(times(pt(-4 * sign), pi), sd)}};
            out.push_back(std::move(d));
        }
        return {out, possible, certain & possible, 0};
    }
    auto contrast = circle_contrast(box);
    I variance = prefix->gram(contrast, contrast);
    // Convexity of the norm bounds its maximum by the four corner norms.
    double corner_upper = 0, radius2 = 0;
    double um = (box.u0 + box.u1) / 2, vm = (box.v0 + box.v1) / 2;
    I tc = coordinate(a, um, um), rc = coordinate(b, vm, vm);
    auto center = circle_contrast({box.j, box.k, box.depth, um, um, vm, vm});
    I center2 = prefix->gram(center, center);
    for (double u : {box.u0, box.u1})
        for (double v : {box.v0, box.v1}) {
            auto w = circle_contrast({box.j, box.k, box.depth, u, u, v, v});
            corner_upper = std::max(corner_upper, prefix->gram(w, w).hi);
            int ae = 0, be = 0;
            bool shared = circle_shared_endpoint(box.j, box.k, ae, be);
            I dq = shared ? times(minus(a.right, a.left), minus(pt(u), pt(um)))
                          : minus(coordinate(a, u, u), tc);
            I dr = shared ? times(minus(b.right, b.left), minus(pt(v), pt(vm)))
                          : minus(coordinate(b, v, v), rc);
            w = prefix->displacement(a.circle_index, b.circle_index, dq, dr);
            radius2 = std::max(radius2, prefix->gram(w, w).hi);
        }
    variance = {std::max(0.0, variance.lo), std::min(variance.hi, corner_upper)};
    if (variance.lo > variance.hi)
        stop("inconsistent compact variance enclosure");
    bool local_excluded = false;
    // At a shared endpoint, centered contrasts are exactly s*y+t*z,
    // s,t>=0, while delta=s*my+t*mz+s^2-t^2. If both coefficients
    // are below -tau/sqrt(n) times their ray norms, the entire Local
    // wedge is infeasible. This avoids cancellation as s,t approach zero.
    I threshold = divide(root(tau2), root(n));
    for (int ae = 0; ae < 2; ++ae) for (int be = 0; be < 2; ++be) {
        int ai = ae ? a.right_id : a.left_id, bi = be ? b.right_id : b.left_id;
        if (ai < 0 || ai != bi) continue;
        I endpoint_a = ae ? a.right : a.left, endpoint_b = be ? b.right : b.left;
        I my = times_point(plus(times_point(endpoint_a, 2), a.mb), ae ? -1 : 1);
        I mz = times_point(plus(times_point(endpoint_b, 2), b.mb), be ? 1 : -1);
        I offset = circle_endpoint_offset(a, box.u0, box.u1, ae);
        I distance = ae ? neg(offset) : offset;
        I q_coefficient = plus(plus(my, pt(std::max(0.0, distance.hi))), times(threshold, root(a.vb)));
        I r_coefficient = plus(mz, times(threshold, root(b.vb)));
        if (q_coefficient.hi < 0 && r_coefficient.hi < 0) local_excluded = true;
    }
    // Re-evaluate remaining small cells when prefix cancellation obscures variance.
    if (!local_excluded && variance.lo <= 0 &&
        (box.depth >= 8 || (box.u0 == box.u1 && box.v0 == box.v1)))
        return evaluate(box, upper_geometry, false);
    I delta = circle_risk_difference(box);
    int possible = select(delta, variance, anchor, false);
    int certain = select(delta, variance, anchor, true);
    if (local_excluded) { possible &= 1; certain &= 1; }
    if (!possible)
        return {{}, 0, 0, 0};
    Direction d;
    d.mask = possible;
    d.sphere = variance.lo <= 0;
    if (!d.sphere) {
        d.compact = std::move(contrast);
        I sd = root(variance);
        for (auto &term : d.compact)
            term.coefficient = divide(term.coefficient, sd);
    }
    if (upper_geometry && box.j != box.k) {
        d.compact_basis = &compressed_basis(box.j, box.k);
        if (d.compact_basis->empty()) {
            d.compact_basis = nullptr;
            d.basis = &basis(box.j, box.k);
        }
        // Adjacent-chart Local halfspaces remain valid away from the vertex.
        for (int ae = 0; ae < 2; ++ae)
            for (int be = 0; be < 2; ++be) {
                int ai = ae ? a.right_id : a.left_id, bi = be ? b.right_id : b.left_id;
                if (ai >= 0 && ai == bi && (possible & 2)) {
                    d.local = local_limit(box.j, box.k, ae, be);
                    if (d.local) {
                        I dq = circle_endpoint_offset(a, box.u0, box.u1, ae);
                        I dr = circle_endpoint_offset(b, box.v0, box.v1, be);
                        I remainder = plus(
                            times_point(d.local->dual_y, std::max(0.0, (ae ? neg(dq) : dq).hi)),
                            times_point(d.local->dual_z, std::max(0.0, (be ? neg(dr) : dr).hi)));
                        d.local_eta = plus(divide(root(tau2), root(n)), remainder).hi;
                    }
                }
            }
    }
    // ||normalize(w)-normalize(center)|| <= 2||w-center||/||center||.
    double width = center2.lo > 0
                       ? std::min(2.0, divide(times_point(root({0, radius2}), 4), root(center2)).hi)
                       : 2;
    return {{std::move(d)}, possible, certain, width};
}
} // namespace cone
List circle_charts_cpp(NumericVector y);

// Internal accumulator: targets are newly allocated, private matrices owned by
// piecewise_bounds. Intentionally update them in place (never pass user
// arrays here). No summation or floating-point operation order is changed.

void merge_piecewise_bounds_cpp(NumericMatrix lower, NumericMatrix upper, NumericMatrix chunk_lower,
                                NumericMatrix chunk_upper, int offset) {
    int n = chunk_lower.nrow(), B = lower.nrow();
    if (offset < 0 || offset > B - n || upper.nrow() != B || chunk_upper.nrow() != n ||
        lower.ncol() != 2 || upper.ncol() != 2 || chunk_lower.ncol() != 2 ||
        chunk_upper.ncol() != 2)
        stop("invalid accumulator dimensions");
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < n; ++i) {
            // Match pmax including missing values; ordinary certified draws have none.
            double a = chunk_lower(i, j), b = chunk_upper(i, j);
            lower(offset + i, j) = ISNAN(a) ? a : std::max(lower(offset + i, j), a);
            upper(offset + i, j) = ISNAN(b) ? b : std::max(upper(offset + i, j), b);
        }
}
} // namespace piecewise_kernel
namespace circle_kernel {

// Opt-in fixed-cutoff size audit. No multiplier calibration is performed.
// Binary64 inputs denote exact real numbers; pi is enclosed, not rounded away.
// Shared outward-rounded chart polynomials for Coverage membership and EL.

using namespace Rcpp;
namespace {
using V = fc::Interval;
V pt(double x) { return {x, x}; }
V neg(V x) { return {-x.hi, -x.lo}; }
V ad(V x, V y) { return fc::add(x, y); }
V su(V x, V y) { return fc::sub(x, y); }
V mu(V x, V y) { return fc::mul(x, y); }
V sq(V x) { return fc::square(x); }
V divi(V x, V y) {
    if (y.lo <= 0 && y.hi >= 0)
        return fc::whole();
    V inverse{fc::down(1 / y.hi), fc::up(1 / y.lo)};
    return mu(x, inverse);
}
V scale(V x, double y) { return mu(x, pt(y)); }
double mag(V x) { return std::max(std::abs(x.lo), std::abs(x.hi)); }
const double P = 3.141592653589793238462643383279502884;
V pi() { return {fc::down(P), fc::up(P)}; }
V interp(V a, V b, double t) { return ad(scale(a, 1 - t), scale(b, t)); }
V span(V a, V b) { return {a.lo, b.hi}; }

// log via convergent atanh series after exact binary exponent reduction.
// No libm log accuracy assumption; the finite-series remainder is enclosed.
V log_unit(V z) {
    V power = z, sum = pt(0), z2 = sq(z);
    for (int j = 0; j < 24; ++j) {
        sum = ad(sum, divi(power, pt(2 * j + 1)));
        power = mu(power, z2);
    }
    double a = mag(z);
    V absz{a, a}, power49 = pt(mag(power)); // The loop already encloses z^49.
    V remainder = divi(scale(power49, 2), scale(su(pt(1), sq(absz)), 49));
    V result = scale(sum, 2);
    return ad(result, {-remainder.hi, remainder.hi});
}
V log_point(double x) {
    if (!(x > 0) || !std::isfinite(x))
        return fc::whole();
    int exponent = 0;
    double mantissa = std::frexp(x, &exponent);
    V z = divi(su(pt(mantissa), pt(1)), ad(pt(mantissa), pt(1)));
    static const V log2 = log_unit(divi(pt(1), pt(3)));
    return ad(log_unit(z), scale(log2, exponent));
}
V logarithm(V x) {
    if (!(x.lo > 0))
        return fc::whole();
    V a = log_point(x.lo), b = log_point(x.hi);
    return {a.lo, b.hi};
}
int choose(int n, int k) {
    static const int c[5][5] = {
        {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}};
    return c[n][k];
}
std::vector<V> bernstein_coefficients(const std::vector<V> &coefficients, V left, V right) {
    int d = coefficients.size() - 1;
    std::vector<V> lp(d + 1, pt(1)), wp(d + 1, pt(1)), power(d + 1, pt(0));
    V width = su(right, left);
    for (int i = 1; i <= d; ++i) {
        lp[i] = mu(lp[i - 1], left);
        wp[i] = mu(wp[i - 1], width);
    }
    for (int k = 0; k <= d; ++k)
        for (int j = k; j <= d; ++j)
            if (mag(coefficients[j]) != 0)
                power[k] = ad(power[k], scale(mu(mu(coefficients[j], lp[j - k]), wp[k]), choose(j, k)));
    std::vector<V> answer;
    for (int j = 0; j <= d; ++j) {
        V b = pt(0);
        for (int k = 0; k <= j; ++k)
            b = ad(b, mu(power[k], divi(pt(choose(j, k)), pt(choose(d, k)))));
        answer.push_back(b);
    }
    return answer;
}
using Matrix = std::vector<std::vector<V>>;
Matrix q_bernstein(const Matrix &coefficients, V q) {
    int degree = coefficients.size() - 1;
    Matrix out(degree + 1, std::vector<V>(degree + 1, pt(0)));
    for (int k = 0; k <= degree; ++k) {
        std::vector<V> a;
        for (int j = 0; j <= degree; ++j)
            a.push_back(coefficients[j][k]);
        auto b = bernstein_coefficients(a, pt(q.lo), pt(q.hi));
        for (int j = 0; j <= degree; ++j)
            out[j][k] = b[j];
    }
    return out;
}
Matrix r_bernstein(const Matrix &qcoeff, V left, V right) {
    Matrix out;
    for (const auto &row : qcoeff) out.push_back(bernstein_coefficients(row, left, right));
    return out;
}
V coefficient_range(const Matrix &coeff) {
    V out{INFINITY, -INFINITY};
    for (const auto &row : coeff) for (V b : row) {
        out.lo = std::min(out.lo, b.lo); out.hi = std::max(out.hi, b.hi);
    }
    return out;
}
std::pair<Matrix, Matrix> bisect_r(const Matrix &coeff) {
    Matrix left = coeff, right = coeff;
    for (size_t i = 0; i < coeff.size(); ++i) {
        auto work = coeff[i];
        int degree = work.size() - 1;
        for (int k = 1; k <= degree; ++k) {
            for (int j = 0; j <= degree - k; ++j) work[j] = scale(ad(work[j], work[j + 1]), .5);
            left[i][k] = work[0]; right[i][degree - k] = work[degree - k];
        }
    }
    return {std::move(left), std::move(right)};
}
struct Arc {
    V left, right, mean, second, variance, origin, local_mean, local_second;
    std::vector<V> lifts, local_lifts;
    int cut;
};
struct Pair {
    std::array<V, 6> variance;
    Matrix risk, negative_square;
};
struct QueryPair {
    Matrix risk, variance, negative_square;
    bool ready = false;
};
struct Data {
    int n;
    std::vector<double> y;
    std::vector<int> counts;
    std::vector<V> weights;
    std::vector<Arc> arcs;
    std::vector<std::vector<Pair>> cache;
    std::vector<bool> ready;
    std::vector<int> prefix_counts;
    std::vector<std::array<V, 2>> prefix;
    std::vector<QueryPair> query_pairs;
    int query_home = -1;
    int dense_pairs = 0;
    V query_q{NAN, NAN};
    V mean(const std::vector<V> &x) {
        V s = pt(0);
        for (size_t i = 0; i < x.size(); ++i)
            s = ad(s, mu(weights[i], x[i]));
        return s;
    }
    Data(NumericVector input) : n(input.size()) {
        if (n < 3)
            stop("at least three observations required");
        std::vector<double> sorted = as<std::vector<double>>(input);
        std::sort(sorted.begin(), sorted.end());
        for (double x : sorted) {
            if (!std::isfinite(x) || std::abs(x) >= pi().lo)
                stop("angles must be in (-pi,pi)");
            if (y.empty() || x != y.back()) {
                y.push_back(x);
                counts.push_back(1);
            } else
                ++counts.back();
        }
        for (int c : counts)
            weights.push_back(divi(pt(c), pt(n)));
        std::vector<int> order;
        for (size_t i = 0; i < y.size(); ++i)
            if (y[i] >= 0)
                order.push_back(i);
        for (size_t i = 0; i < y.size(); ++i)
            if (y[i] < 0)
                order.push_back(i);
        std::vector<V> lifts;
        for (double x : y)
            lifts.push_back(ad(pt(x), scale(pi(), x >= 0 ? -2 : 0)));
        prefix_counts.push_back(0);
        prefix.push_back({pt(0), pt(0)});
        for (int i : order) {
            prefix_counts.push_back(prefix_counts.back() + counts[i]);
            prefix.push_back({ad(prefix.back()[0], scale(lifts[i], counts[i])),
                              ad(prefix.back()[1], scale(sq(lifts[i]), counts[i]))});
        }
        V left = neg(pi());
        for (size_t j = 0; j <= order.size(); ++j) {
            V right = j < order.size() ? ad(pt(y[order[j]]), scale(pi(), y[order[j]] >= 0 ? -1 : 1))
                                       : pi();
            // The sole identically zero first arc is caused by an observation at 0.
            bool zero = j == 0 && j < order.size() && y[order[j]] == 0;
            if (!zero) {
                Arc a;
                a.cut = j;
                a.left = left;
                a.right = right;
                a.lifts = lifts;
                a.mean = mean(lifts);
                std::vector<V> seconds, centered;
                for (V x : lifts) {
                    seconds.push_back(sq(x));
                    centered.push_back(sq(su(x, a.mean)));
                }
                a.second = mean(seconds);
                a.variance = mean(centered);
                a.origin = interp(left, right, .5);
                std::vector<V> local_seconds;
                for (V x : lifts) {
                    V v = su(x, a.origin);
                    a.local_lifts.push_back(v);
                    local_seconds.push_back(sq(v));
                }
                a.local_mean = mean(a.local_lifts);
                a.local_second = mean(local_seconds);
                arcs.push_back(a);
            }
            if (j < order.size())
                lifts[order[j]] = ad(lifts[order[j]], scale(pi(), 2));
            left = right;
        }
        cache.resize(arcs.size());
        ready.assign(arcs.size(), false);
        query_pairs.resize(arcs.size());
    }
    void prepare(int home) {
        if (ready[home])
            return;
        const Arc &q = arcs[home];
        for (const Arc &r : arcs) {
            Pair z;
            z.variance.fill(pt(0));
            // In antipode order both lifts are x + 2*pi*1{i<cut}.
            // Their centered contrast coefficients are affine in x on at
            // most three segments: the x^2 terms cancel algebraically.
            std::array<int, 4> cuts{0, std::min(q.cut, r.cut), std::max(q.cut, r.cut), int(y.size())};
            for (int k = 0; k < 3; ++k) {
                int l = cuts[k], h = cuts[k + 1];
                if (l == h) continue;
                V s0 = divi(pt(prefix_counts[h] - prefix_counts[l]), pt(n)),
                  s1 = divi(su(prefix[h][0], prefix[l][0]), pt(n)),
                  s2 = divi(su(prefix[h][1], prefix[l][1]), pt(n));
                V aq = su(scale(pi(), l < q.cut ? 2 : 0), q.origin),
                  ar = su(scale(pi(), l < r.cut ? 2 : 0), r.origin);
                std::array<V, 2> a{scale(su(aq, q.local_mean), -2), pt(-2)},
                                 b{scale(su(ar, r.local_mean), 2), pt(2)},
                                 c{su(su(sq(aq), sq(ar)), su(q.local_second, r.local_second)),
                                   scale(su(aq, ar), 2)};
                auto product = [&](const auto &u, const auto &v) {
                    return ad(ad(mu(mu(u[0], v[0]), s0),
                                 mu(ad(mu(u[0], v[1]), mu(u[1], v[0])), s1)),
                              mu(mu(u[1], v[1]), s2));
                };
                std::array<V, 6> terms{product(a, a), scale(product(a, b), 2), product(b, b),
                                      scale(product(a, c), 2), scale(product(b, c), 2), product(c, c)};
                for (int t = 0; t < 6; ++t) z.variance[t] = ad(z.variance[t], terms[t]);
            }
            bool ill_conditioned = false;
            for (V x : z.variance)
                ill_conditioned |= x.hi - x.lo > 1e-6 * std::max(1e-24, mag(x));
            if (&q != &r && ill_conditioned) {
                ++dense_pairs;
                std::array<std::vector<V>, 6> terms;
                for (size_t i = 0; i < y.size(); ++i) {
                    V a = scale(su(q.local_lifts[i], q.local_mean), -2),
                      b = scale(su(r.local_lifts[i], r.local_mean), 2),
                      c = su(su(sq(q.local_lifts[i]), q.local_second),
                             su(sq(r.local_lifts[i]), r.local_second));
                    std::array<V, 6> values{sq(a),scale(mu(a,b),2),sq(b),scale(mu(a,c),2),scale(mu(b,c),2),sq(c)};
                    for (int t = 0; t < 6; ++t) terms[t].push_back(values[t]);
                }
                for (int t = 0; t < 6; ++t) {
                    V dense = mean(terms[t]);
                    z.variance[t] = {std::max(z.variance[t].lo,dense.lo),std::min(z.variance[t].hi,dense.hi)};
                    if (!fc::valid(z.variance[t])) stop("inconsistent circle variance enclosure");
                }
            }
            for (int t : {0,2,5}) z.variance[t].lo = std::max(0.0,z.variance[t].lo);
            z.risk = Matrix(3, std::vector<V>(3, pt(0)));
            z.risk[0][0] = su(q.local_second, r.local_second);
            z.risk[1][0] = scale(q.local_mean, -2);
            z.risk[2][0] = pt(1);
            z.risk[0][1] = scale(r.local_mean, 2);
            z.risk[0][2] = pt(-1);
            z.negative_square = Matrix(5, std::vector<V>(5, pt(0)));
            for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
                for (int c = 0; c < 3; ++c) for (int d = 0; d < 3; ++d)
                    if (mag(z.risk[a][b]) && mag(z.risk[c][d]))
                        z.negative_square[a+c][b+d] = su(z.negative_square[a+c][b+d],
                            scale(mu(z.risk[a][b], z.risk[c][d]), n));
            cache[home].push_back(z);
        }
        ready[home] = true;
    }
    // 1=all interior IN, -1=all interior OUT, 0=unknown. Arc endpoints have measure zero.
    int classify(int home, V q, double cutoff, int r_budget, int r_depth, uint64_t &visits) {
        const Arc &a = arcs[home];
        V c2 = sq(pt(cutoff));
        V diagonal = su(mu(c2, a.variance), scale(sq(su(q, a.mean)), n));
        if (diagonal.hi < 0)
            return -1;
        q = su(q, a.origin);
        bool unknown = diagonal.lo < 0;
        prepare(home);
        if (query_home != home || query_q.lo != q.lo || query_q.hi != q.hi) {
            query_home = home; query_q = q;
            for (auto &p : query_pairs) p.ready = false;
        }
        for (size_t j = 0; j < arcs.size(); ++j) {
            if ((int)j == home)
                continue;
            const Pair &p = cache[home][j];
            const auto &v = p.variance;
            auto &qp = query_pairs[j];
            if (!qp.ready) {
                Matrix variance(5, std::vector<V>(5, pt(0)));
                variance[2][0] = v[0]; variance[1][1] = v[1]; variance[0][2] = v[2];
                variance[1][0] = v[3]; variance[0][1] = v[4]; variance[0][0] = v[5];
                V left = su(arcs[j].left, arcs[j].origin), right = su(arcs[j].right, arcs[j].origin);
                qp.risk = r_bernstein(q_bernstein(p.risk, q), left, right);
                qp.variance = r_bernstein(q_bernstein(variance, q), left, right);
                qp.negative_square = r_bernstein(q_bernstein(p.negative_square, q), left, right);
                qp.ready = true;
            }
            const Matrix &riskq = qp.risk;
            Matrix marginq = qp.negative_square;
            for (int k = 0; k < 5; ++k) for (int l = 0; l < 5; ++l)
                marginq[k][l] = ad(marginq[k][l], mu(c2, qp.variance[k][l]));
            struct Node {
                int depth;
                Matrix risk, margin;
            };
            std::vector<Node> stack;
            stack.push_back({0, riskq, std::move(marginq)});
            int local = 0;
            while (!stack.empty()) {
                Node u = std::move(stack.back());
                stack.pop_back();
                if (local >= r_budget) {
                    unknown = true;
                    break;
                }
                ++local; ++visits;
                V pr = coefficient_range(u.risk), kr = coefficient_range(u.margin);
                if (pr.hi <= 0 || kr.lo >= 0)
                    continue;
                auto rs = bisect_r(u.risk), ms = bisect_r(u.margin);
                V pm{INFINITY, -INFINITY}, km{INFINITY, -INFINITY};
                for (const auto &row : rs.first) pm.lo = std::min(pm.lo, row.back().lo);
                for (const auto &row : ms.first) km.hi = std::max(km.hi, row.back().hi);
                if (pm.lo > 0 && km.hi < 0)
                    return -1;
                if (u.depth >= r_depth) {
                    unknown = true;
                    continue;
                }
                stack.push_back({u.depth + 1, std::move(rs.first), std::move(ms.first)});
                stack.push_back({u.depth + 1, std::move(rs.second), std::move(ms.second)});
            }
        }
        return unknown ? 0 : 1;
    }
    // Enclose the ordinary chart mean EL statistic; every dual bracket is validated.
    V el(const Arc &arc, V q, double cutoff = INFINITY) {
        std::vector<V> z;
        double maximum = -INFINITY, minimum = INFINITY;
        for (V x : arc.lifts) {
            V v = su(x, q);
            z.push_back(v);
            maximum = std::max(maximum, v.hi);
            minimum = std::min(minimum, v.lo);
        }
        double min_lower = INFINITY, max_upper = -INFINITY;
        for (V v : z) {
            min_lower = std::min(min_lower, v.lo);
            max_upper = std::max(max_upper, v.hi);
        }
        if (min_lower > 0 || max_upper < 0)
            return {INFINITY, INFINITY};
        if (!(maximum > 0 && minimum < 0))
            return fc::whole();
        auto gradient = [&](double t) {
            V sum = pt(0);
            for (size_t i = 0; i < z.size(); ++i) {
                V v = z[i], den = ad(pt(1), scale(v, t));
                if (den.lo <= 0)
                    return fc::whole();
                sum = ad(sum, scale(divi(v, den), counts[i]));
            }
            return sum;
        };
        double lo = (-1 / maximum) * (1 - 1e-10), hi = (-1 / minimum) * (1 - 1e-10);
        // Domain safety is checked independently of these floating proposals.
        V gl = gradient(lo), gh = gradient(hi);
        if (!(gl.lo > 0 && gh.hi < 0))
            return fc::whole();
        auto statistic = [&]() {
            V answer = pt(0), root{lo, hi};
            for (size_t i = 0; i < z.size(); ++i) {
                V den = ad(pt(1), mu(root, z[i]));
                if (den.lo <= 0) return fc::whole();
                answer = ad(answer, scale(logarithm(den), 2 * counts[i]));
            }
            answer.lo = std::max(0.0, answer.lo);
            return answer;
        };
        for (int it = 0; it < 70; ++it) {
            double m = lo + (hi - lo) / 2;
            if (m == lo || m == hi)
                break;
            V g = gradient(m);
            if (g.lo > 0)
                lo = m;
            else if (g.hi < 0)
                hi = m;
            else {
                // A bracketed root can land exactly at m. Bound its distance using
                // a verified uniform lower bound on the derivative magnitude.
                V slope = pt(0);
                for (size_t i = 0; i < z.size(); ++i) {
                    V den = ad(pt(1), mu(V{lo, hi}, z[i]));
                    slope = ad(slope, scale(divi(pt(sq(z[i]).lo), sq(pt(den.hi))), counts[i]));
                }
                if (!(slope.lo > 0))
                    break;
                double radius = fc::up(mag(g) / slope.lo), a = std::max(lo, fc::down(m - radius)),
                       b = std::min(hi, fc::up(m + radius));
                if (!(a <= b) || (a == lo && b == hi))
                    break;
                lo = a;
                hi = b;
            }
            if (it == 15 && std::isfinite(cutoff)) {
                V value = statistic();
                if (value.hi <= cutoff || value.lo > cutoff) return value;
            }
        }
        return statistic();
    }
    // Within a cut chart, mean EL is convex in q and minimized at the chart mean.
    // Endpoints certify inclusion; the endpoint nearer that mean certifies exclusion.
    int classify_el(const Arc &arc, V left, V right, double cutoff) {
        V a = el(arc, left, cutoff), b = el(arc, right, cutoff);
        if (a.hi <= cutoff && b.hi <= cutoff)
            return 1;
        if (right.hi < arc.mean.lo && b.lo > cutoff)
            return -1;
        if (left.lo > arc.mean.hi && a.lo > cutoff)
            return -1;
        return 0;
    }
};
} // namespace

// True-pi, outward-rounded ordinary mean-score EL. No G/L optimization here.

NumericMatrix circle_population_mean_bounds_cpp(double lambda) {
    fc::check_platform();
    if (!R_finite(lambda) || lambda > 1 || lambda <= -1)
        stop("invalid circle lambda");
    if (lambda >= 0) return NumericMatrix(1, 2);
    V target = ad(pt(1), divi(scale(pt(lambda), 3), scale(su(pt(1), pt(lambda)), 4)));
    // On [0,1], alternating Taylor terms decrease. Enclose the degree-24
    // polynomial and the next term; no libm acos value is used as evidence.
    auto cosine = [](double x) {
        V term = pt(1), sum = term, x2 = sq(pt(x));
        for (int j = 1; j <= 12; ++j) {
            term = neg(divi(mu(term, x2), pt((2 * j - 1) * (2 * j))));
            sum = ad(sum, term);
        }
        V remainder = divi(mu(term, x2), pt(25 * 26));
        return V{fc::down(sum.lo - remainder.hi), sum.hi};
    };
    if (target.lo <= cosine(1).hi) stop("mean enclosure requires an angle below one radian");
    double lo = 0, hi = 1;
    for (int j = 0; j < 60; ++j) {
        double mid = lo + (hi - lo) / 2;
        if (mid == lo || mid == hi) break;
        V c = cosine(mid);
        if (c.lo > target.hi) lo = mid;
        else if (c.hi < target.lo) hi = mid;
        else break;
    }
    NumericMatrix out(2, 2);
    out(0, 0) = -hi; out(0, 1) = -lo;
    out(1, 0) = lo; out(1, 1) = hi;
    return out;
}

List circle_el_bounds_cpp(NumericVector y, NumericVector queries, double cutoff,
                          Nullable<NumericVector> query_upper = R_NilValue) {
    fc::check_platform();
    if (!R_finite(cutoff) || cutoff < 0)
        stop("invalid EL cutoff");
    Data d(y);
    NumericVector ends = query_upper.isNotNull() ? NumericVector(query_upper) : queries;
    if (ends.size() != queries.size()) stop("query enclosure size mismatch");
    NumericVector lower(queries.size()), upper(queries.size());
    CharacterVector decision(queries.size());
    for (int k = 0; k < queries.size(); ++k) {
        V q{queries[k], ends[k]};
        if (!R_finite(q.lo) || !R_finite(q.hi) || q.lo > q.hi ||
            q.lo <= -pi().lo || q.hi >= pi().lo)
            stop("EL query must be inside (-pi,pi)");
        V e = fc::whole();
        bool zero = true;
        for (double a : y)
            zero = zero && a == q.lo && a == q.hi;
        if (zero)
            e = pt(0);
        else
            for (const Arc &a : d.arcs)
                if (a.left.hi < q.lo && q.hi < a.right.lo) {
                    e = d.el(a, q);
                    break;
                }
        lower[k] = std::max(0.0, e.lo);
        upper[k] = e.hi;
        decision[k] = e.hi <= cutoff ? "IN" : e.lo > cutoff ? "OUT" : "UNRESOLVED";
    }
    return List::create(_["lower"] = lower, _["upper"] = upper, _["decision"] = decision,
                        _["assurance"] =
                            "TRUE_PI_INTERVAL_EL_STATISTIC_AT_ENCLOSED_QUERY_AND_FIXED_BINARY64_CUTOFF");
}

List circle_membership_cpp(NumericVector y, double q, double cutoff, double q_upper = NA_REAL) {
    fc::check_platform();
    if (R_IsNA(q_upper)) q_upper = q;
    if (!R_finite(q) || !R_finite(q_upper) || q > q_upper || !R_finite(cutoff) || cutoff < 0)
        stop("invalid circle query");
    Data d(y);
    uint64_t visits = 0;
    int answer = 0;
    for (size_t h = 0; h < d.arcs.size(); ++h)
        if (d.arcs[h].left.hi < q && q_upper < d.arcs[h].right.lo) {
            answer = d.classify(h, {q, q_upper}, cutoff, 4096, 24, visits);
            break;
        }
    return List::create(_["decision"] = answer == 1    ? "IN"
                                        : answer == -1 ? "OUT"
                                                       : "UNRESOLVED",
                        _["r_cells"] = double(visits));
}
} // namespace circle_kernel
