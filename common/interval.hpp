// [[Rcpp::plugins(cpp17)]]
#define FC_INLINE_BINARY64_STEPS
#include <Rcpp.h>
#include <array>
#include <chrono>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#ifndef FC_INTERVAL_NEXTAFTER_HPP
#define FC_INTERVAL_NEXTAFTER_HPP
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
namespace fc {
// The encoding order of finite NORMAL binary64 numbers is monotone in the
// unsigned magnitude, and reversed for negative values. Move exactly one bit
// pattern. Boundary, subnormal, zero, infinity and NaN cases use the standard
// routine, retaining its exceptional-value and floating-environment behavior.
// Unrecognized encodings also retain the standard routine.
inline double binary64_step(double x, bool toward_positive) {
#if defined(__BYTE_ORDER__) && defined(__FLOAT_WORD_ORDER__) &&                                    \
    __BYTE_ORDER__ == __FLOAT_WORD_ORDER__
    if constexpr (sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits == 53) {
        std::uint64_t bits;
        std::memcpy(&bits, &x, sizeof x);
        std::uint64_t magnitude = bits & UINT64_C(0x7fffffffffffffff);
        if (magnitude > UINT64_C(0x0010000000000000) && magnitude < UINT64_C(0x7fefffffffffffff)) {
            if (toward_positive == (x > 0))
                ++bits;
            else
                --bits;
            double out;
            std::memcpy(&out, &bits, sizeof out);
            return out;
        }
    }
#endif
    return std::nextafter(x, toward_positive ? std::numeric_limits<double>::infinity()
                                             : -std::numeric_limits<double>::infinity());
}
} // namespace fc
#endif

#ifndef FC_INTERVAL_FINITE_HPP
#define FC_INTERVAL_FINITE_HPP
// Elementary binary64 enclosures for a finite target with supplied loss intervals.
// Requires IEC 60559 binary64, round-to-nearest, gradual underflow, no fast-math.
// Nonfinite intermediate results cause a whole interval, never an unsafe decision.
#ifdef __FAST_MATH__
#error "Certified interval kernel must not be built with fast-math"
#endif
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0
#error "Certified interval kernel requires infinity and NaN semantics"
#endif
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
namespace fc {
struct Interval {
    double lo, hi;
};
inline Interval whole() {
    return {-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
}
#ifdef FC_INLINE_BINARY64_STEPS
inline double down(double x) { return binary64_step(x, false); }
inline double up(double x) { return binary64_step(x, true); }
#else
inline double down(double x) { return std::nextafter(x, -std::numeric_limits<double>::infinity()); }
inline double up(double x) { return std::nextafter(x, std::numeric_limits<double>::infinity()); }
#endif
inline bool valid(Interval x) { return !std::isnan(x.lo) && !std::isnan(x.hi) && x.lo <= x.hi; }
inline Interval add(Interval a, Interval b) {
    double lo = a.lo + b.lo, hi = a.hi + b.hi;
    if (!std::isfinite(lo) || !std::isfinite(hi))
        return whole();
    return {down(lo), up(hi)};
}
inline Interval sub(Interval a, Interval b) { return add(a, {-b.hi, -b.lo}); }
inline Interval mul(Interval a, Interval b) {
    if ((a.lo == 0 && a.hi == 0) || (b.lo == 0 && b.hi == 0))
        return {0, 0};
    const double x[4] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    for (double v : x)
        if (!std::isfinite(v))
            return whole();
    return {down(*std::min_element(x, x + 4)), up(*std::max_element(x, x + 4))};
}
inline Interval square(Interval a) {
    if (a.lo == 0 && a.hi == 0)
        return {0, 0};
    const double l = a.lo * a.lo, h = a.hi * a.hi;
    if (!std::isfinite(l) || !std::isfinite(h))
        return {0, std::numeric_limits<double>::infinity()};
    double lo = (a.lo <= 0 && a.hi >= 0) ? 0 : std::max(0.0, down(std::min(l, h)));
    return {lo, up(std::max(l, h))};
}
inline Interval ratio_positive(double a, double b) {
    if (!(a > 0 && b > 0) || !std::isfinite(a) || !std::isfinite(b))
        throw std::invalid_argument("invalid positive ratio");
    const double x = a / b;
    if (!std::isfinite(x))
        return whole();
    return {std::max(0.0, down(x)), up(x)};
}
inline void check_platform() {
    if (!std::numeric_limits<double>::is_iec559 || std::numeric_limits<double>::digits != 53 ||
        std::fegetround() != FE_TONEAREST)
        throw std::runtime_error("binary64 round-to-nearest is required");
    // This check is necessary but does not replace a no-fast-math build contract.
    volatile double tiny = std::numeric_limits<double>::denorm_min();
    volatile double same = tiny * 1.0;
    if (same == 0)
        throw std::runtime_error("gradual underflow is required");
}
enum class Decision { In, Out, Unresolved };
struct PairDecision {
    Decision decision;
    Interval delta;
    Interval margin;
};
inline PairDecision pair_decision(const std::vector<Interval> &q, const std::vector<Interval> &r,
                                  const std::vector<double> &counts, double cutoff) {
    check_platform();
    if (q.empty() || q.size() != r.size() || q.size() != counts.size() || !std::isfinite(cutoff) ||
        cutoff < 0)
        throw std::invalid_argument("invalid query");
    double n = 0;
    bool same_exact_loss = true;
    for (std::size_t i = 0; i < q.size(); ++i) {
        if (!valid(q[i]) || !valid(r[i]) || !std::isfinite(q[i].lo) || !std::isfinite(q[i].hi) ||
            !std::isfinite(r[i].lo) || !std::isfinite(r[i].hi) || !(counts[i] > 0) ||
            !std::isfinite(counts[i]) || std::floor(counts[i]) != counts[i])
            throw std::invalid_argument("invalid loss intervals or counts");
        n += counts[i];
        if (n > 4503599627370496.0)
            throw std::invalid_argument("integer count bound exceeded");
        same_exact_loss =
            same_exact_loss && q[i].lo == q[i].hi && r[i].lo == r[i].hi && q[i].lo == r[i].lo;
    }
    if (same_exact_loss)
        return {Decision::In, {0, 0}, {0, 0}};
    std::vector<Interval> d(q.size()), p(q.size());
    Interval mean{0, 0};
    for (std::size_t i = 0; i < q.size(); ++i) {
        d[i] = sub(q[i], r[i]);
        p[i] = ratio_positive(counts[i], n);
        mean = add(mean, mul(p[i], d[i]));
    }
    if (mean.hi <= 0)
        return {Decision::In, mean, whole()};
    Interval variance{0, 0};
    for (std::size_t i = 0; i < q.size(); ++i)
        variance = add(variance, mul(p[i], square(sub(d[i], mean))));
    variance.lo = std::max(0.0, variance.lo);
    Interval margin = sub(mul(square({cutoff, cutoff}), variance), mul({n, n}, square(mean)));
    if (margin.lo >= 0)
        return {Decision::In, mean, margin};
    if (mean.lo > 0 && margin.hi < 0)
        return {Decision::Out, mean, margin};
    return {Decision::Unresolved, mean, margin};
}
} // namespace fc
#endif
