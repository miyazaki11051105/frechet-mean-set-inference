"""Exact arithmetic for calibration, membership and cutoff audits."""


def index():
    import sys
    from fractions import Fraction

    value = (int(sys.argv[1]) + 1) * (1 - Fraction(float(sys.argv[2])))
    print(-(-value.numerator // value.denominator))


def euclidean():
    """Exact binary64-input support algebra; standard library only, no RNG."""
    import csv
    import sys
    import math
    from fractions import Fraction as F

    def support(y):
        """Independent columns of centered Y, proved using exact row differences."""
        d = len(y[0])
        a = [[v - u for v, u in zip(row, y[0])] for row in y[1:]]
        pivots, row = [], 0
        for col in range(d):
            pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
            if pivot is None:
                continue
            a[row], a[pivot] = a[pivot], a[row]
            divisor = a[row][col]
            a[row] = [v / divisor for v in a[row]]
            for i in range(row + 1, len(a)):
                ratio = a[i][col]
                a[i] = [v - ratio * u for v, u in zip(a[i], a[row])]
            pivots.append(col)
            row += 1
            if row == len(a):
                break
        return pivots

    def float_bounds(value):
        """Tight binary64 enclosure of one exact rational (no decimal truncation)."""
        x = float(value)
        lo = x if F(x) <= value else math.nextafter(x, -math.inf)
        hi = x if F(x) >= value else math.nextafter(x, math.inf)
        return lo, hi

    def basis(y):
        """Exact centered column space, with rational Gram-Schmidt preconditioning.

        Scaling each nonzero vector by its largest exact component prevents tiny
        residual directions from being lost to cancellation or underflow. The
        exported component intervals enclose this EXACT basis, including all ranks.
        """
        n = len(y)
        vectors = []
        for j in support(y):
            mean = sum(row[j] for row in y) / n
            v = [row[j] - mean for row in y]
            for q in vectors:
                coefficient = sum(a * b for a, b in zip(v, q)) / sum(a * a for a in q)
                v = [a - coefficient * b for a, b in zip(v, q)]
            scale = max(abs(a) for a in v)
            if not scale:
                raise ArithmeticError("independent support column became zero")
            vectors.append([a / scale for a in v])
        return vectors

    def euclidean_query(y, q, cutoff):
        """Wald inversion including exact affine support and rank zero."""
        n, d = len(y), len(y[0])
        mean = [sum(row[j] for row in y) / n for j in range(d)]
        z = [[v - m for v, m in zip(row, mean)] for row in y]
        delta = [v - m for v, m in zip(q, mean)]
        ids = support(y)
        if not ids:
            return "IN" if not any(delta) else "OUT"
        # Solve covariance[:, ids] beta = delta on independent coordinates;
        # then check ALL coordinates for affine-support membership.
        cov = [[sum(row[i] * row[j] for row in z) / n for j in ids] for i in range(d)]
        a = [cov[i][:] + [delta[i]] for i in ids]
        for j in range(len(ids)):
            divisor = a[j][j]
            a[j] = [v / divisor for v in a[j]]
            for i in range(len(ids)):
                if i != j:
                    ratio = a[i][j]
                    a[i] = [v - ratio * u for v, u in zip(a[i], a[j])]
        beta = [row[-1] for row in a]
        if any(
            sum(v * b for v, b in zip(row, beta)) != x for row, x in zip(cov, delta)
        ):
            return "OUT"
        statistic2 = n * sum(delta[i] * b for i, b in zip(ids, beta))
        return "IN" if statistic2 <= cutoff * cutoff else "OUT"

    def main():
        with open(sys.argv[1], newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        names = [k for k in rows[0] if k.startswith("y")]
        y = [[F(float(row[k])) for k in names] for row in rows]
        if sys.argv[2] == "support":
            print("column")
            for j in support(y):
                print(j + 1)
        elif sys.argv[2] == "query":
            q = [F(float(rows[0]["q" + str(j)])) for j in range(len(names))]
            print("decision")
            print(euclidean_query(y, q, F(float(rows[0]["cutoff"]))))
        elif sys.argv[2] == "basis":
            vectors = basis(y)
            print(
                ",".join(
                    f"{end}{j}" for j in range(len(vectors)) for end in ("lo", "hi")
                )
            )
            for i in range(len(y)):
                print(
                    ",".join(f"{x:.17g}" for q in vectors for x in float_bounds(q[i]))
                )
        else:
            raise ValueError("unknown exact operation")

    if __name__ == "__main__":
        main()


def selection():
    """Exact selection at rational representative points, including equalities."""
    from fractions import Fraction as F
    import csv
    import sys

    def masks(model, points, eps, tau):
        pieces = {}
        for row in model:
            pieces.setdefault(int(row["piece"]), []).append(row)
        parsed = {}
        for j, rows in pieces.items():
            exact_counts = [F(r["count"]) for r in rows]
            if any(c.denominator != 1 or c <= 0 for c in exact_counts):
                raise ValueError("counts must be positive integers")
            counts = [int(c) for c in exact_counts]
            n = sum(counts)
            parsed[j] = (
                F(float(rows[0]["left"])),
                F(float(rows[0]["right"])),
                [F(float(r["b"])) for r in rows],
                [F(float(r["c"])) for r in rows],
                [F(c, n) for c in counts],
            )

        def risk(piece, t):
            _, _, b, c, p = piece
            return sum(w * (t * t + x * t + y) for w, x, y in zip(p, b, c))

        minima = []
        for piece in parsed.values():
            left, right, b, c, p = piece
            t = min(right, max(left, -sum(w * x for w, x in zip(p, b)) / 2))
            minima.append(risk(piece, t))
        minimum = min(minima)
        out = []
        for row in points:
            a, b = parsed[int(row["j"])], parsed[int(row["k"])]
            u, w = F(float(row["u"])), F(float(row["v"]))
            t, s = a[0] + u * (a[1] - a[0]), b[0] + w * (b[1] - b[0])
            d = [
                t * t + x * t + y - s * s - z * s - h
                for x, y, z, h in zip(a[2], a[3], b[2], b[3])
            ]
            delta = sum(p * x for p, x in zip(a[4], d))
            variance = sum(p * (x - delta) ** 2 for p, x in zip(a[4], d))
            mask = 1 if variance > 0 else 0
            if (
                mask
                and risk(a, t) - minimum <= eps
                and (delta >= 0 or n * delta**2 <= tau**2 * variance)
            ):
                mask |= 2
            out.append(mask)
        return out

    if __name__ == "__main__":
        with open(sys.argv[1], newline="") as f:
            model = list(csv.DictReader(f))
        with open(sys.argv[2], newline="") as f:
            points = list(csv.DictReader(f))
        values = masks(model, points, *(F(float(x)) for x in sys.argv[3:5]))
        print("mask")
        for value in values:
            print(value)


def circle():
    """Certified circle membership with rational interval arithmetic and true pi.

    No floating polynomial roots. Bernstein bounds cover every comparison arc.
    The arc containing q is handled by exact symbolic cancellation of (r-q)^2.
    Finite budgets return UNRESOLVED, never a tolerance-based acceptance.
    """
    from fractions import Fraction as F
    from math import comb
    import csv
    import sys

    class I:
        def __init__(self, lo=0, hi=None):
            self.lo = F(lo)
            self.hi = self.lo if hi is None else F(hi)
            if self.lo > self.hi:
                raise ValueError("reversed interval")

        def __add__(self, other):
            other = interval(other)
            return I(self.lo + other.lo, self.hi + other.hi)

        __radd__ = __add__

        def __neg__(self):
            return I(-self.hi, -self.lo)

        def __sub__(self, other):
            return self + -interval(other)

        def __rsub__(self, other):
            return interval(other) + -self

        def __mul__(self, other):
            other = interval(other)
            values = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
            return I(min(values), max(values))

        __rmul__ = __mul__

        def square(self):
            return I(
                0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                max(self.lo**2, self.hi**2),
            )

    def interval(x):
        return x if isinstance(x, I) else I(x)

    def pi_interval(bits=160):
        def atan_inverse(k):
            total, j = F(0), 0
            while True:
                term = F((-1) ** j, (2 * j + 1) * k ** (2 * j + 1))
                total += term
                remainder = F((-1) ** (j + 1), (2 * j + 3) * k ** (2 * j + 3))
                if abs(remainder) < F(1, 2**bits):
                    return I(
                        min(total, total + remainder), max(total, total + remainder)
                    )
                j += 1

        return 16 * atan_inverse(5) - 4 * atan_inverse(239)

    def polynomial_product(a, b):
        c = [I() for _ in range(len(a) + len(b) - 1)]
        for j, x in enumerate(a):
            for k, y in enumerate(b):
                c[j + k] = c[j + k] + x * y
        return c

    def bernstein_range(coeff, left, right):
        """Enclose a polynomial on an interval, including uncertain arc endpoints."""
        degree = len(coeff) - 1
        # Power coefficients under x=left+(right-left)*u, u in [0,1].
        powers_l, powers_w = [I(1)], [I(1)]
        for _ in range(degree):
            powers_l.append(powers_l[-1] * left)
            powers_w.append(powers_w[-1] * (right - left))
        power = [
            sum(
                (
                    coeff[j] * comb(j, k) * powers_l[j - k] * powers_w[k]
                    for j in range(k, degree + 1)
                ),
                I(),
            )
            for k in range(degree + 1)
        ]
        bernstein = [
            sum((power[k] * F(comb(j, k), comb(degree, k)) for k in range(j + 1)), I())
            for j in range(degree + 1)
        ]
        return I(min(x.lo for x in bernstein), max(x.hi for x in bernstein))

    def query(y, counts, q, cutoff, max_cells=20000, max_depth=40, bits=160):
        if not y or len(y) != len(counts) or any(c <= 0 for c in counts) or cutoff < 0:
            raise ValueError("invalid circle query")
        pi = pi_interval(bits)
        if any(abs(x) >= pi.lo for x in y + [q]):
            raise ValueError("signed angles in (-pi,pi) required")
        if len(set(y)) == 1:
            return ("IN" if q == y[0] else "OUT"), 0
        n = sum(counts)
        weights = [F(c, n) for c in counts]
        mean = lambda xs: sum((w * x for w, x in zip(weights, xs)), I())
        winding = [-1 if x >= 0 else 0 for x in y]
        lift = [I(x) + 2 * k * pi for x, k in zip(y, winding)]
        order = sorted(range(len(y)), key=lambda j: (y[j] < 0, y[j]))
        cuts = [I(y[j]) - pi if y[j] >= 0 else I(y[j]) + pi for j in order]
        arcs, left = [], -pi
        for j in range(len(y) + 1):
            right = cuts[j] if j < len(y) else pi
            arcs.append((left, right, lift[:]))
            if j < len(y):
                index = order[j]
                winding[index] += 1
                lift[index] = I(y[index]) + 2 * winding[index] * pi
            left = right
        containing = [j for j, (l, r, _) in enumerate(arcs) if l.hi < q < r.lo]
        if len(containing) != 1:
            return "UNRESOLVED", 0
        home = containing[0]
        qlift = arcs[home][2]
        lq = [(q - x).square() for x in qlift]
        visited, unresolved = 0, False
        for j, (left, right, lift) in enumerate(arcs):
            if j == home:
                mu = mean(lift)
                var = mean([(x - mu).square() for x in lift])
                margin = cutoff**2 * var - n * (q - mu).square()
                if margin.hi < 0:
                    return "OUT", visited
                if margin.lo < 0:
                    unresolved = True
                continue
            a = [x - z.square() for x, z in zip(lq, lift)]
            b = [2 * z for z in lift]
            ma, mb = mean(a), mean(b)
            ac, bc = [x - ma for x in a], [x - mb for x in b]
            p = [ma, mb, I(-1)]
            variance = [
                mean([x.square() for x in ac]),
                2 * mean([x * z for x, z in zip(ac, bc)]),
                mean([z.square() for z in bc]),
            ]
            pp = polynomial_product(p, p)
            k = [
                (cutoff**2 * variance[i] if i < 3 else I()) - n * pp[i]
                for i in range(5)
            ]
            stack = [(F(0), F(1), 0)]
            while stack:
                visited += 1
                if visited > max_cells:
                    return "UNRESOLVED", visited
                lo, hi, depth = stack.pop()
                l, r = (1 - lo) * left + lo * right, (1 - hi) * left + hi * right
                pr, kr = bernstein_range(p, l, r), bernstein_range(k, l, r)
                if pr.hi <= 0 or kr.lo >= 0:
                    continue
                mid = (lo + hi) / 2
                t = (1 - mid) * left + mid * right
                pm, km = bernstein_range(p, t, t), bernstein_range(k, t, t)
                if pm.lo > 0 and km.hi < 0:
                    return "OUT", visited
                if depth >= max_depth:
                    unresolved = True
                else:
                    stack.extend([(lo, mid, depth + 1), (mid, hi, depth + 1)])
        return ("UNRESOLVED" if unresolved else "IN"), visited

    if __name__ == "__main__":
        with open(sys.argv[1], newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        exact = lambda x: F(float(x))
        answer, cells = query(
            [exact(r["y"]) for r in rows],
            [int(r["count"]) for r in rows],
            exact(rows[0]["q"]),
            exact(rows[0]["cutoff"]),
        )
        print("decision,cells")
        print(f"{answer},{cells}")


def graph():
    """Exact full-edge inequality inversion; Python standard library only.

    CSV coefficients are integers; coordinates and cutoff denote exact binary64
    numbers. Sturm variations with one-sided endpoint signs count distinct roots
    in open intervals, including repeated roots. No tolerance or plotting mesh.
    """
    from fractions import Fraction as F

    def trim(a):
        while len(a) > 1 and a[-1] == 0:
            a.pop()
        return a

    def add(a, b):
        c = [F(0)] * max(len(a), len(b))
        for i, x in enumerate(a):
            c[i] += x
        for i, x in enumerate(b):
            c[i] += x
        return trim(c)

    def scale(a, s):
        return trim([x * s for x in a])

    def mul(a, b):
        c = [F(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] += x * y
        return trim(c)

    def deriv(a):
        return [i * a[i] for i in range(1, len(a))] or [F(0)]

    def evaluate(a, x):
        value = F(0)
        for y in reversed(a):
            value = value * x + y
        return value

    def remainder(a, b):
        a = a[:]
        while len(a) >= len(b) and a != [0]:
            k = len(a) - len(b)
            a = add(a, [F(0)] * k + scale(b, -a[-1] / b[-1]))
        return a

    def sturm(a):
        out = [trim(a[:]), deriv(a)]
        if out[1] == [0]:
            return out[:1]
        while True:
            r = scale(remainder(out[-2], out[-1]), -1)
            if r == [0]:
                return out
            # Positive normalization limits integer growth and preserves signs.
            out.append(scale(r, 1 / abs(r[-1])))

    def side_sign(a, x, side):
        order = 0
        while a != [0]:
            v = evaluate(a, x)
            if v:
                return (1 if v > 0 else -1) * (side**order)
            a = deriv(a)
            order += 1
        return 0

    def variations(chain, x, side):
        signs = [side_sign(a, x, side) for a in chain]
        signs = [s for s in signs if s]
        return sum(a != b for a, b in zip(signs, signs[1:]))

    def violation(p, k, left, right, budget=10000):
        """Existence of P>0 and K<0 on a CLOSED rational interval.

        Every interior component with the required signs borders a root of P*K
        or an interval endpoint. With <=1 interior root, the two inward endpoint
        signs therefore suffice. Otherwise bisect, counting roots exactly.
        """
        p, k = trim(p[:]), trim(k[:])
        if p == [0] or k == [0]:
            return False
        chain = sturm(mul(p, k))
        stack = [(left, right)]
        visited = 0
        while stack:
            visited += 1
            if visited > budget:
                raise RuntimeError("exact root isolation budget exceeded")
            lo, hi = stack.pop()
            for x, side in ((lo, 1), (hi, -1)):
                if evaluate(p, x) > 0 and evaluate(k, x) < 0:
                    return True
                if lo < hi and side_sign(p, x, side) > 0 and side_sign(k, x, side) < 0:
                    return True
            if lo == hi:
                continue
            roots = variations(chain, lo, 1) - variations(chain, hi, -1)
            if roots > 1:
                mid = (lo + hi) / 2
                stack.extend(((lo, mid), (mid, hi)))
        return False

    def query_prepared(prepared, cutoff):
        for piece, p, q, negative_mean_square, left, right in prepared:
            k = add(scale(q, cutoff * cutoff), negative_mean_square)
            if violation(p, k, left, right):
                return "OUT", piece
        return "IN", ""

    return mul, scale, query_prepared


def ranking():
    """Evaluation-only rational oracle for vertex-supported S4 laws.

    The R inference path never imports population weights or this oracle. Reuse
    the existing exact polynomial inversion, with genuinely rational coordinates.
    """
    import csv
    import sys
    from fractions import Fraction as F

    mul, scale, query_prepared = graph()

    def read(path):
        with open(path, newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    def oracle(rows):
        out = []
        for row in rows:
            w, a, b = (F(row[k]) for k in ("W", "A", "B"))
            if w <= 0:
                raise ValueError("positive total weight required")
            t = min(F(1), max(F(0), -b / (2 * w)))
            risk = a / w + b * t / w + t * t
            out.append((row, t, risk))
        best = min(x[2] for x in out)
        writer = csv.writer(sys.stdout)
        writer.writerow(
            [
                "edge",
                "u",
                "v",
                "numerator",
                "denominator",
                "risk_numerator",
                "risk_denominator",
                "minimum",
            ]
        )
        for row, t, risk in out:
            writer.writerow(
                [
                    row["edge"],
                    row["u"],
                    row["v"],
                    t.numerator,
                    t.denominator,
                    risk.numerator,
                    risk.denominator,
                    int(risk == best),
                ]
            )

    def prepare(edges, point):
        groups = {}
        for row in edges:
            groups.setdefault(row["edge"], []).append(row)
        first = next(iter(groups.values()))
        weights = [int(r["count"]) for r in first]
        n = sum(weights)
        if any(x <= 0 for x in weights):
            raise ValueError("positive observed counts required")
        weights = [F(x, n) for x in weights]
        t = F(int(point["numerator"]), int(point["denominator"]))
        if not 0 <= t <= 1:
            raise ValueError("query outside edge")
        query = groups[point["edge"]]
        losses = [t * t + F(r["b"]) * t + F(r["c"]) for r in query]
        mean = lambda x: sum(w * y for w, y in zip(weights, x))
        prepared = []
        for edge, records in groups.items():
            if [r["count"] for r in records] != [r["count"] for r in first]:
                raise ValueError("inconsistent observed counts")
            a = [x - F(r["c"]) for x, r in zip(losses, records)]
            b = [-F(r["b"]) for r in records]
            ma, mb = mean(a), mean(b)
            ac, bc = [x - ma for x in a], [x - mb for x in b]
            p = [ma, mb, F(-1)]
            variance = [
                mean([x * x for x in ac]),
                2 * mean([x * y for x, y in zip(ac, bc)]),
                mean([x * x for x in bc]),
            ]
            prepared.append((edge, p, variance, scale(mul(p, p), -n), F(0), F(1)))
        return prepared

    def evaluate(edges, points, cutoffs):
        # Prepare each true mean once; share it across G/L and stop a method as
        # soon as any mean is excluded. All 36 comparator edges remain present.
        if [r["method"] for r in cutoffs] != ["G", "L"]:
            raise ValueError("G/L cutoffs required")
        active = {r["method"]: F(float(r["cutoff"])) for r in cutoffs}
        if any(x < 0 for x in active.values()):
            raise ValueError("negative cutoff")
        result = {}
        for i, point in enumerate(points, 1):
            prepared = prepare(edges, point)
            for method, cutoff in list(active.items()):
                decision, edge = query_prepared(prepared, cutoff)
                if decision == "OUT":
                    result[method] = (decision, i, edge)
                    del active[method]
            if not active:
                break
        for method in active:
            result[method] = ("IN", "", "")
        writer = csv.writer(sys.stdout)
        writer.writerow(["method", "decision", "witness_point", "witness_edge"])
        for row in cutoffs:
            writer.writerow([row["method"], *result[row["method"]]])

    if __name__ == "__main__":
        if sys.argv[1] == "oracle":
            oracle(read(sys.argv[2]))
        elif sys.argv[1] == "evaluate":
            evaluate(read(sys.argv[2]), read(sys.argv[3]), read(sys.argv[4]))
        else:
            raise ValueError("unknown operation")


def cutoff():
    """Independent exact-rational audit of fixed-B ranks, nesting and decimal gap.

    The native interval engine certifies the all-domain draw enclosures. This audit
    checks their aggregation; it does not re-prove the continuous optimization.
    Exit 2 means a valid enclosure whose requested gap is unresolved, not failure.
    """
    import json
    import math
    import sys
    from fractions import Fraction

    def number(value, allow_infinite=False):
        # jsonlite serializes a valid unbounded enclosure as "Inf". Missing and
        # NaN values must never be confused with such an unresolved upper bound.
        if allow_infinite and value == "Inf":
            return math.inf
        assert isinstance(value, (int, float, str)) and not isinstance(value, bool)
        value = float(value)
        assert not math.isnan(value) and (
            math.isfinite(value) or (allow_infinite and value == math.inf)
        )
        return value

    def verify(proof):
        lo = [[number(x) for x in row] for row in proof["lower"]]
        hi = [[number(x, True) for x in row] for row in proof["upper"]]
        B, k = proof["B"], proof["k"]
        assert type(B) is int and B > 0 and type(k) is int
        alpha = Fraction.from_float(number(proof["alpha"]))
        assert 0 < alpha < 1
        assert k == B + 1 - (alpha * (B + 1)).__floor__() and k <= B
        assert len(lo) == len(hi) == B
        for a, b in zip(lo, hi):
            assert len(a) == len(b) == 2
            assert all(math.isfinite(x) and x >= 0 for x in a)
            assert all(not math.isnan(x) and x >= 0 for x in b)
            assert all(x <= y for x, y in zip(a, b))
            assert a[0] >= a[1] and b[0] >= b[1]
        target = Fraction(proof["requested_decimal_gap"])
        assert target >= 0
        cutoff = {
            key: [number(x, key != "lower") for x in proof["cutoff"][key]]
            for key in ("lower", "upper", "used", "gap")
        }
        assert all(len(value) == 2 for value in cutoff.values())
        attained = True
        for j, method in enumerate(("G", "L")):
            lower = sorted(a[j] for a in lo)[k - 1]
            upper = sorted(a[j] for a in hi)[k - 1]
            # R named vectors serialize as arrays; the ordering is G, L.
            assert lower == cutoff["lower"][j]
            assert upper == cutoff["upper"][j]
            assert upper == cutoff["used"][j]
            if math.isfinite(upper):
                width = Fraction.from_float(upper) - Fraction.from_float(lower)
                assert math.isfinite(cutoff["gap"][j])
                assert Fraction.from_float(cutoff["gap"][j]) >= width
                ok = width <= target
                print(
                    f"{method}: exact binary64 endpoint gap {float(width):.17g}; decimal target {target}; attained={ok}"
                )
            else:
                assert cutoff["gap"][j] == math.inf
                ok = False
                print(f"{method}: infinite upper bound; unresolved")
            attained &= ok
        return attained

    if __name__ == "__main__":
        with open(sys.argv[1], encoding="utf-8") as stream:
            proof = json.load(stream)
        sys.exit(0 if verify(proof) else 2)


def lengths():
    """Independent full-circle angular audit; all published summaries are degrees."""
    import csv, json, math, sys
    from fractions import Fraction as F
    from pathlib import Path
    from collections import defaultdict
    directory=Path(sys.argv[1]);reps=int(sys.argv[2]);tolerance=F(float.fromhex(sys.argv[3]))
    def atan(k):
        total=F(0)
        for j in range(200):
            total+=F((-1)**j,(2*j+1)*k**(2*j+1))
            remainder=F((-1)**(j+1),(2*j+3)*k**(2*j+3))
            if abs(remainder)<F(1,2**180):return min(total,total+remainder),max(total,total+remainder)
        raise ArithmeticError('pi enclosure did not converge')
    a,b=atan(5);c,d=atan(239);pi_lo,pi_hi=16*a-4*d,16*b-4*c
    def rational(row,key):return F(float.fromhex(row[key]))
    def outward(value,upper):
        x=float(value)
        if (upper and F(x)<value) or (not upper and F(x)>value):
            x=math.nextafter(x,math.inf if upper else -math.inf)
        return x
    def decimal_out(value,upper):
        v=value*10**4
        integer=-(-v.numerator//v.denominator) if upper else v.numerator//v.denominator
        return f'{integer//10**4}.{integer%10**4:04d}'
    with (directory/'raw.csv').open(newline='',encoding='utf-8') as stream:rows=list(csv.DictReader(stream))
    methods=('G','L','EL')
    expected={(f'circle-p{p}-n{n}-r{r:06d}',m) for p in ('-0.125','0','0.5','1')
        for n in (50,100,200) for r in range(1,reps+1) for m in methods}
    assert len(rows)==36*reps and {(x['id'],x['method']) for x in rows}==expected
    groups=defaultdict(list);pairs=defaultdict(dict);attained=0
    for row in rows:
        assert row['id']==f"circle-p{float(row['parameter']):.17g}-n{int(row['n'])}-r{int(row['replication']):06d}"
        lo,hi=rational(row,'radian_lower_hex'),rational(row,'radian_upper_hex')
        inside,outside,band,M=(int(row[k]) for k in ('inside','outside','band','grid_size'))
        assert M==2**24 and min(inside,outside,band)>=0 and inside+outside+band<=M
        assert 0<=lo<=2*pi_lo*inside/M and hi>=2*pi_hi*(M-outside)/M and lo<=hi
        dlo,dhi=F(360*inside,M),F(360*(M-outside),M)
        assert rational(row,'degrees_lower_hex')==dlo and rational(row,'degrees_upper_hex')==dhi
        assert F(float(row['degrees_lower']))==dlo and F(float(row['degrees_upper']))==dhi
        assert F(float(row['degrees_width']))==dhi-dlo
        declared=row['status']=='ANGLE_PRECISION_ATTAINED'
        assert declared==(dhi-dlo<=tolerance)
        assert 0<=dlo<=dhi<=360 and (row['method']!='EL' or band==0)
        attained+=declared
        groups[float(row['parameter']),int(row['n']),row['method']].append((dlo,dhi,declared))
        pairs[row['id']][row['method']]=(dlo,dhi)
    summaries=[]
    for (p,n,m),values in sorted(groups.items()):
        assert len(values)==reps
        lo=sum((v[0] for v in values),F(0))/reps;hi=sum((v[1] for v in values),F(0))/reps
        summaries.append(dict(parameter=p,n=n,method=m,replications=reps,
            mean_degrees_lower=outward(lo,False),mean_degrees_upper=outward(hi,True),
            mean_degrees_lower_exact=str(lo),mean_degrees_upper_exact=str(hi),
            mean_degrees_lower_tex=decimal_out(lo,False),mean_degrees_upper_tex=decimal_out(hi,True),
            precision_attained=sum(v[2] for v in values),valid_wide=sum(not v[2] for v in values)))
    diffs=defaultdict(list);by_id={r['id']:r for r in rows}
    for id,z in pairs.items():
        assert z['G'][0]>=z['L'][0] and z['G'][1]>=z['L'][1], ('G/L nesting',id)
        r=by_id[id]
        for m,k in (('G','L'),('G','EL'),('L','EL')):
            lo=z[m][0]-z[k][1];hi=z[m][1]-z[k][0]
            if (m,k)==('G','L'):lo=max(F(0),lo)
            diffs[float(r['parameter']),int(r['n']),m+'_minus_'+k].append((lo,hi))
    paired=[]
    for (p,n,contrast),values in sorted(diffs.items()):
        lo=sum((v[0] for v in values),F(0))/reps;hi=sum((v[1] for v in values),F(0))/reps
        paired.append(dict(parameter=p,n=n,replications=reps,contrast=contrast,
            mean_degrees_lower=outward(lo,False),mean_degrees_upper=outward(hi,True),
            mean_degrees_lower_exact=str(lo),mean_degrees_upper_exact=str(hi)))
    assert len(summaries)==36 and len(paired)==36
    for name,data in [('summary.csv',summaries),('paired.csv',paired)]:
        with (directory/name).open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    report=dict(status='PASS',datasets=12*reps,rows=len(rows),precision_attained=attained,
        valid_wide=len(rows)-attained,summary_groups=len(summaries),paired_groups=len(paired),
        display_units='degrees',tolerance_degrees_binary64=sys.argv[3],
        checks=['complete G/L/EL grid','true pi via rational Machin bounds','integer partition measures',
            'exact angular width comparison','G/L nesting','rational means and paired differences'],
        scope='Numerical angular enclosures conditional on frozen observations, tuning and B draws; no population coverage or Monte Carlo CI claim.')
    (directory/'audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'Circle angles: {len(rows)} rows, {attained} precision attained; exact arithmetic PASS')


def angle_table():
    """Bind the published degree table, including TeX, to exact mean bounds."""
    import csv,json,sys
    from pathlib import Path
    from fractions import Fraction as F
    directory=Path(sys.argv[1])
    def read(name):
        with (directory/name).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
    summary=read('summary.csv');table=read('table.csv');text=(directory/'table.tex').read_text(encoding='utf-8')
    assert len(summary)==36 and len(table)==12
    grouped={(float(x['parameter']),int(x['n']),x['method']):x for x in summary}
    assert len(grouped)==36
    seen=set()
    for row in table:
        key=float(row['lambda']),int(row['n']);assert key not in seen;seen.add(key)
        values=[]
        for method,column in [('G','Global'),('L','Local'),('EL','EL')]:
            z=grouped[key+(method,)];assert int(z['replications'])==int(row['replications'])
            lo,hi=z['mean_degrees_lower_tex'],z['mean_degrees_upper_tex']
            assert F(lo)<=F(z['mean_degrees_lower_exact'])<=F(z['mean_degrees_upper_exact'])<=F(hi)
            value=f'[{lo}, {hi}]';assert row[column]==value;values.append(value)
        line=f'{key[0]:g} & {key[1]} & '+' & '.join(values)+r' \\'
        assert text.splitlines().count(line)==1, ('TeX row mismatch',key)
    if int(table[0]['replications'])!=1000:assert 'PILOT' in text
    audit=json.loads((directory/'audit.json').read_text(encoding='utf-8'))
    audit['checks'].append('published CSV and TeX table match outward rational means')
    (directory/'audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
    print('Angular CSV/TeX table: independent exact audit PASS')


if __name__ == "__main__":
    import sys

    command = sys.argv.pop(1)
    handlers = {
        "index": index,
        "euclidean": euclidean,
        "selection": selection,
        "circle": circle,
        "ranking": ranking,
        "cutoff": cutoff,
        "angles": lengths,
        "angle-table": angle_table,
    }
    handlers[command]()
