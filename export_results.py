#!/usr/bin/env python3
"""Render the bundled reference figure and tables without R.

Bundled CSVs are final aggregates, not observations or a new simulation.  The
--run option instead reads the same aggregates from a completed production run.
"""
import argparse
import csv
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import tempfile

# The read-only entry point must not leave an import cache in the checkout.
_bytecode = sys.dont_write_bytecode
try:
    sys.dont_write_bytecode = True
    from verify_release import csv_fingerprint
finally:
    sys.dont_write_bytecode = _bytecode

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT / "reference"
SAMPLES = (50, 100, 200)
PARAMETERS = (Fraction(-1, 8), Fraction(0), Fraction(1, 2), Fraction(1))
COVERAGE_FIELDS = ["family", "parameter", "n", "geometry", "method", "replications",
                   "coverage_lower", "coverage_upper", "MC95_lower", "MC95_upper",
                   "membership_unresolved"]
RANKING_FIELDS = ["geometry", "n", "replications", "G_total", "L_total",
                  "G_only_total", "L_only_total", "L_strict_subset_G_count"]
CIRCLE_FIELDS = ["parameter", "n", "method", "replications",
                 "mean_degrees_lower_exact", "mean_degrees_upper_exact"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def normalized(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_csv(path, fields):
    require(path.is_file(), f"Missing CSV: {path.name}")
    before = path.stat()
    require(before.st_size < 2_000_000, f"Oversized CSV: {path.name}")
    data = normalized(path)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
            f"CSV changed while reading: {path.name}")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames),
            f"Missing or duplicate CSV columns: {path.name}")
    require(reader.fieldnames == fields, f"Unexpected CSV columns: {path.name}")
    rows = list(reader)
    require(rows and all(None not in row and all(v is not None for v in row.values()) for row in rows),
            f"Malformed CSV rows: {path.name}")
    return rows


def number(value):
    require(isinstance(value, str) and re.fullmatch(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?(?:/\d+)?", value),
            f"Invalid finite numeric value: {value!r}")
    return Fraction(value)


def integer(value):
    x = number(value)
    require(x.denominator == 1, f"Noninteger value: {value}")
    return x.numerator


def rounded(value, places):
    """Exact round-to-nearest, ties-to-even, with fixed decimal places."""
    scale = 10 ** places
    negative = value < 0
    q, r = divmod(abs(value.numerator) * scale, value.denominator)
    if 2 * r > value.denominator or (2 * r == value.denominator and q % 2):
        q += 1
    sign = "-" if negative and q else ""
    if places == 0:
        return f"{sign}{q}"
    return f"{sign}{q // scale}.{q % scale:0{places}d}"


def complete(rows, key, expected, name):
    keys = [key(row) for row in rows]
    require(len(keys) == len(set(keys)), f"Duplicate {name} condition")
    require(set(keys) == expected, f"Missing or unexpected {name} conditions")
    require(all(integer(row["replications"]) == 1000 for row in rows),
            f"{name} requires 1000 completed replications per condition")


def coverage_data(rows):
    expected = set()
    for n in SAMPLES:
        for p in (Fraction(0), Fraction(1, 2), Fraction(1)):
            expected.update(("euclidean", p, n, "euclidean", m) for m in ("G", "WALD_CHI"))
        for g in ("euclidean", "kendall", "footrule", "graph"):
            expected.update(("graph" if g == "graph" else "ranking", Fraction(0), n, g, m) for m in ("G", "L"))
        for p in PARAMETERS:
            expected.update(("circle", p, n, "circle", m) for m in ("G", "L", "EL_S1_RAW", "WALD_S1_KDE_CENTERED"))
    key = lambda r: (r["family"], number(r["parameter"]), integer(r["n"]), r["geometry"], r["method"])
    complete(rows, key, expected, "coverage")
    for row in rows:
        lo, hi, mc_lo, mc_hi = [number(row[k]) for k in COVERAGE_FIELDS[6:10]]
        require(0 <= mc_lo <= lo <= hi <= mc_hi <= 1, "Invalid coverage or Monte Carlo enclosure")
        require(lo == hi and integer(row["membership_unresolved"]) == 0,
                "Reference results require resolved coverage membership")
        require((lo * 1000).denominator == 1, "Coverage is incompatible with 1000 datasets")
    return {key(r): r for r in rows}


def ranking_data(rows):
    key = lambda r: (r["geometry"], integer(r["n"]))
    complete(rows, key, {(g, n) for g in ("euclidean", "kendall", "footrule") for n in SAMPLES}, "ranking")
    for row in rows:
        g, local, difference, local_only, reduced = [integer(row[k]) for k in RANKING_FIELDS[3:]]
        require(0 <= local <= g <= 24_000 and difference == g - local and local_only == 0,
                "Invalid resolved ranking cardinalities or nesting")
        require(0 <= reduced <= 1000 and reduced <= difference <= 24 * reduced,
                "Invalid ranking strict-subset count")
    return {key(r): r for r in rows}


def angle_display(lo, hi):
    require(0 <= lo <= hi <= 360, "Invalid exact angular enclosure")
    display = rounded((lo + hi) / 2, 1)
    error = max(abs(number(display) - lo), abs(hi - number(display)))
    require(error < Fraction(3, 10), "Displayed angular mean cannot certify error < 0.3 degrees")
    return display


def circle_data(rows):
    key = lambda r: (number(r["parameter"]), integer(r["n"]), r["method"])
    complete(rows, key, {(p, n, m) for p in PARAMETERS for n in SAMPLES for m in ("G", "L", "EL")}, "circle")
    for row in rows:
        angle_display(number(row[CIRCLE_FIELDS[-2]]), number(row[CIRCLE_FIELDS[-1]]))
    return {key(r): r for r in rows}


def select_coverage(rows):
    return [{k: row[k] for k in COVERAGE_FIELDS} for row in rows
            if not (row["family"] == "euclidean" and row["method"] in ("L", "WALD_F_GAUSSIAN_REFERENCE"))]


def source_ranking(rows):
    reduced = []
    for row in rows:
        out = {k: row[k] for k in RANKING_FIELDS[:3]}
        for prefix in ("G", "L", "G_only", "L_only"):
            lo, hi = number(row[prefix + "_mean_lower"]), number(row[prefix + "_mean_upper"])
            require(lo == hi and integer(row[prefix + "_count_unresolved"]) == 0,
                    "Unresolved ranking cardinality")
            total = integer(rounded(lo * 1000, 0))
            # R writes binary64 mean values as decimal; recover their integer sum.
            require(abs(lo * 1000 - total) < Fraction(1, 10**9), "Ranking mean is not an integer-count average")
            out[prefix + "_total"] = str(total)
        require(row["L_strict_subset_G_certified"] == row["L_strict_subset_G_possible"],
                "Unresolved ranking strict-subset count")
        out["L_strict_subset_G_count"] = row["L_strict_subset_G_certified"]
        reduced.append(out)
    return reduced


def source_circle(rows):
    for row in rows:
        require(integer(row["precision_attained"]) == 1000 and integer(row["valid_wide"]) == 0,
                "Incomplete Circle precision certification")
    return [{k: row[k] for k in CIRCLE_FIELDS} for row in rows]


def source_inputs(run, manifest):
    require(run.is_dir(), "Completed run directory does not exist")
    paths = {"coverage_figure_data": run / "figures" / "figure-data.csv",
             "ranking_summary": run.with_name(run.name + "-ranking-candidates") / "summary.csv",
             "angle_summary": run / "angles" / "summary.csv"}
    sources = {}
    require(set(manifest["source_exports"]) == set(paths), "Unexpected reference source export set")
    for name, spec in manifest["source_exports"].items():
        sources[name] = []
        require(csv_fingerprint(paths[name], spec["keys"], spec["excluded_execution_columns"], sources[name]) == spec["scientific_fingerprint"],
                f"Run scientific aggregate differs from the frozen reference: {spec['path']}")
    return (select_coverage(sources["coverage_figure_data"]),
            source_ranking(sources["ranking_summary"]), source_circle(sources["angle_summary"]))


def load_bundle(bundle=BUNDLE):
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    require(set(manifest) == {"schema", "scope", "files", "mapping", "source_exports", "rendered_lf_sha256"},
            "Unexpected reference manifest fields")
    require(manifest["schema"] == "frechet-reference-results-v1", "Unknown reference manifest")
    require(set(manifest["source_exports"]) == {"coverage_figure_data", "ranking_summary", "angle_summary"},
            "Unexpected reference source export set")
    for spec in manifest["source_exports"].values():
        require(set(spec) == {"keys", "excluded_execution_columns", "location", "path", "scientific_fingerprint"},
                "Unexpected reference source fields")
    require(set(manifest["files"]) == {"coverage.csv", "ranking.csv", "circle.csv",
            "figure-coverage.template.tex", "table-ranking.template.tex", "table-circle.template.tex", "plot-style.tex"},
            "Unexpected reference bundle file set")
    for name, expected in manifest["files"].items():
        require(digest(normalized(bundle / name)) == expected, f"Reference bundle hash mismatch: {name}")
    rows = (read_csv(bundle / "coverage.csv", COVERAGE_FIELDS),
            read_csv(bundle / "ranking.csv", RANKING_FIELDS),
            read_csv(bundle / "circle.csv", CIRCLE_FIELDS))
    return manifest, rows


def render(inputs, bundle=BUNDLE):
    cov, ranks, angles = coverage_data(inputs[0]), ranking_data(inputs[1]), circle_data(inputs[2])
    figure = (bundle / "figure-coverage.template.tex").read_text(encoding="utf-8")
    panels = [("euclidean", p, "euclidean") for p in (Fraction(0), Fraction(1, 2), Fraction(1))]
    panels += [("graph" if g == "graph" else "ranking", Fraction(0), g) for g in ("euclidean", "kendall", "footrule", "graph")]
    panels += [("circle", p, "circle") for p in PARAMETERS] * 2
    styles = {"paperg": "G", "paperl": "L", "paperel": "EL_S1_RAW", "paperwald": "WALD_S1_KDE_CENTERED"}

    def coverage_block(match):
        panel, style = int(match[1]), match[2]
        family, parameter, geometry = panels[panel - 1]
        method = "WALD_CHI" if family == "euclidean" and style == "paperwald" else styles[style]
        lines = []
        for n in SAMPLES:
            row = cov[(family, parameter, n, geometry, method)]
            if 8 <= panel <= 11 and 100 * number(row["MC95_upper"]) < 88:
                continue
            if panel >= 12 and 100 * number(row["MC95_lower"]) > 88:
                continue
            offset = {"G": -.09, "L": -.03, "EL_S1_RAW": .03, "WALD_S1_KDE_CENTERED": .09}[method] if family == "circle" else (-.0325 if method == "G" else .0325)
            # Match the numerical renderer's binary64 plotting coordinates and
            # six-place decimal serialization; exact statistical values stay in CSV.
            y = float(row["coverage_lower"])
            values = (n * math.pow(2, offset), 100 * y,
                      100 * (y - float(row["MC95_lower"])),
                      100 * (float(row["MC95_upper"]) - y))
            lines.append(" ".join(f"{v:.6f}" for v in values) + "\\\\")
        require(lines, "Current layout cannot represent this changed coverage series")
        return "\n".join(lines)

    figure, count = re.subn(r"@@COVERAGE_(\d+)_(\w+)@@", coverage_block, figure)
    require(count == 32, "Invalid coverage template")
    require(len(re.findall(r"^\d+\.\d{6} \d+\.\d{6} \d+\.\d{6} \d+\.\d{6}\\\\$", figure, re.M)) == 92,
            "Unexpected plotted coverage row count")
    ranking_lines = []
    for geometry, label in (("euclidean", "Euclidean"), ("kendall", "Kendall"), ("footrule", "footrule")):
        for n in SAMPLES:
            row = ranks[(geometry, n)]
            values = [rounded(Fraction(integer(row[k]), 1000), 3) for k in ("G_total", "L_total", "G_only_total")]
            values.append(rounded(Fraction(integer(row["L_strict_subset_G_count"]), 10), 1))
            ranking_lines.append(f"{label if n == 50 else ''} & {n} & " + " & ".join(values) + " \\\\")
        ranking_lines.append("\\hline")
    circle_lines = []
    for parameter, label in zip(PARAMETERS, ("$-1/8$", "$0$", "$1/2$", "$1$")):
        for n in SAMPLES:
            values = []
            for method in ("G", "L", "EL"):
                row = angles[(parameter, n, method)]
                values.append(angle_display(number(row[CIRCLE_FIELDS[-2]]), number(row[CIRCLE_FIELDS[-1]])))
            circle_lines.append(f"{label if n == 50 else ''} & {n} & " + " & ".join(values) + " \\\\")
        circle_lines.append("\\hline")
    outputs = {"figure-coverage.tex": figure,
               "table-ranking.tex": (bundle / "table-ranking.template.tex").read_text(encoding="utf-8").replace("@@ROWS@@", "\n".join(ranking_lines)),
               "table-circle.tex": (bundle / "table-circle.template.tex").read_text(encoding="utf-8").replace("@@ROWS@@", "\n".join(circle_lines)),
               "plot-style.tex": (bundle / "plot-style.tex").read_text(encoding="utf-8")}
    require(all("@@" not in text for text in outputs.values()), "Unfilled reference template marker")
    return {name: text.encode("utf-8") for name, text in outputs.items()}


def write_outputs(output, outputs):
    output = output.resolve()
    require(output != ROOT and output not in ROOT.parents, "Output may not replace the source checkout")
    if ROOT in output.parents:
        require(output.relative_to(ROOT).parts[0] == "build", "Within source checkout use build/ for generated outputs")
    if output.exists():
        require(output.is_dir() and {p.name for p in output.iterdir()} == set(outputs), "Choose a new or identical results output directory")
        require(all(not (output / name).is_symlink() and (output / name).read_bytes() == data for name, data in outputs.items()),
                "Existing reference outputs differ; choose a new directory")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".results-pending-", dir=output.parent) as pending:
        temporary = Path(pending) / "results"
        temporary.mkdir()
        for name, data in outputs.items():
            (temporary / name).write_bytes(data)
        temporary.rename(output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check", help="Check bundled aggregates and exact reference render hashes")
    check.add_argument("--run", type=Path, help="Also verify a completed run's aggregate render against the reference")
    check.add_argument("--output", type=Path, help="Also verify existing exported results without writing")
    export = sub.add_parser("export", help="Export reference TeX fragments and figure style")
    export.add_argument("--run", type=Path, help="Read aggregates from a completed 1000-replication run instead of the bundle")
    export.add_argument("--output", type=Path, default=ROOT / "build" / "results")
    args = parser.parse_args(argv)
    manifest, bundled = load_bundle()
    reference = render(bundled)
    require({name: digest(data) for name, data in reference.items()} == manifest["rendered_lf_sha256"], "Rendered outputs differ from the fixed reference")
    if args.command == "check":
        if args.run is not None:
            require(render(source_inputs(args.run.resolve(), manifest)) == reference, "Run aggregates do not reproduce the fixed reference")
        if args.output is not None:
            directory = args.output.resolve()
            require(directory.is_dir() and {p.name for p in directory.iterdir()} == set(reference), "Unexpected exported results files")
            require(all(not (directory / name).is_symlink() and (directory / name).read_bytes() == data for name, data in reference.items()),
                    "Existing export differs from the fixed reference")
        print(json.dumps({"status": "PASS", "coverage_input_rows": 90, "coverage_plotted_rows": 92,
                          "ranking_rows": 9, "circle_rows": 36,
                          "angle_error_degrees": "strictly less than 0.3 for all 36 displayed means",
                          "reference_checked": True}, sort_keys=True))
    else:
        outputs = render(source_inputs(args.run.resolve(), manifest)) if args.run is not None else reference
        require(outputs == reference, "Run aggregates do not reproduce the fixed reference; no files written")
        write_outputs(args.output, outputs)
        print(json.dumps({"status": "PASS", "input": "completed run" if args.run else "bundled final aggregates",
                          "matches_reference": outputs == reference,
                          "files": {name: digest(data) for name, data in outputs.items()}}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, UnicodeError, ZeroDivisionError) as error:
        print(f"export_results.py: {error}", file=sys.stderr)
        sys.exit(1)
