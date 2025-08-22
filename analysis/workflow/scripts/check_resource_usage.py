#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd

# # basic
# python parse_logs_vs_bench.py --logs /path/to/snakemake_run.log

# # if your benchmark paths in the log are relative to a project root:
# python parse_logs_vs_bench.py --logs /path/to/snakemake_run.log --root /my/project/root

# # multiple logs + export CSV
# python parse_logs_vs_bench.py --logs log_a.txt log_b.txt --out bench_compare.csv


# ---------- Bench TSV reader ----------
def read_benchmark_metrics(bench_path: Path) -> Tuple[Optional[float], Optional[float]]:
    """
    Return (peak_rss_mb, elapsed_min) from a Snakemake benchmark TSV.
    - Tries several common column names.
    - Assumes RSS is KB if the value is large (converts to MB).
    - If multiple rows exist (retries), uses the max per metric.
    """
    if not bench_path.exists():
        return None, None
    try:
        df = pd.read_csv(bench_path, sep="\t")
    except Exception:
        try:
            df = pd.read_csv(bench_path, sep="\t", engine="python")
        except Exception:
            return None, None

    df.columns = [c.strip().lower() for c in df.columns]

    # elapsed time (seconds) -> minutes
    elapsed_s = None
    for cand in ("s", "seconds", "elapsed", "elapsed_s"):
        if cand in df.columns:
            v = pd.to_numeric(df[cand], errors="coerce").max()
            if pd.notna(v):
                elapsed_s = float(v)
                break

    # peak RSS (prefer KB and convert)
    peak_rss_mb = None
    for cand in ("max_rss", "peak_rss", "rss", "maxresidentset", "max_resident_set"):
        if cand in df.columns:
            v = pd.to_numeric(df[cand], errors="coerce").max()
            if pd.notna(v):
                v = float(v)
                peak_rss_mb = v / 1024.0 if v > 8192 else v  # heuristic: large => KB
                break

    elapsed_min = (elapsed_s / 60.0) if elapsed_s is not None else None
    return peak_rss_mb, elapsed_min

# ---------- Log parser ----------
RULE_BLOCK_START = re.compile(r"^\[[^\]]+\]\s*?\n?")  # lines like: [Fri Aug 22 12:32:09 2025]
RULE_LINE = re.compile(r"^rule\s+([A-Za-z0-9_./-]+):\s*$")
WILDCARDS_LINE = re.compile(r"^\s*wildcards:\s*(.+?)\s*$")
RESOURCES_LINE = re.compile(r"^\s*resources:\s*(.+?)\s*$")
BENCHMARK_LINE = re.compile(r"^\s*benchmark:\s*([^\s].*?)\s*$")

def parse_resources(s: str) -> Dict[str, str]:
    """
    Parse 'mem_mb=1191, mem_mib=1136, runtime=10, ...' -> dict
    Handles key=value items split by commas. Values may include spaces or '=' in rare cases, we split on first '='.
    """
    out = {}
    for part in re.split(r",\s*", s.strip()):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
        else:
            out[part.strip()] = ""
    return out

def parse_log(path: Path) -> List[Dict]:
    """
    Parse a Snakemake text log and return a list of entries:
    {rule, wildcards_str, mem_mb, runtime, benchmark_file}
    """
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        curr_rule = None
        curr_wc = ""
        curr_res_s = ""
        curr_bench = ""

        for line in f:
            # detect start of a new block by 'rule ...:' line (we don't strictly need timestamp boundary)
            m_rule = RULE_LINE.match(line)
            if m_rule:
                # flush previous (if benchmark/resources present)
                if curr_rule and curr_bench:
                    res = parse_resources(curr_res_s) if curr_res_s else {}
                    rows.append({
                        "rule": curr_rule,
                        "wildcards": curr_wc,
                        "alloc_mem_mb": _as_float_or_none(res.get("mem_mb")),
                        "alloc_runtime_min": _as_float_or_none(res.get("runtime")),
                        "benchmark_file": curr_bench
                    })
                # start new block
                curr_rule = m_rule.group(1)
                curr_wc = ""
                curr_res_s = ""
                curr_bench = ""
                continue

            m_wc = WILDCARDS_LINE.match(line)
            if m_wc:
                curr_wc = m_wc.group(1)
                continue

            m_res = RESOURCES_LINE.match(line)
            if m_res:
                curr_res_s = m_res.group(1)
                continue

            m_bench = BENCHMARK_LINE.match(line)
            if m_bench:
                curr_bench = m_bench.group(1)
                continue

        # flush last block
        if curr_rule and curr_bench:
            res = parse_resources(curr_res_s) if curr_res_s else {}
            rows.append({
                "rule": curr_rule,
                "wildcards": curr_wc,
                "alloc_mem_mb": _as_float_or_none(res.get("mem_mb")),
                "alloc_runtime_min": _as_float_or_none(res.get("runtime")),
                "benchmark_file": curr_bench
            })

    return rows

def _as_float_or_none(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(str(x).strip())
    except Exception:
        return None

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Compare Snakemake log-declared resources vs benchmarked usage.")
    ap.add_argument("--logs", nargs="+", required=True, help="Paths to Snakemake text logs to parse")
    ap.add_argument("--root", default=".", help="Optional project root to resolve relative benchmark paths")
    ap.add_argument("--out", default="/dev/stdout", help="Optional TSV output path")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    entries: List[Dict] = []
    for lp in args.logs:
        entries.extend(parse_log(Path(lp)))

    if not entries:
        print("No jobs with benchmark lines were found in the provided logs.")
        return

    # Deduplicate by (rule, wildcards, benchmark_file)
    seen = set()
    uniq = []
    for e in entries:
        key = (e["rule"], e["wildcards"], e["benchmark_file"])
        if key not in seen:
            seen.add(key)
            uniq.append(e)

    rows = []
    for e in uniq:
        bench_path = (root / e["benchmark_file"]).resolve() if not Path(e["benchmark_file"]).is_absolute() else Path(e["benchmark_file"])
        peak_rss_mb, elapsed_min = read_benchmark_metrics(bench_path)

        mem_util = None
        mem_over = None
        if e["alloc_mem_mb"] is not None and peak_rss_mb is not None:
            mem_util = 100.0 * peak_rss_mb / e["alloc_mem_mb"]
            mem_over = peak_rss_mb > e["alloc_mem_mb"]

        time_util = None
        time_over = None
        if e["alloc_runtime_min"] is not None and elapsed_min is not None:
            time_util = 100.0 * elapsed_min / e["alloc_runtime_min"]
            time_over = elapsed_min > e["alloc_runtime_min"]

        rows.append({
            "rule": e["rule"],
            "wildcards": e["wildcards"],
            "benchmark_file": str(bench_path),
            "alloc_mem_mb": e["alloc_mem_mb"],
            "peak_rss_mb": None if peak_rss_mb is None else round(peak_rss_mb, 1),
            "mem_util_%": None if mem_util is None else round(mem_util, 1),
            "mem_over": mem_over,
            "alloc_runtime_min": e["alloc_runtime_min"],
            "elapsed_min": None if elapsed_min is None else round(elapsed_min, 1),
            "time_util_%": None if time_util is None else round(time_util, 1),
            "time_over": time_over,
        })

    df = pd.DataFrame(rows).sort_values(["rule", "benchmark_file"]).reset_index(drop=True)

    if args.out:
        df.to_csv(args.out, index=False, sep="\t")
        print(f"\nWrote CSV: {args.out}")

if __name__ == "__main__":
    main()
