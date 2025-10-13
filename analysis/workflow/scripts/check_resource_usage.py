#!/usr/bin/env python3
import re
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
from datetime import datetime

# Usage
# -----
# My personal favorite
#     python workflow/scripts/check_resource_usage.py --logs log | tee -a bench.tsv | column -t -s $'\t' | less
# Without appending to an existing file
#     python workflow/scripts/check_resource_usage.py --logs log | column -t -s $'\t' | less
# If your benchmark paths in the log are relative to a project root (default: ".")
#     python workflow/scripts/check_resource_usage.py --logs log --root /my/project/root
# Multiple logs + export to file
#     python workflow/scripts/check_resource_usage.py --logs log_a.txt log_b.txt --out bench.tsv


# ---------- Bench TSV reader ----------
def read_benchmark_metrics(bench_path: Path) -> Tuple[Optional[float], Optional[float]]:
    """
    Return (peak_rss_mb, elapsed_min) from a Snakemake benchmark TSV.
    - Tries several common column names.
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
                peak_rss_mb = v
                break

    elapsed_min = (elapsed_s / 60.0) if elapsed_s is not None else None
    return peak_rss_mb, elapsed_min


# ---------- Log parser ----------
RULE_LINE = re.compile(r"^rule\s+([A-Za-z0-9_./-]+):\s*$")
ERROR_IN_RULE_LINE = re.compile(r"^Error in rule\s+([A-Za-z0-9_./-]+):\s*$")
WILDCARDS_LINE = re.compile(r"^\s*wildcards:\s*(.+?)\s*$")
RESOURCES_LINE = re.compile(r"^\s*resources:\s*(.+?)\s*$")
BENCHMARK_LINE = re.compile(r"^\s*benchmark:\s*([^\s].*?)\s*$")
LOG_LINE = re.compile(r"^\s*log:\s*(.+?)\s*$")
SUBMIT_LINE = re.compile(r"Job\s+\d+\s+has been submitted with SLURM jobid\s+\d+\s*\(log:\s*(.*?)\)", re.IGNORECASE)


def _clean_log_part(part: str) -> str:
    return part.split(" (", 1)[0].strip()


def parse_resources(s: str) -> Dict[str, str]:
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
    {rule, wildcards, alloc_mem_mb, alloc_runtime_min, benchmark_file, job_logs[]}
    """
    rows = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        curr_rule = None
        curr_wc = ""
        curr_res_s = ""
        curr_bench = ""
        curr_logs: List[str] = []

        for line in f:
            m_rule = RULE_LINE.match(line)
            m_err = ERROR_IN_RULE_LINE.match(line)
            if m_rule or m_err:
                if curr_rule and curr_bench:
                    res = parse_resources(curr_res_s) if curr_res_s else {}
                    rows.append({
                        "rule": curr_rule,
                        "wildcards": curr_wc,
                        "alloc_mem_mb": _as_float_or_none(res.get("mem_mb")),
                        "alloc_runtime_min": _as_float_or_none(res.get("runtime")),
                        "benchmark_file": curr_bench,
                        "job_logs": curr_logs.copy(),
                    })
                curr_rule = (m_rule.group(1) if m_rule else m_err.group(1))
                curr_wc = ""
                curr_res_s = ""
                curr_bench = ""
                curr_logs = []
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

            m_log = LOG_LINE.match(line)
            if m_log:
                parts = [_clean_log_part(p) for p in m_log.group(1).split(",") if p.strip()]
                curr_logs.extend(parts)
                continue

            for sm in SUBMIT_LINE.finditer(line):
                path_part = sm.group(1).strip()
                for p in [x.strip() for x in path_part.split(",") if x.strip()]:
                    curr_logs.append(_clean_log_part(p))

        if curr_rule and curr_bench:
            res = parse_resources(curr_res_s) if curr_res_s else {}
            rows.append({
                "rule": curr_rule,
                "wildcards": curr_wc,
                "alloc_mem_mb": _as_float_or_none(res.get("mem_mb")),
                "alloc_runtime_min": _as_float_or_none(res.get("runtime")),
                "benchmark_file": curr_bench,
                "job_logs": curr_logs.copy(),
            })
    return rows


def _as_float_or_none(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(str(x).strip())
    except Exception:
        return None


def determine_fail_reason(log_paths: List[Path]) -> str:
    """
    For failed jobs (bench file missing), check logs for cause.
    """
    if not log_paths:
        return "other"
    for lp in log_paths:
        if not lp.exists():
            continue
        try:
            text = lp.read_text(errors="ignore")
        except Exception:
            continue
        if "DUE TO TIME LIMIT" in text:
            return "time"
        if (" Killed " in text) or ("Some of the step tasks have been OOM Killed." in text):
            return "memory"
    return "other"


def extract_jobid_node_submit(log_paths: List[Path]) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Extract SLURM job ID from filename, node name from 'host:' line,
    and submission time (normalized to ISO) from the first '[...]' timestamp.
    """
    jobid = None
    node = None
    submit_time = None
    for lp in log_paths:
        # Jobid from filename
        m = re.search(r"(\d+)\.log$", str(lp))
        if m:
            jobid = m.group(1)
        # Parse contents
        if lp.exists():
            try:
                with lp.open("r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        # Submission timestamp
                        if line.startswith("[") and line.endswith("]\n") and submit_time is None:
                            raw = line.strip("[]\n")
                            try:
                                dt = datetime.strptime(raw, "%a %b %d %H:%M:%S %Y")
                                submit_time = dt.strftime("%Y-%m-%d %H:%M:%S")
                            except Exception:
                                submit_time = raw  # fallback
                        # Node
                        if line.lower().startswith("host:"):
                            node = line.split(":", 1)[1].strip()
            except Exception:
                continue
    return jobid, node, submit_time


# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="Compare Snakemake log-declared resources vs benchmarked usage.")
    ap.add_argument("--logs", nargs="+", required=True, help="Paths to Snakemake text logs to parse")
    ap.add_argument("--root", default=".", help="Optional project root to resolve relative benchmark/log paths")
    ap.add_argument("--out", default="/dev/stdout", help="Optional TSV output path")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    entries: List[Dict] = []
    for lp in args.logs:
        entries.extend(parse_log(Path(lp)))

    if not entries:
        print("No jobs with benchmark lines were found in the provided logs.")
        return

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

        job_log_paths: List[Path] = []
        for p in e.get("job_logs", []):
            if not p:
                continue
            pp = Path(p)
            if not pp.is_absolute():
                pp = (root / pp).resolve()
            job_log_paths.append(pp)

        peak_rss_mb, elapsed_min = read_benchmark_metrics(bench_path)

        mem_util = mem_over = None
        if e["alloc_mem_mb"] is not None and peak_rss_mb is not None:
            mem_util = 100.0 * peak_rss_mb / e["alloc_mem_mb"]
            mem_over = peak_rss_mb > e["alloc_mem_mb"]

        time_util = time_over = None
        if e["alloc_runtime_min"] is not None and elapsed_min is not None:
            time_util = 100.0 * elapsed_min / e["alloc_runtime_min"]
            time_over = elapsed_min > e["alloc_runtime_min"]

        # Success/failure detection
        if bench_path.exists():
            fail_reason = "success"
        else:
            fail_reason = determine_fail_reason(job_log_paths)

        jobid, node, submit_time = extract_jobid_node_submit(job_log_paths)

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
            "fail_reason": fail_reason,
            "jobid": jobid,
            "node": node,
            "submit_time": submit_time,
        })

    df = pd.DataFrame(rows).sort_values(["rule", "benchmark_file"]).reset_index(drop=True)

    if args.out:
        df.to_csv(args.out, index=False, sep="\t")


if __name__ == "__main__":
    main()
