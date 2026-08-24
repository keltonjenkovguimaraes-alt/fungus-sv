#!/usr/bin/env python3
"""
Corrected synthetic-benchmark evaluator for FUNGUS-SV.

Fixes the bug in valid_sv/benchmarks/run_calibration.py's evaluate():
that function matches each consensus SV to the FIRST truth SV within
2kb of the same type, without removing matched truth SVs from the
pool -- so two consensus calls near the same truth SV could both be
counted as true positives, and per-type breakdowns computed after the
fact (by filtering a separately-run pass) could disagree with a pooled
total, because filtering-then-counting and counting-then-summing are
not guaranteed to agree once double-matching is possible.

This script instead:
  1. Performs ONE-TO-ONE greedy matching (nearest distance first,
     same type, <=2kb), removing each truth SV from the pool the
     moment it's matched, so no truth SV can produce two TPs.
  2. Builds a single list of match records (each tagged with its SV
     type) as the one and only source of truth.
  3. Derives BOTH the pooled stats and the per-type stats by filtering
     that same list -- so pooled and per-type are mathematically
     guaranteed to reconcile (sum of per-type tp/fn always equals
     pooled tp/fn), which the original ad-hoc table9 was not.

Usage (single replicate):
    python3 evaluate_synthetic_replicates.py \
        --consensus-vcf results/rep1/consensus_svs.vcf \
        --truth-vcf results/benchmarks/synthetic_truth_rep1.vcf \
        --validation-json results/rep1/validation/validation_results.json

Usage (pool multiple replicates):
    python3 evaluate_synthetic_replicates.py \
        --consensus-vcf results/rep1/consensus_svs.vcf results/rep2/consensus_svs.vcf results/rep3/consensus_svs.vcf \
        --truth-vcf results/benchmarks/synthetic_truth_rep1.vcf results/benchmarks/synthetic_truth_rep2.vcf results/benchmarks/synthetic_truth_rep3.vcf \
        --validation-json results/rep1/validation/validation_results.json results/rep2/validation/validation_results.json results/rep3/validation/validation_results.json \
        --output results/synthetic_benchmark_reconciled.json
"""

import argparse
import json
import re
from collections import defaultdict


def parse_truth(truth_vcf):
    truth = {}
    with open(truth_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            m = re.search(r"SVTYPE=(\w+)", parts[7])
            svtype = m.group(1) if m else "UNK"
            truth[parts[2]] = {"chrom": parts[0], "pos": int(parts[1]), "type": svtype}
    return truth


def parse_consensus(consensus_vcf):
    consensus = []
    with open(consensus_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            m = re.search(r"SVTYPE=(\w+)", parts[7])
            svtype = m.group(1) if m else "UNK"
            consensus.append({"id": parts[2], "chrom": parts[0], "pos": int(parts[1]), "type": svtype})
    return consensus


def match_one_to_one(consensus, truth, max_dist=2000):
    """
    Greedy nearest-distance, one-to-one matching. Each truth SV can be
    claimed by at most one consensus call. Returns:
      matches: list of {consensus_id, truth_id, type, distance}
      unmatched_consensus: consensus calls that matched no truth (candidate FPs)
    """
    # Build all candidate (distance, consensus_idx, truth_id) triples within range/type/chrom
    candidates = []
    for ci, c in enumerate(consensus):
        for tid, t in truth.items():
            if t["chrom"] == c["chrom"] and t["type"] == c["type"]:
                dist = abs(t["pos"] - c["pos"])
                if dist <= max_dist:
                    candidates.append((dist, ci, tid))

    # Greedy: claim nearest pairs first, each truth and each consensus used at most once
    candidates.sort(key=lambda x: x[0])
    used_truth, used_consensus = set(), set()
    matches = []
    for dist, ci, tid in candidates:
        if tid in used_truth or ci in used_consensus:
            continue
        used_truth.add(tid)
        used_consensus.add(ci)
        matches.append({
            "consensus_id": consensus[ci]["id"],
            "truth_id": tid,
            "type": truth[tid]["type"],
            "distance": dist,
        })

    unmatched_consensus = [c for ci, c in enumerate(consensus) if ci not in used_consensus]
    return matches, unmatched_consensus


def evaluate_replicate(consensus_vcf, truth_vcf, validation_json=None, max_dist=2000):
    truth = parse_truth(truth_vcf)
    consensus = parse_consensus(consensus_vcf)
    matches, unmatched_consensus = match_one_to_one(consensus, truth, max_dist)

    tscores = {}
    if validation_json:
        with open(validation_json) as f:
            vdata = json.load(f)
        tscores = {r["sv_id"]: r["t_score"] for r in vdata.get("results", [])}

    return {
        "truth": truth,
        "consensus": consensus,
        "matches": matches,
        "unmatched_consensus": unmatched_consensus,
        "tscores": tscores,
    }


def summarize(replicates):
    """
    replicates: list of evaluate_replicate() outputs.
    Returns pooled stats and per-type stats derived from the SAME
    underlying match lists, so they are guaranteed to reconcile.
    """
    all_truth_by_type = defaultdict(int)
    all_matches_by_type = defaultdict(int)
    all_detected = 0
    all_fp = 0

    for rep in replicates:
        for t in rep["truth"].values():
            all_truth_by_type[t["type"]] += 1
        for m in rep["matches"]:
            all_matches_by_type[m["type"]] += 1
        all_detected += len(rep["consensus"])
        all_fp += len(rep["unmatched_consensus"])

    types = sorted(all_truth_by_type.keys())
    per_type = {}
    for t in types:
        truth_n = all_truth_by_type[t]
        tp = all_matches_by_type[t]
        fn = truth_n - tp
        per_type[t] = {
            "truth_n": truth_n, "tp": tp, "fn": fn,
            "recall": round(tp / truth_n, 4) if truth_n else 0.0,
        }

    pooled_truth_n = sum(all_truth_by_type.values())
    pooled_tp = sum(all_matches_by_type.values())
    pooled_fn = pooled_truth_n - pooled_tp
    pooled_detected = all_detected
    pooled_fp = all_fp

    pooled = {
        "truth_n": pooled_truth_n, "detected_n": pooled_detected,
        "tp": pooled_tp, "fp": pooled_fp, "fn": pooled_fn,
        "recall": round(pooled_tp / pooled_truth_n, 4) if pooled_truth_n else 0.0,
        "precision": round(pooled_tp / pooled_detected, 4) if pooled_detected else 0.0,
    }

    # Self-check: this should ALWAYS pass by construction, since both
    # pooled and per-type are filtered from the same match list.
    assert sum(v["tp"] for v in per_type.values()) == pooled_tp, \
        "Internal inconsistency -- this should be mathematically impossible; report as a bug."
    assert sum(v["fn"] for v in per_type.values()) == pooled_fn, \
        "Internal inconsistency -- this should be mathematically impossible; report as a bug."

    return {"pooled": pooled, "per_type": per_type}


def main():
    ap = argparse.ArgumentParser(description="Reconciled synthetic-benchmark evaluator (pooled + per-type, guaranteed consistent)")
    ap.add_argument("--consensus-vcf", nargs="+", required=True)
    ap.add_argument("--truth-vcf", nargs="+", required=True)
    ap.add_argument("--validation-json", nargs="+", default=None)
    ap.add_argument("--max-dist", type=int, default=2000)
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    n = len(args.consensus_vcf)
    assert len(args.truth_vcf) == n, "Number of --truth-vcf must match --consensus-vcf"
    val_jsons = args.validation_json if args.validation_json else [None] * n
    assert len(val_jsons) == n, "Number of --validation-json must match --consensus-vcf, if given"

    replicates = []
    for i in range(n):
        rep = evaluate_replicate(args.consensus_vcf[i], args.truth_vcf[i], val_jsons[i], args.max_dist)
        replicates.append(rep)
        print(f"Replicate {i+1}: {len(rep['truth'])} truth SVs, {len(rep['consensus'])} consensus calls, "
              f"{len(rep['matches'])} matched (one-to-one)")

    summary = summarize(replicates)

    print(f"\n{'='*60}")
    print("  POOLED (from one-to-one matching, no double-counting)")
    print(f"{'='*60}")
    p = summary["pooled"]
    print(f"  truth_n={p['truth_n']} detected_n={p['detected_n']} tp={p['tp']} fp={p['fp']} fn={p['fn']}")
    print(f"  recall={p['recall']*100:.1f}%  precision={p['precision']*100:.1f}%")

    print(f"\n  PER-TYPE (derived from the SAME matches -- guaranteed to sum to pooled above)")
    for t, v in summary["per_type"].items():
        print(f"  {t}: truth_n={v['truth_n']} tp={v['tp']} fn={v['fn']} recall={v['recall']*100:.1f}%")
    print(f"{'='*60}\n")

    if args.output:
        with open(args.output, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"Written: {args.output}")


if __name__ == "__main__":
    main()
