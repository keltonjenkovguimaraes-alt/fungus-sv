#!/usr/bin/env python3
"""SV Report Card Generator with SMaHT confidence mapping."""

# SMaHT (Zhang et al. 2025) confidence designation mapping:
# TRIPLE_TRIANGULATED ≈ SMaHT HighConf (multi-layer, cross-validation)
# DOUBLE_CONFIRMED ≈ SMaHT LowConf (two layers agree)
# SINGLE_LINE ≈ SMaHT LowConf (single evidence line)
# WEAK / CONTRADICTED ≈ SMaHT LikelyArtifact (no cross-evidence)


def generate_report_card(result):
    """Generate per-SV report card."""
    return f"""
=== VALIDATION REPORT: {result.sv_id} ===
Type: {result.sv_type}
T-Score: {result.t_score:.3f}
Confidence: {result.tier.name}
"""


def generate_summary_table(results):
    """Generate summary with SMaHT-style confidence tiers."""
    summary = "SV Validation Summary\n"
    summary += "=" * 50 + "\n"
    summary += f"{'ID':<30} {'Type':<8} {'T-Score':<10} {'Confidence':<15}\n"
    summary += "-" * 50 + "\n"
    
    # Sort by T-score descending
    sorted_results = sorted(results, key=lambda r: r.t_score, reverse=True)
    for r in sorted_results:
        summary += f"{r.sv_id:<30} {r.sv_type:<8} {r.t_score:<10.3f} {r.tier.name:<15}\n"
    
    # Add size-stratified summary
    summary += "\n" + generate_size_stratified_summary(results)
    
    return summary


def generate_size_stratified_summary(results):
    """Size-stratified summary (Liu et al. 2024: SVs cluster 50-400bp)."""
    if not results:
        return ""
    
    bins = {
        '50-100 bp': [],
        '100-500 bp': [],
        '500-5000 bp': [],
        '>5000 bp': []
    }
    
    for r in results:
        size = getattr(r, 'sv_size', 0)
        if size < 100:
            bins['50-100 bp'].append(r)
        elif size < 500:
            bins['100-500 bp'].append(r)
        elif size < 5000:
            bins['500-5000 bp'].append(r)
        else:
            bins['>5000 bp'].append(r)
    
    lines = ["Size-Stratified Summary:"]
    lines.append(f"  {'Bin':<15} {'N':>5} {'Mean T':>8} {'HighConf':>10} {'LowConf':>10}")
    lines.append(f"  {'-'*50}")
    
    for bin_name, svs in bins.items():
        if not svs:
            continue
        n = len(svs)
        mean_t = sum(r.t_score for r in svs) / n
        high = sum(1 for r in svs if r.t_score >= 0.6)
        low = n - high
        lines.append(f"  {bin_name:<15} {n:>5} {mean_t:>8.3f} {high:>10} {low:>10}")
    
    lines.append("\n  SMaHT-style confidence (Zhang et al. 2025):")
    lines.append("    HighConf ≈ T ≥ 0.6, LowConf ≈ T < 0.6, LikelyArtifact ≈ T < 0.2")
    
    return "\n".join(lines)
