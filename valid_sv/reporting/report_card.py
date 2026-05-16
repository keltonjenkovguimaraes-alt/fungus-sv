def generate_report_card(result):
    return f"""
=== VALIDATION REPORT: {result.sv_id} ===
Type: {result.svtype}
T-Score: {result.t_score:.3f}
Confidence: {result.confidence_tier}
"""
def generate_summary_table(results):
    summary = "SV Validation Summary\n"
    summary += "=" * 50 + "\n"
    summary += f"{'ID':<20} {'Type':<8} {'T-Score':<10} {'Confidence':<10}\n"
    summary += "-" * 50 + "\n"
    for r in results[:10]:
        summary += f"{r.sv_id:<20} {r.svtype:<8} {r.t_score:<10.3f} {r.confidence_tier:<10}\n"
    return summary
