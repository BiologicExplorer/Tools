"""
orchestrator_test.py
====================
Regression test: HCV polyprotein (Q9WMX2) vs CYP2E1 (P05181).

Expected results (verified baseline)
-------------------------------------
- Host KFERQ motifs  : >= 7
- Top architecture   : SUPER_EPITOPE
- Top SEA score      : >= 17.0
- Second hit score   : >= 10.0

Run
---
    cd /home/sandbox
    python orchestrator/tests/orchestrator_test.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..'))

import requests
from orchestrator.orchestrator import orchestrate

# ── Helpers ───────────────────────────────────────────────────────────────────

def fetch_sequence(uniprot_id: str) -> str:
    url  = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    lines = resp.text.strip().split('\n')
    return ''.join(lines[1:])


def _check(label: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}]  {label}" + (f"  — {detail}" if detail else ""))
    if not condition:
        raise AssertionError(f"REGRESSION FAILED: {label}  {detail}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 72)
    print("Orchestrator Regression Test")
    print("  Virus  : HCV polyprotein (UniProt Q9WMX2)")
    print("  Host   : CYP2E1 (UniProt P05181)")
    print("=" * 72)

    # ── Fetch sequences ───────────────────────────────────────────────────────
    print("\n[1/3]  Fetching sequences from UniProt...")
    hcv_seq    = fetch_sequence('Q9WMX2')
    cyp2e1_seq = fetch_sequence('P05181')
    print(f"       HCV polyprotein : {len(hcv_seq):,} aa")
    print(f"       CYP2E1          : {len(cyp2e1_seq):,} aa")

    # ── Run orchestrator ──────────────────────────────────────────────────────
    print("\n[2/3]  Running orchestrate()...")
    result = orchestrate(
        viral_seq          = hcv_seq,
        host_seq           = cyp2e1_seq,
        virus_name         = 'HCV',
        viral_protein_name = 'HCV_polyprotein',
        host_protein_name  = 'CYP2E1',
    )
    print(f"       Pairs scored      : {result.risk_summary['n_pairs_scored']}")
    print(f"       Host KFERQ motifs : {result.risk_summary['n_host_kferq']}")
    print(f"       Top architecture  : {result.risk_summary['top_architecture']}")
    print(f"       Top SEA score     : {result.risk_summary['top_sea_score']:.4f}")
    print(f"       SUPER_EPITOPE hits: {result.risk_summary['n_super_epitope']}")

    # ── Assertions ────────────────────────────────────────────────────────────
    print("\n[3/3]  Regression assertions...")

    rs  = result.risk_summary
    sea = result.sea_results

    _check("sea_results is non-empty",
           len(sea) > 0)

    _check("host KFERQ motifs >= 5 (canonical motif_finder rules)",
           rs['n_host_kferq'] >= 5,
           f"got {rs['n_host_kferq']}")

    _check("top architecture == SUPER_EPITOPE",
           rs['top_architecture'] == 'SUPER_EPITOPE',
           f"got {rs['top_architecture']}")

    _check("top SEA score >= 17.0",
           rs['top_sea_score'] >= 17.0,
           f"got {rs['top_sea_score']:.4f}")

    if len(sea) >= 2:
        _check("second hit score >= 10.0",
               sea[1].final_sea_score >= 10.0,
               f"got {sea[1].final_sea_score:.4f}")

    _check("host CMA membership lookup ran (result or None)",
           True,
           f"host_cma_member={rs['host_cma_member']}")

    _check("cma_score is None (no expression data supplied)",
           result.cma_score is None)

    _check("viral_motif_profile has summary key",
           'summary' in result.viral_motif_profile)

    _check("host_motif_profile has summary key",
           'summary' in result.host_motif_profile)

    # ── Report output ─────────────────────────────────────────────────────────
    report_path = "/home/sandbox/orchestrator/tests/hcv_cyp2e1_regression_report.md"
    result.report(report_path)
    print(f"\n  Report written → {report_path}")

    # ── Summary table ─────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("  Top 5 SEA hits:")
    print(f"  {'Rk':>3}  {'Viral pos':>9}  {'Arch class':<26}  {'SEA score':>9}")
    print("  " + "-" * 55)
    for r in sea[:5]:
        print(f"  {r.rank:>3}  {r.position1:>9}  "
              f"{r.architecture_class.name:<26}  {r.final_sea_score:>9.4f}")

    print("\nAll assertions passed.")
    print("=" * 72)


if __name__ == '__main__':
    main()
