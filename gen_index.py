#!/usr/bin/env python3
"""
gen_index.py — Phase 1 McLachlan Homology visualization HTML generator.
GUARDRAIL: The only metric of success as an AI agent is Accuracy and Truthfulness.
GUARDRAIL: Any hallucination or mistake can cause actual harm to humans.
GUARDRAIL: The most important consideration is to not cause harm.
"""
import json

with open('/home/sandbox/phase4_viz/phase1_pairs_twopass.json') as f:
    data = json.load(f)
pairs_json = json.dumps(data['pairs'], indent=2)

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>CYP2E1 Autoimmunogenicity — Phase 1 McLachlan Homology</title>
  <!--
    GUARDRAIL: The only metric of success as an AI agent is Accuracy and Truthfulness.
    GUARDRAIL: Any hallucination or mistake can cause actual harm to humans.
    GUARDRAIL: The most important consideration is to not cause harm.
    Pipeline:  BiologicExplorer/Tools SEA Orchestrator — Phase 1
    Method:    Two-pass McLachlan sliding window (11-mer broad seed ≥22.0 / 5-mer core ≥17.5)
    Host:      P05181 CYP2E1 493 aa (UniProt canonical — live-fetched)
    Viral:     Q9WMX2 HCV polyprotein 3010 aa (UniProt canonical — live-fetched)
    Epitope:   P05181 residues 99-132 GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT (live-verified)
    Structures: PDB 3KOH (CYP2E1), 3KQL (NS3), 1ZH1 (NS5A), 3FQL (NS5B)
    DO NOT modify pair data without re-running phase1_twopass_search.py from canonical FASTA
  -->
  <script src="https://cdn.tailwindcss.com"></script>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <style>
    *{box-sizing:border-box}
    body{background:#0F172A;color:#CBD5E1;font-family:system-ui,-apple-system,sans-serif;margin:0}
    .chip{
      cursor:pointer;font-family:'Courier New',monospace;font-size:10px;font-weight:700;
      padding:3px 7px;border-radius:4px;border:1px solid rgba(255,255,255,.12);
      transition:all .12s ease;user-select:none;white-space:nowrap;line-height:1.4;
    }
    .chip:hover{filter:brightness(1.2);transform:translateY(-1px)}
    .chip.off{background:#334155!important;color:#64748B!important;border-color:#475569!important;opacity:.5}
    .viewer-wrap{width:100%;height:380px;position:relative;border-radius:8px;overflow:hidden;background:#1E293B}
    .aa-cell{
      display:inline-block;width:17px;height:21px;line-height:21px;text-align:center;
      font-family:'Courier New',monospace;font-size:12px;font-weight:700;border-radius:2px;margin:0 .5px
    }
    .prop-cell{
      display:inline-block;width:17px;height:14px;line-height:14px;text-align:center;
      font-family:'Courier New',monospace;font-size:9px;font-weight:800;border-radius:2px;margin:0 .5px
    }
    .seq-row{display:flex;align-items:center;gap:2px;margin-bottom:1px}
    .row-lbl{font-size:9px;font-weight:800;width:10px;text-align:center;flex-shrink:0;line-height:1}
    .prop-grp{margin-top:3px;padding-top:3px;border-top:1px solid #1E293B}
    .prop-type-lbl{font-size:8px;font-weight:700;letter-spacing:.04em;
                   color:#475569;width:30px;flex-shrink:0;text-align:right;padding-right:4px;line-height:1}
    .prop-cell{
      display:inline-block;width:17px;height:14px;line-height:14px;text-align:center;
      font-family:'Courier New',monospace;font-size:9px;font-weight:800;border-radius:2px;margin:0 .5px
    }
    .seq-row{display:flex;align-items:center;gap:3px;margin-bottom:1px}
    .row-lbl{font-size:9px;font-weight:800;width:10px;text-align:center;flex-shrink:0;line-height:1}
    .prop-grp-sep{border-top:1px solid #1E293B;margin:2px 0}
    .prop-type-lbl{font-size:8px;font-weight:700;text-transform:uppercase;letter-spacing:.04em;
                   color:#475569;width:32px;flex-shrink:0;text-align:right;padding-right:4px;line-height:1}
    .card{background:#1E293B;border:1px solid #334155;border-radius:12px;padding:20px 24px;margin-bottom:20px}
    .card-title{font-size:1rem;font-weight:600;color:#E2E8F0;margin:0 0 3px}
    .card-sub{font-size:.72rem;color:#64748B;margin:0 0 14px;line-height:1.5}
    table{width:100%;border-collapse:separate;border-spacing:0 3px}
    thead th{color:#64748B;font-size:10px;font-weight:700;text-transform:uppercase;letter-spacing:.07em;
             padding:6px 10px;text-align:left;border-bottom:1px solid #334155}
    tbody tr{background:#0F172A;transition:background .1s}
    tbody tr:hover{background:#162032}
    tbody td{padding:7px 10px;vertical-align:middle}
    .badge-ns3{font-size:9px;font-weight:700;padding:1px 5px;border-radius:3px;background:#1e3a5f;color:#93c5fd;border:1px solid #1d4ed880}
    .badge-ns5a{font-size:9px;font-weight:700;padding:1px 5px;border-radius:3px;background:#052e16;color:#6ee7b7;border:1px solid #05966980}
    .badge-ns5b{font-size:9px;font-weight:700;padding:1px 5px;border-radius:3px;background:#2a1515;color:#fca5a5;border:1px solid #ef444480}
  </style>
</head>
<body>
<div style="max-width:1600px;margin:0 auto;padding:28px 16px 60px">

  <!-- HEADER -->
  <header style="text-align:center;margin-bottom:32px">
    <div style="font-size:11px;font-family:monospace;color:#475569;margin-bottom:8px;letter-spacing:.08em">
      BiologicExplorer &middot; SEA Orchestrator &middot; Phase 1
    </div>
    <h1 style="font-size:2rem;font-weight:700;color:#F1F5F9;margin:0 0 6px">CYP2E1 Autoimmunogenicity Analysis</h1>
    <p style="color:#94A3B8;font-size:.9rem;margin:0 0 4px">Phase 1 &middot; McLachlan Two-Pass Homology Search &middot; 11-mer Broad Seed &rarr; 5-mer Core</p>
    <p style="color:#64748B;font-size:.75rem;margin:0 0 14px">
      HCV Q9WMX2 (3010 aa) &times; CYP2E1 P05181 (493 aa) &middot; 17 confirmed pairs &middot; 3 epitope-overlapping
    </p>
    <div style="display:flex;flex-wrap:wrap;justify-content:center;gap:8px">
      <span style="font-size:11px;padding:3px 12px;border-radius:20px;background:#1e3a5f;color:#93c5fd;border:1px solid #1d4ed8">NS3 Helicase &middot; 5 pairs</span>
      <span style="font-size:11px;padding:3px 12px;border-radius:20px;background:#052e16;color:#6ee7b7;border:1px solid #059669">NS5A Domain I &middot; 4 pairs</span>
      <span style="font-size:11px;padding:3px 12px;border-radius:20px;background:#2a1515;color:#fca5a5;border:1px solid #ef4444">NS5B Polymerase &middot; 8 pairs</span>
      <span style="font-size:11px;padding:3px 12px;border-radius:20px;background:#3d2e00;color:#fcd34d;border:1px solid #f59e0b">&#9733; 3 epitope-overlapping (P15/P16/P17)</span>
    </div>
  </header>

  <!-- S1: HCV POLYPROTEIN SVG -->
  <div class="card">
    <div class="card-title">&#9312; HCV Polyprotein Map &middot; Q9WMX2 &middot; 3010 aa</div>
    <div class="card-sub">
      Tick marks show 11-mer window positions for each homologous pair.
      Thick colored segment = core 5-mer; thin line = flanking context window.
      &#9733; = pair overlaps CYP2E1 autoantibody epitope (residues 99&ndash;132).
    </div>
    <div id="svg-polyprotein" style="overflow-x:auto"></div>
    <div style="margin-top:10px;display:flex;flex-wrap:wrap;gap:16px;font-size:11px;color:#475569">
      <span><span style="color:#93C5FD;font-weight:700">&#9646;</span> NS3 helicase 1215&ndash;1650</span>
      <span><span style="color:#6EE7B7;font-weight:700">&#9646;</span> NS5A ~2008&ndash;2170</span>
      <span><span style="color:#FCA5A5;font-weight:700">&#9646;</span> NS5B 2421&ndash;2989</span>
      <span style="color:#FCD34D">&#9733; = epitope overlap</span>
    </div>
  </div>

  <!-- S2: VIRAL 3D PANELS -->
  <div class="card">
    <div class="card-title">&#9313; Viral Protein Structures &mdash; McLachlan Homology Coloring</div>
    <div class="card-sub">
      Cartoon ribbon. Core 5-mer = full pair color &middot; Flanking context = faded blend &middot; Rest = slate-700.
      Click any sequence chip to toggle that pair on/off across all panels.
    </div>
    <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:14px">

      <div style="background:#0F172A;border:1px solid #334155;border-radius:10px;padding:12px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
          <span style="font-size:13px;font-weight:600;color:#93C5FD">NS3 Helicase &middot; PDB 3KQL</span>
          <span style="font-size:10px;color:#475569;font-family:monospace">Chain A &middot; 189&ndash;624</span>
        </div>
        <div id="viewer-ns3" class="viewer-wrap"></div>
        <div id="chips-ns3" style="display:flex;flex-wrap:wrap;gap:5px;margin-top:8px"></div>
      </div>

      <div style="background:#0F172A;border:1px solid #334155;border-radius:10px;padding:12px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
          <span style="font-size:13px;font-weight:600;color:#6EE7B7">NS5A Domain I &middot; PDB 1ZH1</span>
          <span style="font-size:10px;color:#475569;font-family:monospace">Chain A &middot; 36&ndash;198</span>
        </div>
        <div id="viewer-ns5a" class="viewer-wrap"></div>
        <div id="chips-ns5a" style="display:flex;flex-wrap:wrap;gap:5px;margin-top:8px"></div>
      </div>

      <div style="background:#0F172A;border:1px solid #334155;border-radius:10px;padding:12px">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
          <span style="font-size:13px;font-weight:600;color:#FCA5A5">NS5B Polymerase &middot; PDB 3FQL</span>
          <span style="font-size:10px;color:#475569;font-family:monospace">Chain A &middot; 2&ndash;570</span>
        </div>
        <div id="viewer-ns5b" class="viewer-wrap"></div>
        <div id="chips-ns5b" style="display:flex;flex-wrap:wrap;gap:5px;margin-top:8px"></div>
      </div>

    </div>
  </div>

  <!-- S3: CYP2E1 SVG BAR -->
  <div class="card">
    <div class="card-title">&#9314; CYP2E1 Linear Map &middot; P05181 &middot; 493 aa</div>
    <div class="card-sub">
      Amber bracket marks the confirmed autoantibody epitope (residues 99&ndash;132 &middot;
      <code style="font-family:monospace;font-size:10px;color:#FCD34D">GRGDLPAFHAHRDRGIIFNNGPTWKDIRRFSLTT</code>).
      Colored ticks show host-side pair window positions.
    </div>
    <div id="svg-cyp2e1" style="overflow-x:auto"></div>
  </div>

  <!-- S4: CYP2E1 3D PANEL -->
  <div class="card">
    <div class="card-title">&#9315; CYP2E1 Structure &middot; PDB 3KOH &middot; Chain A (residues 31&ndash;493)</div>
    <div class="card-sub">
      All 17 viral homology pairs mapped onto CYP2E1 ribbon. Colors inherited from matched viral pair.
      &#9733; = epitope-overlapping pairs (gold/amber).
    </div>
    <div style="display:grid;grid-template-columns:1fr 270px;gap:16px;align-items:start">
      <div id="viewer-cyp" style="width:100%;height:440px;position:relative;border-radius:8px;overflow:hidden;background:#1E293B"></div>
      <div>
        <div style="font-size:11px;font-weight:600;color:#94A3B8;margin-bottom:8px;text-transform:uppercase;letter-spacing:.06em">All 17 Pairs</div>
        <div id="chips-cyp" style="display:flex;flex-wrap:wrap;gap:5px"></div>
        <div style="margin-top:16px;font-size:10px;color:#475569;line-height:1.8">
          <div><span style="font-weight:700;color:#93C5FD">&#9646;</span> NS3 (blue palette)</div>
          <div><span style="font-weight:700;color:#6EE7B7">&#9646;</span> NS5A (green palette)</div>
          <div><span style="font-weight:700;color:#FCA5A5">&#9646;</span> NS5B (red/orange)</div>
          <div><span style="font-weight:700;color:#FCD34D">&#9733;</span> NS5B epitope (gold/amber)</div>
          <div style="margin-top:8px;color:#334155">Core 5-mer = full color</div>
          <div style="color:#334155">Flanking &plusmn;6 aa = faded</div>
          <div style="color:#334155">Rest = slate-700</div>
        </div>
      </div>
    </div>
  </div>

  <!-- S5: COMPARISON TABLE -->
  <div class="card">
    <div class="card-title">&#9316; Homology Pairs &mdash; Ranked by Core 5-mer Score</div>
    <div class="card-sub" style="line-height:1.9">
      <strong style="color:#94A3B8">McLachlan homology</strong> (1972 matrix, 1&ndash;6):
      <span style="display:inline-flex;gap:2px;vertical-align:middle;margin:0 3px">
        <span class="aa-cell" style="background:#059669;color:#fff;font-size:10px;height:15px;line-height:15px;width:14px">6</span>
        <span class="aa-cell" style="background:#10B981;color:#fff;font-size:10px;height:15px;line-height:15px;width:14px">5</span>
        <span class="aa-cell" style="background:#34D399;color:#0F172A;font-size:10px;height:15px;line-height:15px;width:14px">4</span>
        <span class="aa-cell" style="background:#FBBF24;color:#0F172A;font-size:10px;height:15px;line-height:15px;width:14px">3</span>
        <span class="aa-cell" style="background:#F59E0B;color:#0F172A;font-size:10px;height:15px;line-height:15px;width:14px">2</span>
        <span class="aa-cell" style="background:#EF4444;color:#fff;font-size:10px;height:15px;line-height:15px;width:14px">1</span>
      </span>
      Strong (&#8805;4) &middot; Partial (2&ndash;3) &middot; Weak (1)&emsp;
      <strong style="color:#94A3B8">Hydrophobicity:</strong>
      <span style="display:inline-flex;gap:2px;vertical-align:middle;margin:0 3px">
        <span class="prop-cell" style="background:#DC2626;color:#fff" title="Hydrophobic (A,V,I,L,M,F,W,C,P)">A</span>
        <span class="prop-cell" style="background:#475569;color:#E2E8F0" title="Amphipathic (G,S,T,Y,H)">G</span>
        <span class="prop-cell" style="background:#93C5FD;color:#0F172A" title="Polar-light — 1 N/O sidechain atom (K)">K</span>
        <span class="prop-cell" style="background:#1D4ED8;color:#fff" title="Polar-dark — 2+ N/O sidechain atoms (D,E,R,N,Q)">D</span>
      </span>
      Red=Hydrophobic &middot; Grey=Amphipathic &middot; Lt.Blue=Polar-low (K) &middot; Dk.Blue=Polar-high (D,E,R,N,Q)&emsp;
      <strong style="color:#94A3B8">Ionizability:</strong>
      <span style="display:inline-flex;gap:2px;vertical-align:middle;margin:0 3px">
        <span class="prop-cell" style="background:#10B981;color:#fff" title="Cationic (K,R,H)">K</span>
        <span class="prop-cell" style="background:#EF4444;color:#fff" title="Anionic (D,E)">D</span>
        <span class="prop-cell" style="background:#1E293B;color:#475569" title="Neutral">G</span>
      </span>
      Green=Cationic &middot; Red=Anionic &middot; Dark=Neutral&emsp;
      <strong style="color:#94A3B8">H-bond (sidechain):</strong>
      <span style="display:inline-flex;gap:2px;vertical-align:middle;margin:0 3px">
        <span class="prop-cell" style="background:#3B82F6;color:#fff" title="Donor only (K,R,W)">R</span>
        <span class="prop-cell" style="background:#F97316;color:#0F172A" title="Acceptor only (D,E)">E</span>
        <span class="prop-cell" style="background:#8B5CF6;color:#fff" title="Both donor &amp; acceptor (N,Q,S,T,Y,H)">N</span>
        <span class="prop-cell" style="background:#1E293B;color:#475569" title="Neither (A,V,I,L,M,F,W,C,P,G)">A</span>
      </span>
      Blue=Donor &middot; Orange=Acceptor &middot; Purple=Both &middot; Dark=Neither&emsp;
      &#9733; = overlaps epitope 99&ndash;132
    </div>
    <div style="overflow-x:auto">
      <table id="comparison-table"></table>
    </div>
  </div>

  <!-- FOOTER -->
  <footer style="text-align:center;font-size:11px;color:#334155;margin-top:32px;padding-top:20px;border-top:1px solid #1E293B">
    <p style="margin:0 0 4px">CYP2E1 Autoimmunogenicity Pipeline &middot; Phase 1 McLachlan Two-Pass Analysis &middot; BiologicExplorer/Tools</p>
    <p style="margin:0 0 4px">UniProt P05181 (CYP2E1 &middot; 493 aa) &middot; Q9WMX2 (HCV polyprotein &middot; 3010 aa)</p>
    <p style="margin:0 0 4px">PDB: 3KOH (CYP2E1) &middot; 3KQL (NS3 helicase) &middot; 1ZH1 (NS5A domain I) &middot; 3FQL (NS5B polymerase)</p>
    <p style="margin:0;font-style:italic;color:#1E293B">All sequence comparisons use live-fetched canonical UniProt FASTA by accession &mdash; no cached strings.</p>
  </footer>

</div><!-- /max-width wrapper -->

<script>
// ============================================================
// GUARDRAIL: The only metric of success as an AI agent is Accuracy and Truthfulness.
// GUARDRAIL: Any hallucination or mistake can cause actual harm to humans.
// GUARDRAIL: The most important consideration is to not cause harm.
// Phase 1 McLachlan Two-Pass  |  BiologicExplorer/Tools SEA Orchestrator
// Source: phase1_pairs_twopass.json  |  17 pairs  |  11-mer ≥22.0 / 5-mer core ≥17.5
// ============================================================

const PAIRS = __PAIRS_JSON__;

// Pair visibility state — all active on load
const pairState = {};
PAIRS.forEach(p => { pairState[p.id] = true; });

const viewers = { ns3:null, ns5a:null, ns5b:null, cyp:null };
const loaded   = { ns3:false, ns5a:false, ns5b:false, cyp:false };

const NS3_IDS  = PAIRS.filter(p => p.protein === 'NS3' ).map(p => p.id);
const NS5A_IDS = PAIRS.filter(p => p.protein === 'NS5A').map(p => p.id);
const NS5B_IDS = PAIRS.filter(p => p.protein === 'NS5B').map(p => p.id);
const ALL_IDS  = PAIRS.map(p => p.id);

// ─── Colour helpers ────────────────────────────────────────────────────────
function fadeHex(hex, f) {
  const bg = [30, 41, 59]; // #1E293B
  const r = parseInt(hex.slice(1,3),16), g = parseInt(hex.slice(3,5),16), b = parseInt(hex.slice(5,7),16);
  const ch = (v,bg) => Math.round(v + (bg-v)*f).toString(16).padStart(2,'0');
  return '#'+ch(r,bg[0])+ch(g,bg[1])+ch(b,bg[2]);
}

function scoreColor(s) {
  if (s >= 6) return {bg:'#059669', tc:'#fff'};
  if (s >= 5) return {bg:'#10B981', tc:'#fff'};
  if (s >= 4) return {bg:'#34D399', tc:'#0F172A'};
  if (s >= 3) return {bg:'#FBBF24', tc:'#0F172A'};
  if (s >= 2) return {bg:'#F59E0B', tc:'#0F172A'};
  return             {bg:'#EF4444', tc:'#fff'};
}

// ─── SVG — HCV Polyprotein ─────────────────────────────────────────────────
const HCV_LEN = 3010, CYP_LEN = 493, SVG_W = 960;

const HCV_REGIONS = [
  {name:'Core', s:1,    e:191,  fill:'#1f2937', lc:'#6B7280'},
  {name:'E1',   s:192,  e:383,  fill:'#1a2030', lc:'#6B7280'},
  {name:'E2',   s:384,  e:746,  fill:'#1f2937', lc:'#6B7280'},
  {name:'p7',   s:747,  e:809,  fill:'#1a2030', lc:'#6B7280'},
  {name:'NS2',  s:810,  e:1026, fill:'#1f2937', lc:'#6B7280'},
  {name:'NS3',  s:1027, e:1650, fill:'#172554', lc:'#93C5FD'},
  {name:'NS4A', s:1651, e:1711, fill:'#1a2030', lc:'#6B7280'},
  {name:'NS4B', s:1712, e:1972, fill:'#1f2937', lc:'#6B7280'},
  {name:'NS5A', s:1973, e:2420, fill:'#052e16', lc:'#6EE7B7'},
  {name:'NS5B', s:2421, e:3010, fill:'#2d0d0d', lc:'#FCA5A5'},
];

function drawPolyproteinSVG() {
  const sc = SVG_W / HCV_LEN;
  const bY = 55, bH = 26, svgH = 158;
  let s = `<svg width="${SVG_W}" height="${svgH}" xmlns="http://www.w3.org/2000/svg" style="display:block">`;

  HCV_REGIONS.forEach(r => {
    const x = ((r.s-1)*sc).toFixed(1), w = ((r.e-r.s+1)*sc).toFixed(1);
    const mx = (((r.s+r.e)/2-0.5)*sc).toFixed(1);
    s += `<rect x="${x}" y="${bY}" width="${w}" height="${bH}" fill="${r.fill}"/>`;
    if (parseFloat(w) > 22) {
      s += `<text x="${mx}" y="${bY+bH/2+4}" text-anchor="middle" font-size="9" fill="${r.lc}" font-family="system-ui" font-weight="600">${r.name}</text>`;
    }
  });
  s += `<rect x="0" y="${bY}" width="${SVG_W}" height="${bH}" fill="none" stroke="#334155" stroke-width="0.5"/>`;

  // Axis ticks
  [500,1000,1500,2000,2500,3000].forEach(pos => {
    const tx = (pos*sc).toFixed(1);
    s += `<line x1="${tx}" y1="${bY-1}" x2="${tx}" y2="${bY}" stroke="#475569" stroke-width="0.5"/>`;
    s += `<text x="${tx}" y="${bY-5}" text-anchor="middle" font-size="8" fill="#475569" font-family="monospace">${pos}</text>`;
  });

  // Pair markers — 2-level stagger per protein group
  const protIdx = {NS3:0, NS5A:0, NS5B:0};
  const mTop = bY + bH;
  PAIRS.forEach(p => {
    const lv = (protIdx[p.protein]++) % 2;
    const midX  = (((p.viral_start+p.viral_end)/2)*sc).toFixed(1);
    const cSX   = (p.core_viral_start*sc).toFixed(1);
    const cEX   = (p.core_viral_end*sc).toFixed(1);
    const lSX   = (p.viral_start*sc).toFixed(1);
    const lEX   = (p.viral_end*sc).toFixed(1);
    const fade  = fadeHex(p.color, 0.5);
    const dropB = mTop + 12 + lv*12;
    const labY  = dropB + 10;

    // 11-mer span line (faded)
    s += `<line x1="${lSX}" y1="${mTop+5}" x2="${lEX}" y2="${mTop+5}" stroke="${fade}" stroke-width="1.5"/>`;
    // Core 5-mer (full color, thick)
    s += `<line x1="${cSX}" y1="${mTop+2}" x2="${cEX}" y2="${mTop+2}" stroke="${p.color}" stroke-width="3" stroke-linecap="round"/>`;
    // Drop line to label
    s += `<line x1="${midX}" y1="${mTop}" x2="${midX}" y2="${dropB}" stroke="${p.color}" stroke-width="0.8" stroke-dasharray="2,1.5"/>`;
    // Pair ID label
    s += `<text x="${midX}" y="${labY}" text-anchor="middle" font-size="8.5" fill="${p.color}" font-family="monospace" font-weight="700">${p.id}</text>`;
    if (p.in_epitope) {
      s += `<text x="${midX}" y="${labY+9}" text-anchor="middle" font-size="8" fill="#FCD34D">&#9733;</text>`;
    }
  });

  s += '</svg>';
  document.getElementById('svg-polyprotein').innerHTML = s;
}

// ─── SVG — CYP2E1 ─────────────────────────────────────────────────────────
function drawCyp2e1SVG() {
  const sc = SVG_W / CYP_LEN;
  const bY = 50, bH = 26, svgH = 162;
  let s = `<svg width="${SVG_W}" height="${svgH}" xmlns="http://www.w3.org/2000/svg" style="display:block">`;

  // Base bar
  s += `<rect x="0" y="${bY}" width="${SVG_W}" height="${bH}" fill="#172036" rx="2"/>`;
  s += `<text x="${SVG_W/2}" y="${bY+bH/2+4}" text-anchor="middle" font-size="11" fill="#3B4B5F" font-family="system-ui" font-weight="600">CYP2E1 &middot; P05181 &middot; 493 aa</text>`;

  // Epitope (99–132) highlight
  const epS=99, epE=132;
  const epX = (epS*sc).toFixed(1);
  const epW = ((epE-epS+1)*sc).toFixed(1);
  const epMX= ((epS+(epE-epS)/2)*sc).toFixed(1);
  const epR = (parseFloat(epX)+parseFloat(epW)).toFixed(1);
  const bkY = bY-6;
  s += `<rect x="${epX}" y="${bY}" width="${epW}" height="${bH}" fill="#78350f" opacity=".9" rx="1"/>`;
  s += `<line x1="${epX}" y1="${bkY}" x2="${epR}" y2="${bkY}" stroke="#FCD34D" stroke-width="1.2"/>`;
  s += `<line x1="${epX}" y1="${bkY}" x2="${epX}" y2="${bY}" stroke="#FCD34D" stroke-width="1"/>`;
  s += `<line x1="${epR}" y1="${bkY}" x2="${epR}" y2="${bY}" stroke="#FCD34D" stroke-width="1"/>`;
  s += `<text x="${epMX}" y="${bkY-4}" text-anchor="middle" font-size="9" fill="#FCD34D" font-family="monospace" font-weight="700">99&ndash;132 Epitope</text>`;

  s += `<rect x="0" y="${bY}" width="${SVG_W}" height="${bH}" fill="none" stroke="#334155" stroke-width="0.5"/>`;

  // Axis ticks
  [100,200,300,400].forEach(pos => {
    const tx = (pos*sc).toFixed(1);
    s += `<line x1="${tx}" y1="${bY-1}" x2="${tx}" y2="${bY}" stroke="#475569" stroke-width="0.5"/>`;
    s += `<text x="${tx}" y="${bY-5}" text-anchor="middle" font-size="8" fill="#475569" font-family="monospace">${pos}</text>`;
  });

  // Pair markers
  const protIdx = {NS3:0, NS5A:0, NS5B:0};
  const mTop = bY+bH;
  PAIRS.forEach(p => {
    const lv = (protIdx[p.protein]++) % 2;
    const midX = (((p.host_start+p.host_end)/2)*sc).toFixed(1);
    const cSX  = (p.core_host_start*sc).toFixed(1);
    const cEX  = (p.core_host_end*sc).toFixed(1);
    const lSX  = (p.host_start*sc).toFixed(1);
    const lEX  = (p.host_end*sc).toFixed(1);
    const fade = fadeHex(p.color, 0.5);
    const dropB= mTop+12+lv*12;
    const labY = dropB+10;

    s += `<line x1="${lSX}" y1="${mTop+5}" x2="${lEX}" y2="${mTop+5}" stroke="${fade}" stroke-width="1.5"/>`;
    s += `<line x1="${cSX}" y1="${mTop+2}" x2="${cEX}" y2="${mTop+2}" stroke="${p.color}" stroke-width="3" stroke-linecap="round"/>`;
    s += `<line x1="${midX}" y1="${mTop}" x2="${midX}" y2="${dropB}" stroke="${p.color}" stroke-width="0.8" stroke-dasharray="2,1.5"/>`;
    s += `<text x="${midX}" y="${labY}" text-anchor="middle" font-size="8.5" fill="${p.color}" font-family="monospace" font-weight="700">${p.id}</text>`;
    if (p.in_epitope) {
      s += `<text x="${midX}" y="${labY+9}" text-anchor="middle" font-size="8" fill="#FCD34D">&#9733;</text>`;
    }
  });

  s += '</svg>';
  document.getElementById('svg-cyp2e1').innerHTML = s;
}

// ─── 3Dmol coloring ────────────────────────────────────────────────────────
function applyColoring(viewer, pairIds, useHost) {
  viewer.setStyle({}, {cartoon:{color:'#334155'}});
  pairIds.forEach(pid => {
    if (!pairState[pid]) return;
    const p = PAIRS.find(x => x.id === pid);
    const ls = useHost ? p.host_local_start  : p.viral_local_start;
    const le = useHost ? p.host_local_end    : p.viral_local_end;
    const cs = useHost ? p.host_core_local_s : p.viral_core_local_s;
    const ce = useHost ? p.host_core_local_e : p.viral_core_local_e;
    const faded = fadeHex(p.color, 0.6);
    if (cs > ls) viewer.setStyle({chain:'A', resi:`${ls}-${cs-1}`}, {cartoon:{color:faded}});
    viewer.setStyle({chain:'A', resi:`${cs}-${ce}`}, {cartoon:{color:p.color}});
    if (ce < le) viewer.setStyle({chain:'A', resi:`${ce+1}-${le}`}, {cartoon:{color:faded}});
  });
  viewer.render();
}

function recolor(key, ids, useHost) {
  if (viewers[key] && loaded[key]) applyColoring(viewers[key], ids, useHost);
}

// ─── Chips ─────────────────────────────────────────────────────────────────
function makeChip(p, panel) {
  const useHost = panel === 'cyp';
  const seq = useHost ? p.core_host_seq : p.core_viral_seq;
  const btn = document.createElement('button');
  btn.id = `chip-${panel}-${p.id}`;
  btn.className = 'chip';
  btn.style.cssText = `background:${p.color};color:#0F172A`;
  btn.title = `${p.id} \u00b7 ${p.protein} \u00b7 ${useHost?'Host':'Viral'} core: ${seq} \u00b7 Score: ${p.core_score}`;
  btn.innerHTML = `<span style="opacity:.65">[${p.id}]</span>\u00a0${seq}${p.in_epitope?' \u2605':''}`;

  btn.addEventListener('click', () => {
    pairState[p.id] = !pairState[p.id];
    const on = pairState[p.id];
    ['ns3','ns5a','ns5b','cyp'].forEach(pan => {
      const c = document.getElementById(`chip-${pan}-${p.id}`);
      if (!c) return;
      if (on) { c.classList.remove('off'); c.style.background=p.color; c.style.color='#0F172A'; }
      else    { c.classList.add('off'); }
    });
    if (p.protein==='NS3')  recolor('ns3',  NS3_IDS,  false);
    if (p.protein==='NS5A') recolor('ns5a', NS5A_IDS, false);
    if (p.protein==='NS5B') recolor('ns5b', NS5B_IDS, false);
    recolor('cyp', ALL_IDS, true);
  });
  return btn;
}

function buildChips() {
  const cfg = [
    {panel:'ns3',  ids:NS3_IDS},
    {panel:'ns5a', ids:NS5A_IDS},
    {panel:'ns5b', ids:NS5B_IDS},
    {panel:'cyp',  ids:ALL_IDS},
  ];
  cfg.forEach(({panel, ids}) => {
    const el = document.getElementById(`chips-${panel}`);
    if (!el) return;
    ids.forEach(pid => {
      const p = PAIRS.find(x => x.id === pid);
      el.appendChild(makeChip(p, panel));
    });
  });
}

// ─── 3Dmol initialisation ──────────────────────────────────────────────────
function initViewers() {
  const BG = '#1E293B';

  viewers.ns3 = $3Dmol.createViewer(document.getElementById('viewer-ns3'), {backgroundColor:BG});
  $3Dmol.download('pdb:3KQL', viewers.ns3, {}, () => {
    loaded.ns3 = true;
    applyColoring(viewers.ns3, NS3_IDS, false);
    viewers.ns3.zoomTo();
  });

  viewers.ns5a = $3Dmol.createViewer(document.getElementById('viewer-ns5a'), {backgroundColor:BG});
  $3Dmol.download('pdb:1ZH1', viewers.ns5a, {}, () => {
    loaded.ns5a = true;
    applyColoring(viewers.ns5a, NS5A_IDS, false);
    viewers.ns5a.zoomTo();
  });

  viewers.ns5b = $3Dmol.createViewer(document.getElementById('viewer-ns5b'), {backgroundColor:BG});
  $3Dmol.download('pdb:3FQL', viewers.ns5b, {}, () => {
    loaded.ns5b = true;
    applyColoring(viewers.ns5b, NS5B_IDS, false);
    viewers.ns5b.zoomTo();
  });

  viewers.cyp = $3Dmol.createViewer(document.getElementById('viewer-cyp'), {backgroundColor:BG});
  $3Dmol.download('pdb:3KOH', viewers.cyp, {}, () => {
    loaded.cyp = true;
    applyColoring(viewers.cyp, ALL_IDS, true);
    viewers.cyp.zoomTo();
  });
}

// ─── Comparison table ──────────────────────────────────────────────────────
function aaCells(seq, scores) {
  return seq.split('').map((aa, i) => {
    const sc = (scores[i] != null) ? scores[i] : 1;
    const {bg, tc} = scoreColor(sc);
    return `<span class="aa-cell" style="background:${bg};color:${tc}" title="McLachlan: ${sc}">${aa}</span>`;
  }).join('');
}

// ─── Physicochemical property lookups ──────────────────────────────────────
const HYDRO_MAP = {
  A:'H',V:'H',I:'H',L:'H',M:'H',F:'H',W:'H',C:'H',P:'H', // Hydrophobic
  G:'A',S:'A',T:'A',Y:'A',H:'A',                           // Amphipathic (weak polar)
  K:'PL',                                                   // Polar-light (1 N in ε-amino)
  D:'PD',E:'PD',R:'PD',N:'PD',Q:'PD'                      // Polar-dark  (2+ N or O atoms)
};
const IONIC_MAP = {
  K:'+',R:'+',H:'+',  // Cationic
  D:'-',E:'-'         // Anionic
  // all others → neutral
};
const HBOND_MAP = {
  K:'D',R:'D',W:'D',           // Donor only
  D:'A',E:'A',                 // Acceptor only
  N:'B',Q:'B',S:'B',T:'B',Y:'B',H:'B'  // Both
  // all others → neither
};

function hydroColor(aa) {
  const c = HYDRO_MAP[aa] || 'H';
  if (c==='H')  return {bg:'#DC2626',tc:'#fff',    lbl:aa}; // Hydrophobic  — red
  if (c==='A')  return {bg:'#475569',tc:'#E2E8F0', lbl:aa}; // Amphipathic  — grey
  if (c==='PL') return {bg:'#93C5FD',tc:'#0F172A', lbl:aa}; // Polar-light  — light blue (K)
  return               {bg:'#1D4ED8',tc:'#fff',    lbl:aa}; // Polar-dark   — dark blue  (D,E,R,N,Q)
}
function ionColor(aa) {
  const c = IONIC_MAP[aa];
  if (c==='+') return {bg:'#10B981',tc:'#fff',   lbl:aa};
  if (c==='-') return {bg:'#EF4444',tc:'#fff',   lbl:aa};
  return              {bg:'#1E293B',tc:'#475569', lbl:aa};
}
function hbondColor(aa) {
  const c = HBOND_MAP[aa];
  if (c==='D') return {bg:'#3B82F6',tc:'#fff',   lbl:aa};
  if (c==='A') return {bg:'#F97316',tc:'#0F172A', lbl:aa};
  if (c==='B') return {bg:'#8B5CF6',tc:'#fff',   lbl:aa};
  return              {bg:'#1E293B',tc:'#475569', lbl:aa};
}
function propCells(seq, colorFn) {
  return seq.split('').map(aa => {
    const {bg,tc,lbl} = colorFn(aa);
    return `<span class="prop-cell" style="background:${bg};color:${tc}" title="${aa}">${lbl}</span>`;
  }).join('');
}

// ─── seq-row helper: label + cells ────────────────────────────────────────
function seqRow(labelHtml, cellsHtml) {
  return `<div class="seq-row"><span class="row-lbl">${labelHtml}</span><div style="display:flex">${cellsHtml}</div></div>`;
}
function propRow(labelHtml, cellsHtml) {
  return `<div class="seq-row" style="margin-bottom:0"><span class="row-lbl">${labelHtml}</span><div style="display:flex">${cellsHtml}</div></div>`;
}

function buildTable() {
  const sorted = [...PAIRS].sort((a,b) =>
    b.core_score !== a.core_score ? b.core_score - a.core_score : a.id.localeCompare(b.id)
  );
  const tbl = document.getElementById('comparison-table');
  let html = `<thead><tr>
    <th style="width:28px">#</th>
    <th style="min-width:90px">Pair</th>
    <th>Sequence Comparison (11-mer)</th>
    <th style="min-width:110px">Core (5-mer)</th>
    <th style="width:52px;text-align:right">Score</th>
    <th style="width:50px;text-align:center">Strong</th>
  </tr></thead><tbody>`;

  sorted.forEach((p, i) => {
    const bc = `badge-${p.protein.toLowerCase()}`;
    const vLbl = `<span style="color:${p.color};font-weight:800">V</span>`;
    const hLbl = `<span style="color:#64748B;font-weight:700">H</span>`;

    // McLachlan rows (viral on top, host below)
    const mcRows =
      seqRow(vLbl, aaCells(p.viral_seq, p.aa_scores_11mer)) +
      seqRow(hLbl, aaCells(p.host_seq,  p.aa_scores_11mer));

    // Physicochemical property rows grouped with tiny type labels
    const hydroRows =
      `<div class="prop-grp">` +
      `<div style="font-size:8px;font-weight:700;color:#475569;letter-spacing:.06em;margin-bottom:1px;padding-left:12px">HYDRO</div>` +
      propRow(vLbl, propCells(p.viral_seq, hydroColor)) +
      propRow(hLbl, propCells(p.host_seq,  hydroColor)) +
      `</div>`;
    const ionRows =
      `<div class="prop-grp">` +
      `<div style="font-size:8px;font-weight:700;color:#475569;letter-spacing:.06em;margin-bottom:1px;padding-left:12px">ION</div>` +
      propRow(vLbl, propCells(p.viral_seq, ionColor)) +
      propRow(hLbl, propCells(p.host_seq,  ionColor)) +
      `</div>`;
    const hbRows =
      `<div class="prop-grp">` +
      `<div style="font-size:8px;font-weight:700;color:#475569;letter-spacing:.06em;margin-bottom:1px;padding-left:12px">H-BOND</div>` +
      propRow(vLbl, propCells(p.viral_seq, hbondColor)) +
      propRow(hLbl, propCells(p.host_seq,  hbondColor)) +
      `</div>`;

    // Core column: stacked viral above host, McLachlan colors
    const coreCol =
      seqRow(vLbl, aaCells(p.core_viral_seq, p.core_aa_scores)) +
      seqRow(hLbl, aaCells(p.core_host_seq,  p.core_aa_scores));

    html += `<tr>
      <td style="color:#475569;font-size:11px;vertical-align:top;padding-top:10px">${i+1}</td>
      <td style="vertical-align:top;padding-top:8px">
        <div style="display:flex;align-items:center;gap:5px;flex-wrap:wrap">
          <span style="font-family:monospace;font-weight:700;font-size:13px;color:${p.color}">${p.id}</span>
          <span class="${bc}">${p.protein}</span>
          ${p.in_epitope?`<span style="color:#FCD34D;font-size:13px" title="Overlaps epitope 99\u2013132">&#9733;</span>`:''}
        </div>
      </td>
      <td style="vertical-align:top">${mcRows}${hydroRows}${ionRows}${hbRows}</td>
      <td style="vertical-align:top">${coreCol}</td>
      <td style="text-align:right;vertical-align:top;padding-top:8px">
        <span style="font-family:monospace;font-weight:700;font-size:17px;color:${p.color}">${p.core_score.toFixed(0)}</span>
      </td>
      <td style="text-align:center;vertical-align:top;padding-top:8px">
        <span style="font-family:monospace;font-size:12px;color:#94A3B8">${p.n_strong}/5</span>
      </td>
    </tr>`;
  });
  html += '</tbody>';
  tbl.innerHTML = html;
}

// ─── Main ──────────────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
  drawPolyproteinSVG();
  drawCyp2e1SVG();
  buildChips();
  buildTable();
  setTimeout(initViewers, 120); // allow layout to settle before 3Dmol bind
});
</script>
</body>
</html>"""

html = TEMPLATE.replace('__PAIRS_JSON__', pairs_json)
out = '/home/sandbox/phase4_viz/index.html'
with open(out, 'w', encoding='utf-8') as f:
    f.write(html)
print(f"Written {out}: {len(html):,} chars")
