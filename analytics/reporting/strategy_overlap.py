"""Interactive two-set comparisons from existing allele-level summary counts."""

from __future__ import annotations

import json
from pathlib import Path

from analytics.analyses.variant_summary import VariantSummary
from .components import strategy_label


def strategy_set_overlap_view(summary: VariantSummary) -> str:
    overlap = summary.overlap
    if overlap is None or len(overlap.strategies) < 2:
        return ""
    labels = [strategy_label(key) for key in overlap.strategies]
    exclusive = summary.unique_contribution.set_index("Strategy")["Unique To Strategy"]
    if exclusive.index.has_duplicates or set(exclusive.index) != set(labels):
        raise ValueError("Strategy set overlap requires one exclusive allele count per strategy")
    sizes = [int(overlap.intersections[i, i]) for i in range(len(labels))]
    unique = [int(exclusive.loc[label]) for label in labels]
    union = int(summary.unique_variant_count)
    if any(not 0 <= count <= size <= union for count, size in zip(unique, sizes)):
        raise ValueError("Strategy set overlap has inconsistent allele counts")
    payload = {
        "strategies": [
            {"key": key, "label": label, "size": size, "exclusive": count}
            for key, label, size, count in zip(overlap.strategies, labels, sizes, unique)
        ],
        "intersections": overlap.intersections.tolist(),
        "union": union,
    }
    encoded = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    script = Path(__file__).with_suffix(".js").read_text()
    return """
    <div id="strategy-set-overlap">
      <h3>Selected set overlap</h3>
      <div class="set-overlap-controls">
        <label>First set<select data-role="first" aria-label="First set"></select></label>
        <label>Compare with<select data-role="second" aria-label="Compare with"></select></label>
      </div>
      <div data-role="legend" class="set-overlap-legend"></div>
      <div class="set-overlap-tools">
        <button type="button" data-role="out" aria-label="Zoom out">−</button>
        <button type="button" data-role="in" aria-label="Zoom in">+</button>
        <button type="button" data-role="fit">Fit</button>
      </div>
      <svg data-role="diagram" role="img" aria-label="Area-proportional overlap of the selected sets"></svg>
      <div data-role="coverage" class="set-overlap-coverage" aria-live="polite"></div>
      <div data-role="regions" class="set-overlap-regions"></div>
      <div data-role="detail" class="set-overlap-detail" aria-live="polite"></div>
      <p class="set-overlap-note">Each normalized allele is counted once. “Only” is relative to the two
      selected sets. Circle and overlap areas are proportional to variant counts; hover shows exact
      counts. Drag to pan, or use the wheel and buttons to zoom.</p>
      <div data-role="tooltip" class="set-overlap-tooltip" role="tooltip" hidden></div>
      <script type="application/json" data-role="config">""" + encoded + """</script>
    </div>
    <style>
      #strategy-set-overlap {position:relative;max-width:980px;}
      #strategy-set-overlap .set-overlap-controls {display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,260px),1fr));gap:12px;margin:12px 0;}
      #strategy-set-overlap label {color:#52606d;font-size:13px;}
      #strategy-set-overlap select {display:block;box-sizing:border-box;width:100%;margin-top:4px;padding:7px 8px;border:1px solid #cbd2d9;border-radius:4px;background:white;color:#1f2933;}
      #strategy-set-overlap button {padding:7px 12px;border:1px solid #cbd2d9;border-radius:4px;background:white;color:#1f2933;cursor:pointer;font:inherit;}
      #strategy-set-overlap button[aria-pressed="true"] {background:#e8eef5;border-color:#52606d;}
      #strategy-set-overlap .set-overlap-legend,#strategy-set-overlap .set-overlap-tools {display:flex;flex-wrap:wrap;gap:12px;margin:10px 0;}
      #strategy-set-overlap .set-overlap-legend span {overflow-wrap:anywhere;}
      #strategy-set-overlap svg {display:block;width:100%;touch-action:none;overflow:hidden;}
      #strategy-set-overlap svg text {fill:#1f2933;pointer-events:none;}
      #strategy-set-overlap .set-overlap-coverage {display:grid;gap:5px;margin:10px 0;font-variant-numeric:tabular-nums;}
      #strategy-set-overlap .set-overlap-regions {display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,210px),1fr));gap:10px;margin-top:14px;}
      #strategy-set-overlap .set-overlap-detail {min-height:2.5em;margin-top:12px;font-size:13px;font-variant-numeric:tabular-nums;}
      #strategy-set-overlap .set-overlap-note {color:#52606d;font-size:13px;}
      #strategy-set-overlap .set-overlap-tooltip {position:absolute;pointer-events:none;z-index:5;box-sizing:border-box;max-width:300px;background:#fff;color:#1f2933;padding:10px;border:1px solid #cbd2d9;border-radius:4px;font-size:13px;box-shadow:0 2px 8px #0002;}
    </style>
    <script>""" + script + """</script>"""
