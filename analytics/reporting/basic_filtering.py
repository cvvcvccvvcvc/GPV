"""Descriptive, context-specific candidate filtering controls."""

from __future__ import annotations

import json

from analytics.analyses.basic_filtering import (
    BasicFilteringAnalysis,
    FILTER_OPTIONS,
    UNION_STRATEGY,
    CANDIDATE_CONSEQUENCE_OPTIONS as CONSEQUENCE_OPTIONS,
)
from analytics.analyses.conservation_validation import TARGET_CONTEXT_OPTIONS, VARIANT_TYPE_OPTIONS
from .components import dataframe_records, strategy_label


def build_basic_filtering_sections(analysis: BasicFilteringAnalysis) -> list[str]:
    return ["<h2>Basic Filtering</h2>", basic_filtering_view(analysis)]


def basic_filtering_view(analysis: BasicFilteringAnalysis) -> str:
    strategies = (
        analysis.candidate_curves["strategy"].drop_duplicates().astype(str).tolist()
        if not analysis.candidate_curves.empty
        else []
    )
    payload = {
        "filters": [{"key": key, "label": label} for key, label, _column in FILTER_OPTIONS],
        "strategies": [
            {"key": key, "label": strategy_label(key)}
            for key in strategies
            if key != UNION_STRATEGY
        ],
        "variantTypes": [{"key": key, "label": label} for key, label in VARIANT_TYPE_OPTIONS],
        "contexts": [{"key": key, "label": label} for key, label in TARGET_CONTEXT_OPTIONS],
        "consequences": [{"key": key, "label": label} for key, label in CONSEQUENCE_OPTIONS],
        "candidate": dataframe_records(analysis.candidate_curves),
    }
    payload_json = json.dumps(payload, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    return (
        """
    <div class="analysis-controls" id="basic-filtering-global-controls">
      <label>Filter<select data-role="filter"></select></label>
      <label>Variant type<select data-role="variant-type"></select></label>
      <label>Target context<select data-role="target-context"></select></label>
      <label>Consequence (VEP)<select data-role="consequence"></select></label>
      <label>Threshold scale<select data-role="threshold-scale"><option value="log">Logarithmic (N + 1)</option><option value="linear">Linear</option></select></label>
    </div>
    <div class="analysis-controls analysis-controls-single" id="basic-filtering-candidate-controls">
      <label>Strategies<select data-role="candidate-strategy"></select></label>
    </div>
    <div id="basic-filtering-candidate-status" class="analysis-note" hidden></div>
    <div id="basic-filtering-retention" class="analysis-plot"></div>
    <div id="basic-filtering-gnomad" class="analysis-plot"></div>
    <p class="analysis-note">Consequences use candidate VEP terms and may overlap (for example, synonymous and splice region).
    Context and consequence must match the same target-gene row; each allele is counted once per selection.
    gnomAD presence describes population overlap, not clinical benignity. Failed lookups are excluded only from its denominator.
    The logarithmic view uses N + 1 so threshold zero remains visible; ticks and hover show N.</p>
    <script>(() => {
      const config = """
        + payload_json
        + """;
      const roles = ['filter', 'variant-type', 'candidate-strategy', 'target-context', 'consequence', 'threshold-scale'];
      const selects = Object.fromEntries(roles.map(role => [role,
        document.querySelector('#tab-basic-filtering [data-role="' + role + '"]')]));
      const addOptions = (select, options) => options.forEach(item => select.add(new Option(item.label, item.key)));
      const union = {key: 'union', label: 'Union (any strategy)'};
      addOptions(selects.filter, config.filters);
      addOptions(selects['variant-type'], config.variantTypes);
      addOptions(selects['candidate-strategy'], [{key: 'compare', label: 'Compare all strategies'}, ...config.strategies, union]);
      addOptions(selects['target-context'], config.contexts);
      addOptions(selects.consequence, config.consequences);
      const colors = ['#356d8f', '#e1812c', '#499c68', '#9b5c98', '#c24d4d', '#79706e'];
      const labels = Object.fromEntries([...config.strategies, union].map(item => [item.key, item.label]));
      const color = strategy => strategy === 'union' ? '#242424' : colors[config.strategies.findIndex(item => item.key === strategy) % colors.length];
      const count = value => Number(value || 0).toLocaleString('en-US');
      const atMost = () => selects.filter.value === 'aligned_max';
      const symbol = () => atMost() ? '≤' : '≥';
      const supported = () => selects['variant-type'].value === 'snv' || ['ortholog', 'strategy'].includes(selects.filter.value);
      const commonRows = row => row.filter_key === selects.filter.value && row.variant_type === selects['variant-type'].value
        && row.target_context === selects['target-context'].value && row.consequence === selects.consequence.value;
      const logarithmic = () => selects['threshold-scale'].value === 'log';
      const xaxis = (groups) => {
        if (!logarithmic()) return {title: {text: 'Threshold (' + symbol() + ' N)'}, type: 'linear', rangemode: 'tozero'};
        const max = Math.max(1, ...groups.map(group => group.rows.length ? group.rows[group.rows.length - 1].threshold : 0));
        const ticks = [0];
        for (let power = 1; power <= max; power *= 10) {
          for (const multiplier of [1, 2, 5]) if (power * multiplier <= max) ticks.push(power * multiplier);
        }
        return {title: {text: 'Threshold (' + symbol() + ' N; log scale of N + 1)'}, type: 'log',
          tickvals: ticks.map(value => value + 1), ticktext: ticks.map(String)};
      };
      function note(id, text) {
        const element = document.getElementById(id);
        element.textContent = text;
        element.hidden = !text;
      }
      function renderCandidate() {
        const selection = selects['candidate-strategy'].value;
        const strategies = selection === 'compare' ? config.strategies.map(item => item.key) : [selection];
        const groups = strategies.map(strategy => ({strategy, rows: config.candidate.filter(row =>
          supported() && commonRows(row) && row.strategy === strategy).sort((a, b) => a.threshold - b.threshold)}));
        const traces = (gnomad) => groups.filter(group => group.rows.length).map(({strategy, rows}) => ({
          type: 'scatter', mode: 'lines', name: labels[strategy], line: {color: color(strategy), width: 2, shape: 'hv'},
          connectgaps: false, x: rows.map(row => logarithmic() ? row.threshold + 1 : row.threshold),
          y: rows.map(row => gnomad ? row.gnomad_found_fraction : row.retained_fraction),
          customdata: rows.map(row => gnomad
            ? [row.threshold, count(row.gnomad_found_count), count(row.gnomad_eligible_count), count(row.gnomad_lookup_failed_count)]
            : [row.threshold, count(row.retained_variant_count), count(row.total_variant_count)]),
          hovertemplate: 'Threshold ' + symbol() + ' %{customdata[0]}<br>' + (gnomad
            ? 'gnomAD: %{customdata[1]} / %{customdata[2]}<br>Found: %{y:.2%}<br>Failed lookups excluded: %{customdata[3]}'
            : 'Retained: %{customdata[1]} / %{customdata[2]}<br>Fraction: %{y:.2%}') + '<extra>%{fullData.name}</extra>',
        }));
        for (const [id, gnomad, title, ylabel] of [
          ['retention', false, 'Candidate alleles retained', 'Fraction retained'],
          ['gnomad', true, 'Exact alleles found in gnomAD', 'gnomAD found fraction'],
        ]) {
          Plotly.react('basic-filtering-' + id, traces(gnomad), {
            title: {text: title}, template: 'plotly_white', height: 410, margin: {l: 72, r: 25, t: 55, b: 110},
            xaxis: xaxis(groups), yaxis: {title: {text: ylabel}, tickformat: '.0%', range: [0, 1]},
            legend: {orientation: 'h', x: 0, y: -0.3, yanchor: 'top'},
          }, {responsive: true});
        }
        note('basic-filtering-candidate-status', !supported() ? 'This filter is available for SNVs only.'
          : groups.every(group => !group.rows.length) ? 'No candidate alleles for this selection.' : '');
      }
      roles.forEach(role => selects[role].addEventListener('change', renderCandidate));
      renderCandidate();
    })();</script>
    """
    )
