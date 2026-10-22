/* Compact P/LP data only. All controls share one allele-level selection. */
const PathogenicView = (() => {
  function select(config, state) {
    return config.variants.filter(row =>
      (state.type === 'all' || row.type === state.type) &&
      (state.view !== 'union' || (
        (!state.orthologEnabled || (row.support !== null && row.support >= state.orthologMin)) &&
        (!state.strategyEnabled || row.strategies.length >= state.strategyMin)
      )));
  }
  function summarize(config, state) {
    const variants = select(config, state);
    const keys = new Set(variants.map(row => row.key));
    const labels = state.view === 'union' ? ['union'] : config.strategies.map(s => s.key);
    const groups = row => state.view === 'union' ? ['union'] : row.strategies;
    const stars = new Map(), consequences = new Map();
    const increment = (map, group, value) => {
      if (!map.has(group)) map.set(group, new Map());
      const counts = map.get(group);
      counts.set(value, (counts.get(value) || 0) + 1);
    };
    variants.forEach(row => groups(row).forEach(group => increment(stars, group,
      config.starOrder.includes(row.stars) ? row.stars : 'Unmapped')));
    // Context rows are unique allele–gene observations, with their own strategy membership.
    config.contexts.filter(row => keys.has(row.key)).forEach(row =>
      groups(row).forEach(group => increment(consequences, group, row.group)));
    const support = new Map(labels.map(label => [label, new Map()]));
    config.support.filter(row => keys.has(row.variant_key)).forEach(row => {
      const group = state.view === 'union' ? 'union' : row.strategy;
      const observations = support.get(group);
      const previous = observations.get(row.variant_key);
      if (!previous || row.alt_support_ortholog_count > previous.alt_support_ortholog_count ||
          (row.alt_support_ortholog_count === previous.alt_support_ortholog_count &&
           String(row.gene_id).localeCompare(String(previous.gene_id)) < 0)) {
        observations.set(row.variant_key, row);
      }
    });
    return {variants, keys, labels, stars, consequences, support};
  }
  function conditions(config, state, selectedStrategy, backgroundKind) {
    const strategy = state.view === 'union' ? 'union' : selectedStrategy;
    const variants = select(config, state).filter(row => strategy === 'union' || row.strategies.includes(strategy));
    const observed = new Map();
    variants.forEach(row => row.conditions.forEach(condition => {
      if (!observed.has(condition.key)) observed.set(condition.key, {name: condition.name, count: 0});
      observed.get(condition.key).count += 1;
    }));
    const background = new Map(config.background.filter(row =>
      row.cohort === backgroundKind && row.variant_type === state.type &&
      (backgroundKind === 'global' || row.strategy === strategy)
    ).map(row => [row.condition_key, {name: row.condition, count: row.variant_count, total: row.total_variant_count}]));
    return {observed, background, total: variants.length};
  }
  return {select, summarize, conditions};
})();

document.addEventListener('DOMContentLoaded', () => {
  const config = JSON.parse(document.getElementById('pathogenic-data').textContent);
  const controls = document.getElementById('pathogenic-controls');
  const conditionsControls = document.getElementById('pathogenic-condition-controls');
  const control = role => controls.querySelector('[data-role="' + role + '"]');
  const conditionControl = role => conditionsControls.querySelector('[data-role="' + role + '"]');
  const names = new Map(config.strategies.map(s => [s.key, s.label]));
  names.set('union', 'All strategies · Union');
  config.strategies.forEach(s => conditionControl('strategy').add(new Option(s.label, s.key)));
  const number = role => {
    const input = control(role);
    const parsed = Number(input.value);
    const value = Number.isFinite(parsed) && parsed >= 1 ? Math.floor(parsed) : Number(input.defaultValue);
    input.value = value;
    return value;
  };
  const state = () => ({
    view: control('view').value, type: control('type').value,
    orthologEnabled: control('ortholog-enabled').checked, orthologMin: number('ortholog-min'),
    strategyEnabled: control('strategy-enabled').checked, strategyMin: number('strategy-min'),
  });
  const label = (s, count) => {
    const filters = [];
    if (s.view === 'union' && s.orthologEnabled) filters.push('orthologs ≥ ' + s.orthologMin);
    if (s.view === 'union' && s.strategyEnabled) filters.push('strategies ≥ ' + s.strategyMin);
    return (s.view === 'union' ? 'All strategies · Union' : 'By strategy · unfiltered support') +
      ' · ' + ({all: 'SNV + INDEL', snv: 'SNV', indel: 'INDEL'}[s.type]) + ' · ' +
      (filters.join(' AND ') || 'no support filters') + ' · ' + count.toLocaleString() + ' unique P/LP alleles';
  };
  const safeText = text => String(text).replace(/[&<>"']/g, c => ({'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'}[c]));
  const plot = (id, traces, layout, hasData) => Plotly.react(id, traces, {
    template: 'plotly_white', height: 390, margin: {l: 70, r: 25, t: 45, b: 100},
    ...layout,
    annotations: hasData ? [] : [{text: 'No eligible observations for this selection.', x: .5, y: .5, xref: 'paper', yref: 'paper', showarrow: false}],
  }, {responsive: true});
  function bars(id, counts, labels, order, colors, yTitle) {
    const categories = order.filter(category => labels.some(group => counts.get(group)?.has(category)));
    const traces = categories.map(category => ({
      type: 'bar', name: category, marker: {color: colors[category]},
      x: labels.map(group => names.get(group)), y: labels.map(group => counts.get(group)?.get(category) || 0),
      hovertemplate: '%{x}<br>%{fullData.name}: %{y:,}<extra></extra>',
    }));
    plot(id, traces, {barmode: 'stack', xaxis: {title: {text: ''}}, yaxis: {title: {text: yTitle}, rangemode: 'tozero'}}, categories.length > 0);
  }
  function supportPlot(view) {
    const allCounts = [];
    const traces = view.labels.flatMap(group => {
      const rows = [...view.support.get(group).values()];
      if (!rows.length) return [];
      const counts = rows.map(row => row.alt_support_ortholog_count);
      allCounts.push(...counts);
      const points = counts.length <= 250 ? 'all' : false;
      const violin = new Set(counts).size >= 3;
      return [{
        type: violin ? 'violin' : 'box', name: names.get(group) + '<br>n=' + counts.length.toLocaleString(),
        y: counts.map(Math.log10), marker: {color: '#356d8f'}, line: {color: '#356d8f'},
        customdata: rows.map(r => [r.variant_key, r.gene_id, r.alt_support_ortholog_count, r.site_aligned_ortholog_count, r.alt_support_family_count]),
        hovertemplate: '%{customdata[0]}<br>Gene: %{customdata[1]}<br>Supporting orthologs: %{customdata[2]}<br>Site-aligned orthologs: %{customdata[3]}<br>Supporting families: %{customdata[4]}<extra>%{fullData.name}</extra>',
        ...(violin ? {points, box: {visible: true}, meanline: {visible: true}, spanmode: 'hard'} : {boxpoints: points}),
      }];
    });
    const maximum = allCounts.length ? Math.max(...allCounts) : 1;
    const ticks = [];
    for (let power = 1; power <= maximum; power *= 10) {
      [1, 2, 5].forEach(multiplier => {if (multiplier * power <= maximum) ticks.push(multiplier * power);});
    }
    plot('pathogenic-support', traces, {height: 440, showlegend: false, yaxis: {
      title: {text: 'Supporting orthologs (log scale)'}, tickmode: 'array', tickvals: ticks.map(Math.log10), ticktext: ticks.map(String),
    }}, traces.length > 0);
  }
  function renderConditions() {
    const s = state();
    const {observed, background, total} = PathogenicView.conditions(config, s, conditionControl('strategy').value, conditionControl('background').value);
    const fraction = (map, key, denominator) => (denominator ? (map.get(key)?.count || 0) / denominator : 0);
    const backgroundTotal = background.values().next().value?.total || 0;
    const score = key => {
      const a = fraction(observed, key, total), b = fraction(background, key, backgroundTotal);
      if (conditionControl('sort').value === 'clinvar') return b;
      if (conditionControl('sort').value === 'absolute-difference') return Math.abs(a - b);
      return a;
    };
    const name = key => observed.get(key)?.name || background.get(key)?.name || '';
    const query = conditionControl('search').value.trim().toLowerCase();
    const keys = [...new Set([...observed.keys(), ...background.keys()])].filter(key =>
      key && name(key).toLowerCase().includes(query)
    ).sort((a, b) => score(b) - score(a) || name(a).localeCompare(name(b))).slice(0, 15).reverse();
    const traces = [[observed, total, 'GAPH', '#356d8f'], [background, backgroundTotal, 'ClinVar', '#9ca3af']].map(([map, denominator, title, color]) => ({
      type: 'bar', orientation: 'h', name: title, marker: {color},
      y: keys.map(key => safeText(name(key))),
      x: keys.map(key => denominator ? fraction(map, key, denominator) : null),
      customdata: keys.map(key => [map.get(key)?.count || 0, denominator]),
      hovertemplate: '%{y}<br>Alleles: %{customdata[0]} / %{customdata[1]}<br>Fraction: %{x:.2%}<extra>%{fullData.name}</extra>',
    }));
    plot('pathogenic-conditions-plot', traces, {height: 600, barmode: 'group', margin: {l: 270, r: 25, t: 30, b: 70},
      xaxis: {title: {text: 'Fraction of P/LP alleles'}, tickformat: '.0%', rangemode: 'tozero'}, yaxis: {automargin: true},
    }, keys.length > 0 && total > 0);
  }
  function render() {
    const s = state(), view = PathogenicView.summarize(config, s);
    document.getElementById('pathogenic-union-filters').hidden = s.view !== 'union';
    control('ortholog-min').disabled = !s.orthologEnabled;
    control('strategy-min').disabled = !s.strategyEnabled;
    conditionControl('strategy').disabled = s.view === 'union';
    // Display the actual group in the disabled selector, rather than a stale strategy name.
    const selector = conditionControl('strategy');
    const unionOption = [...selector.options].find(o => o.value === 'union');
    if (s.view === 'union') {
      if (!unionOption) selector.add(new Option(names.get('union'), 'union'));
      selector.value = 'union';
    } else if (unionOption) {unionOption.remove(); selector.selectedIndex = 0;}
    const description = label(s, view.variants.length);
    document.getElementById('pathogenic-selection').textContent = description;
    const qcLabel = document.getElementById('pathogenic-qc-selection');
    if (qcLabel) qcLabel.textContent = description;
    const metrics = document.getElementById('pathogenic-metrics');
    metrics.replaceChildren();
    [
      ['Unique P/LP alleles', view.variants.length],
      ['With ≥2 review stars', view.variants.filter(row => Number(row.stars) >= 2).length],
      ['Named conditions', new Set(view.variants.flatMap(row => row.conditions.map(c => c.key))).size],
    ].forEach(([title, count]) => {
      const card = document.createElement('div'); card.className = 'metric-card';
      const caption = document.createElement('div'); caption.className = 'metric-label'; caption.textContent = title;
      const value = document.createElement('div'); value.className = 'metric-value'; value.textContent = count.toLocaleString();
      card.append(caption, value); metrics.append(card);
    });
    bars('pathogenic-stars', view.stars, view.labels, config.starOrder, config.starColors, 'Unique P/LP alleles');
    bars('pathogenic-consequences', view.consequences, view.labels, config.consequenceOrder, config.consequenceColors, 'Allele–target-gene observations');
    supportPlot(view);
    renderConditions();
    document.dispatchEvent(new CustomEvent('pathogenic-selection-change', {detail: {keys: [...view.keys]}}));
  }
  function openTab(tab, anchor) {
    document.querySelector('.tab-button[data-tab="' + tab + '"]')?.click();
    document.getElementById(anchor)?.scrollIntoView({block: 'start', behavior: 'smooth'});
  }
  controls.querySelectorAll('select,input').forEach(input => input.addEventListener('change', render));
  conditionsControls.querySelectorAll('select').forEach(input => input.addEventListener('change', renderConditions));
  conditionControl('search').addEventListener('input', renderConditions);
  document.querySelectorAll('[data-pathogenic-reset]').forEach(button => button.addEventListener('click', () => {
    control('type').value = 'all'; control('ortholog-enabled').checked = false; control('strategy-enabled').checked = false;
    control('ortholog-min').value = 6; control('strategy-min').value = 2; render();
  }));
  document.getElementById('pathogenic-open-qc').addEventListener('click', () => openTab('qc', 'pathogenic-qc-cases'));
  document.getElementById('pathogenic-edit-filters')?.addEventListener('click', () => openTab('pathogenic-clinvar-hits', 'pathogenic-controls'));
  render();
});
