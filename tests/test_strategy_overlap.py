from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from analytics.analyses.variant_summary import (
    VARIANT_USECOLS,
    StrategyOverlap,
    build_variant_summary,
)
from analytics.reporting.document import render_html
from analytics.reporting.strategy_overlap import strategy_set_overlap_view
from analytics.reporting.variant_profile import build_variant_sections


def _config(html: str) -> dict:
    return json.loads(re.search(r'data-role="config">(.*?)</script>', html, re.S)[1])


def _run_js(assertions: str) -> None:
    script = Path('analytics/reporting/strategy_overlap.js').read_text()
    subprocess.run(
        ['node', '-e', 'const document={addEventListener(){}};\n' + script + assertions],
        check=True, capture_output=True, text=True,
    )


def test_set_overlap_from_rows_and_warm_cache(tmp_path: Path) -> None:
    """Distinct allele counts survive gene-context duplication and cached rendering."""
    rows = []
    for pos, gene, strategies in [
        (100, '1', 's1'), (100, '2', 's2'),
        (200, '1', 's1,s2,s3'), (300, '1', 's1'),
        (400, '1', 's2'), (500, '1', 's3'), (600, '1', 's1,s3'),
    ]:
        rows.append({**dict.fromkeys(VARIANT_USECOLS, ''),
            'variant_key': f'1:{pos}:A>G', 'gene_id': gene, 'event_type': 'snv',
            'ref': 'A', 'alt': 'G', 'strategies': strategies, 'lookup_status': 'ok',
            'vep_status': 'ok', 'vep_primary_consequence': 'intron_variant'})
    annotations = tmp_path / 'annotations.tsv.gz'
    pd.DataFrame(rows).to_csv(annotations, sep='\t', index=False)
    evidence = tmp_path / 'ortholog.tsv.gz'
    pd.DataFrame(columns=['strategy', 'target_context', 'taxonomic_scope', 'evidence_unit',
        'site_aligned_count', 'alt_support_count', 'gnomad_found_count',
        'gnomad_not_found_count', 'gnomad_lookup_failed_count']).to_csv(evidence, sep='\t', index=False)
    kwargs = dict(strategy_label=str, ortholog_evidence_summary_path=evidence)
    summary = build_variant_summary(annotations, tmp_path / 'derived', **kwargs)
    cached = build_variant_summary(annotations, tmp_path / 'derived', **kwargs)
    assert cached.cache_hit
    expected = _config(strategy_set_overlap_view(summary))
    assert _config(strategy_set_overlap_view(cached)) == expected
    assert expected['union'] == 6
    assert {s['key']: (s['size'], s['exclusive']) for s in expected['strategies']} == {
        's1': (4, 1), 's2': (3, 1), 's3': (3, 1)}
    _run_js('''
      const assert=require('node:assert/strict');
      const config=''' + json.dumps(expected) + ''';
      const union=StrategySetOverlap.compare(config,'s1','others');
      assert.equal(union.left,4);assert.equal(union.right,5);assert.equal(union.shared,3);
      assert.equal(union.onlyLeft,1);assert.equal(union.onlyRight,2);
      assert.equal(union.leftCoverage,.6);assert.equal(union.rightCoverage,.75);
      const pair=StrategySetOverlap.compare(config,'s1','s2');
      assert.equal(pair.shared,2);assert.equal(pair.union,5);
    ''')
    stats = summary.strategy_stats.assign(**{'Genes with result': 2})
    html = render_html([('overview', 'Overview', ['<p>Smoke test</p>']),
        ('candidates', 'Candidate Profile', build_variant_sections(cached, stats))])
    assert html.count('<h3>Selected set overlap</h3>') == 1
    assert 'Pairwise strategy overlap' in html
    assert 'strategy-overlap-jaccard' in html
    assert 'DOMContentLoaded' in html
    assert 'ResizeObserver' in html
    assert 'map-ont' not in json.dumps(expected)  # No demo strategies or counts in production.
    (tmp_path / 'candidate_profile_smoke.html').write_text(html)


def test_set_overlap_uses_labels_and_rejects_missing_counts() -> None:
    summary = SimpleNamespace(overlap=StrategyOverlap(
        ['nucmer', 'minimap2_asm20'], np.array([[4, 2], [2, 3]]),
        np.array([[4, 5], [5, 3]]), np.array([[1, .4], [.4, 1]])),
        unique_variant_count=5, unique_contribution=pd.DataFrame([
            {'Strategy': 'minimap2 asm20', 'Unique To Strategy': 1},
            {'Strategy': 'nucmer', 'Unique To Strategy': 2}]))
    config = _config(strategy_set_overlap_view(summary))
    assert config['strategies'] == [
        {'key': 'nucmer', 'label': 'nucmer', 'size': 4, 'exclusive': 2},
        {'key': 'minimap2_asm20', 'label': 'minimap2 asm20', 'size': 3, 'exclusive': 1}]
    summary.unique_contribution = summary.unique_contribution.iloc[:1]
    with pytest.raises(ValueError, match='one exclusive allele count per strategy'):
        strategy_set_overlap_view(summary)
    assert strategy_set_overlap_view(SimpleNamespace(overlap=None)) == ''


def test_browser_geometry_and_empty_denominators() -> None:
    _run_js('''
      const assert=require('node:assert/strict');
      for (const [left,right,shared] of [
        [10,20,0],[20,10,10],[10,20,10],[10,10,10],
        [0,20,0],[20,0,0],[0,0,0],[1000000,1,1],[1746987946,870869981,670127426],
      ]) {
        const g=StrategySetOverlap.layout({left,right,shared},736,340);
        if (!left && !right) {assert.equal(g,null);continue;}
        assert(Number.isFinite(g.scale) && Number.isFinite(g.d));
        const area=StrategySetOverlap.intersectionArea(g.r,g.s,g.d);
        assert(Math.abs(area-shared/Math.max(left,right))<1e-12);
        assert(Math.abs(Math.PI*g.r*g.r-left/Math.max(left,right))<1e-12);
        assert(Math.abs(Math.PI*g.s*g.s-right/Math.max(left,right))<1e-12);
      }
      const config={strategies:[
        {key:'a',label:'a',size:0,exclusive:0},{key:'b',label:'b',size:20,exclusive:20}],
        intersections:[[0,0],[0,20]],union:20};
      const forward=StrategySetOverlap.compare(config,'a','others');
      assert.equal(forward.leftCoverage,0);assert.equal(forward.rightCoverage,null);
      const reverse=StrategySetOverlap.compare(config,'b','a');
      assert.equal(reverse.leftCoverage,null);assert.equal(reverse.rightCoverage,0);
      assert.throws(()=>StrategySetOverlap.compare(config,'a','a'));
    ''')


def test_browser_controls_and_hidden_tab_resize() -> None:
    """Execute the shipped controller, including report-tab resize and Plotly selection."""
    _run_js(r'''
      const assert=require('node:assert/strict');
      class Element {
        constructor(tag) {this.tag=tag;this.children=[];this.events={};this.attrs={};
          this.style={};this.clientWidth=0;this.viewBox={baseVal:{width:736,height:340}};}
        appendChild(e) {this.children.push(e);return e;}
        append(...nodes) {this.children.push(...nodes);}
        replaceChildren(...nodes) {this.children=nodes;}
        setAttribute(k,v) {this.attrs[k]=v;if(k==='viewBox') {
          const a=v.split(' ').map(Number);this.viewBox.baseVal={width:a[2],height:a[3]};}}
        addEventListener(k,f) {this.events[k]=f;}
        getBoundingClientRect() {return {left:0,top:0,width:this.clientWidth,height:340};}
        setPointerCapture() {}
      }
      const roles=['first','second','legend','out','in','fit','diagram','coverage',
        'regions','detail','tooltip','config'];
      const nodes=Object.fromEntries(roles.map(r=>[r,new Element(r)]));
      nodes.config.textContent=JSON.stringify({strategies:[
        {key:'a',label:'A',size:4,exclusive:1},
        {key:'b',label:'B',size:3,exclusive:1},
        {key:'c',label:'C',size:3,exclusive:1}],
        intersections:[[4,2,2],[2,3,1],[2,1,3]],union:6});
      document.createElement=t=>new Element(t);
      document.createElementNS=(_,t)=>new Element(t);
      document.createTextNode=text=>({text});
      let resize,click;
      global.ResizeObserver=class {constructor(f){resize=f;}observe(){}};
      document.querySelector=()=>({on:(_,f)=>{click=f;}});
      const root={querySelector:s=>nodes[s.match(/"(.*?)"/)[1]],
        getBoundingClientRect:()=>({left:0,top:0,width:736})};
      StrategySetOverlap.init(root);
      assert.equal(nodes.diagram.children.length,0); // Hidden Candidate Profile.
      nodes.diagram.clientWidth=736;resize();
      const group=()=>nodes.diagram.children.find(e=>e.tag==='g');
      assert(group());
      assert.equal(nodes.coverage.children[0].textContent,'A contains 60.0% of Union of other strategies.');
      const before=group().attrs.transform;nodes.in.events.click();
      assert.notEqual(group().attrs.transform,before);nodes.fit.events.click();
      assert.equal(group().attrs.transform,before);
      nodes.second.events.change({target:{value:'b'}});
      assert.equal(nodes.coverage.children[0].textContent,'A contains 66.7% of B.');
      nodes.regions.children[0].events.click();assert(nodes.detail.textContent.includes('Variants: 2'));
      click({points:[{x:'C',y:'B'}]});
      assert.equal(nodes.first.value,'b');assert.equal(nodes.second.value,'c');
      assert.equal(nodes.coverage.children[0].textContent,'B contains 33.3% of C.');
      nodes.diagram.clientWidth=320;resize();assert.equal(nodes.diagram.viewBox.baseVal.width,320);
      nodes.first.events.change({target:{value:'c'}});
      assert.notEqual(nodes.first.value,nodes.second.value);
    ''')
