"use strict";

const StrategySetOverlap = (() => {
  const colors = ['#3878b6', '#dc8735', '#529867', '#a377bb', '#bc5b61', '#589fa8'];

  function compare(config, first, second) {
    const index = config.strategies.findIndex(s => s.key === first);
    const a = config.strategies[index];
    if (!a) throw new Error('Unknown first strategy: ' + first);
    const other = second === 'others';
    const secondIndex = config.strategies.findIndex(s => s.key === second);
    const b = config.strategies[secondIndex];
    if (!other && (!b || first === second)) throw new Error('Choose two different strategy sets');
    const left = a.size;
    const right = other ? config.union - a.exclusive : b.size;
    const shared = other ? left - a.exclusive : config.intersections[index][secondIndex];
    if (shared < 0 || shared > Math.min(left, right)) throw new Error('Inconsistent strategy intersection');
    return {
      left, right, shared, onlyLeft: left - shared, onlyRight: right - shared,
      union: left + right - shared,
      leftCoverage: right ? shared / right : null,
      rightCoverage: left ? shared / left : null,
      leftLabel: a.label, rightLabel: other ? 'Union of other strategies' : b.label,
      leftColor: colors[index % colors.length],
      rightColor: other ? '#798b9c' : colors[secondIndex % colors.length],
    };
  }

  function intersectionArea(r, s, d) {
    if (d >= r + s) return 0;
    if (d <= Math.abs(r - s)) return Math.PI * Math.min(r, s) ** 2;
    const clamp = x => Math.max(-1, Math.min(1, x));
    return r*r*Math.acos(clamp((d*d+r*r-s*s)/(2*d*r)))
      + s*s*Math.acos(clamp((d*d+s*s-r*r)/(2*d*s)))
      - .5*Math.sqrt(Math.max(0,(-d+r+s)*(d+r-s)*(d-r+s)*(d+r+s)));
  }

  function layout(counts, width, height) {
    const unit = Math.max(counts.left, counts.right);
    if (!unit) return null;
    const r = Math.sqrt(counts.left / unit / Math.PI);
    const s = Math.sqrt(counts.right / unit / Math.PI);
    const target = counts.shared / unit;
    let lo = Math.abs(r-s), hi = r+s;
    if (!target) lo = hi = hi + .05;
    else if (counts.shared === Math.min(counts.left, counts.right)) hi = lo;
    else for (let k=0; k<80; k++) {
      const mid = (lo+hi)/2;
      if (intersectionArea(r,s,mid) > target) lo = mid; else hi = mid;
    }
    const d = (lo+hi)/2;
    return {r,s,d,mid:(d+s-r)/2,cx:width/2,cy:height/2,
      scale:Math.min((width-30)/(d+s+r),(height-25)/(2*Math.max(r,s)))};
  }

  function init(root) {
    const config = JSON.parse(root.querySelector('[data-role="config"]').textContent);
    const $ = role => root.querySelector('[data-role="'+role+'"]');
    const svg = $('diagram'), ns = 'http://www.w3.org/2000/svg';
    const fmt = n => n.toLocaleString('en-US');
    const pct = n => n === null ? 'n/a' : (100*n).toFixed(1)+'%';
    const compact = n => n >= 1e9 ? (n/1e9).toFixed(3)+'B'
      : n >= 1e6 ? (n/1e6).toFixed(1)+'M' : fmt(n);
    let first = config.strategies.reduce((a,b) => b.size > a.size ? b : a).key;
    let second = config.strategies.length > 2 ? 'others' : config.strategies.find(s=>s.key!==first).key;
    let zoom=1,pan={x:0,y:0},geometry=null,group=null,drag=null,selected=null;

    function el(parent,tag,attrs={},value) {
      const node=document.createElementNS(ns,tag);
      for (const [k,v] of Object.entries(attrs)) node.setAttribute(k,v);
      if (value !== undefined) node.textContent=value;
      parent.appendChild(node);
      return node;
    }
    function addOption(select,key,label) {
      const option=document.createElement('option');
      option.value=key;option.textContent=label;select.appendChild(option);
    }
    function choices() {
      $('first').replaceChildren();$('second').replaceChildren();
      for (const s of config.strategies) {
        addOption($('first'),s.key,s.label);
        if (s.key !== first) addOption($('second'),s.key,s.label);
      }
      if (config.strategies.length > 2) addOption($('second'),'others','Union of other strategies');
      if (second===first) second=config.strategies.find(s=>s.key!==first).key;
      $('first').value=first;$('second').value=second;
    }
    function hideTip() { $('tooltip').hidden=true; }
    function tip(event,lines) {
      const box=root.getBoundingClientRect(),node=$('tooltip');node.replaceChildren();
      for (const line of lines) {const e=document.createElement('div');e.textContent=line;node.appendChild(e);}
      const width=Math.min(300,box.width-16);node.style.maxWidth=width+'px';
      node.style.left=Math.max(4,Math.min(event.clientX-box.left+12,box.width-width-8))+'px';
      node.style.top=(event.clientY-box.top+14)+'px';node.hidden=false;
    }
    function regionAt(x,y) {
      if (!geometry) return null;
      const {r,s,d}=geometry;
      const a=r>0 && x*x+y*y<=r*r,b=s>0 && (x-d)**2+y*y<=s*s;
      return a&&b?'shared':a?'left':b?'right':null;
    }
    function hit(event) {
      if (!geometry) return null;
      const box=svg.getBoundingClientRect(),g=geometry;
      const x=(event.clientX-box.left)*svg.viewBox.baseVal.width/box.width;
      const y=(event.clientY-box.top)*svg.viewBox.baseVal.height/box.height;
      return regionAt((x-g.cx-pan.x)/(g.scale*zoom)+g.mid,(y-g.cy-pan.y)/(g.scale*zoom));
    }
    function lines(region) {
      const m=compare(config,first,second);
      if (region==='shared') return ['Shared variants: '+fmt(m.shared),
        pct(m.rightCoverage)+' of '+m.leftLabel,pct(m.leftCoverage)+' of '+m.rightLabel];
      if (region==='left') return ['Present in '+m.leftLabel,'Absent from '+m.rightLabel,
        'Variants: '+fmt(m.onlyLeft),pct(m.left ? m.onlyLeft/m.left : null)+' of first set'];
      return ['Present in '+m.rightLabel,'Absent from '+m.leftLabel,
        'Variants: '+fmt(m.onlyRight),pct(m.right ? m.onlyRight/m.right : null)+' of second set'];
    }
    function transform() {
      if (geometry) {const g=geometry;group.setAttribute('transform',
        `translate(${g.cx+pan.x} ${g.cy+pan.y}) scale(${g.scale*zoom}) translate(${-g.mid} 0)`);}
    }
    function bestLabel(region) {
      const g=geometry;let best=null;
      for (let ix=0;ix<90;ix++) for (let iy=0;iy<65;iy++) {
        const x=-g.r+(g.d+g.s+g.r)*ix/89,y=Math.max(g.r,g.s)*(2*iy/64-1);
        if (regionAt(x,y)!==region) continue;
        const clear=Math.min(Math.abs(g.r-Math.hypot(x,y)),Math.abs(g.s-Math.hypot(x-g.d,y)));
        if (!best || clear>best.clear) best={x,y,clear};
      }
      return best;
    }
    function drawDiagram(m) {
      const width=svg.clientWidth;
      if (!width) return; // Candidate Profile is initially a hidden report tab.
      const height=width<400?270:340;
      svg.replaceChildren();svg.setAttribute('viewBox',`0 0 ${width} ${height}`);svg.style.height=height+'px';
      el(svg,'title',{},`${m.leftLabel}: ${fmt(m.left)}; ${m.rightLabel}: ${fmt(m.right)}; shared: ${fmt(m.shared)}.`);
      geometry=layout(m,width,height);
      if (!geometry) {el(svg,'text',{x:width/2,y:height/2,'text-anchor':'middle'},'No variants in either selected set');return;}
      group=el(svg,'g');const {r,s,d,scale}=geometry,defs=el(group,'defs');
      const clip=el(defs,'clipPath',{id:'strategy-set-overlap-clip'});el(clip,'circle',{cx:0,cy:0,r});
      if (r) el(group,'circle',{cx:0,cy:0,r,fill:m.leftColor,'fill-opacity':.25,stroke:m.leftColor,'stroke-width':1.4/scale});
      if (s) el(group,'circle',{cx:d,cy:0,r:s,fill:m.rightColor,'fill-opacity':.25,stroke:m.rightColor,'stroke-width':1.4/scale});
      if (m.shared) el(group,'circle',{cx:d,cy:0,r:s,fill:'#7968ad','fill-opacity':.3,'clip-path':'url(#strategy-set-overlap-clip)'});
      transform();
      for (const [region,count] of [['left',m.onlyLeft],['shared',m.shared],['right',m.onlyRight]]) {
        if (!count) continue;
        const p=bestLabel(region),value=compact(count);
        if (!p || p.clear*scale<value.length*3.7+4) continue;
        el(group,'text',{x:p.x,y:p.y,'text-anchor':'middle','dominant-baseline':'middle','font-size':12/scale},value);
      }
    }
    function drawRegions(m) {
      const container=$('regions');container.replaceChildren();
      for (const [region,label,count] of [['left','First set only',m.onlyLeft],['shared','Shared',m.shared],['right','Second set only',m.onlyRight]]) {
        const button=document.createElement('button');button.type='button';
        button.textContent=label+': '+compact(count);button.setAttribute('aria-pressed',String(selected===region));
        button.addEventListener('click',()=>{selected=region;drawRegions(compare(config,first,second));});
        button.addEventListener('pointermove',event=>tip(event,lines(region)));
        button.addEventListener('pointerleave',hideTip);container.appendChild(button);
      }
      $('detail').textContent=selected?lines(selected).join(' · '):'Select a region for its exact counts and fractions.';
    }
    function draw() {
      const m=compare(config,first,second);drawDiagram(m);drawRegions(m);
      $('legend').replaceChildren();
      for (const [label,size,color] of [[m.leftLabel,m.left,m.leftColor],[m.rightLabel,m.right,m.rightColor]]) {
        const span=document.createElement('span'),swatch=document.createElement('span');
        swatch.style.cssText=`display:inline-block;width:10px;height:10px;margin-right:6px;background:${color}`;
        span.append(swatch,document.createTextNode(label+' · '+compact(size)));$('legend').appendChild(span);
      }
      $('coverage').replaceChildren();
      for (const [a,b,fraction] of [[m.leftLabel,m.rightLabel,m.leftCoverage],[m.rightLabel,m.leftLabel,m.rightCoverage]]) {
        const line=document.createElement('div');line.textContent=a+' contains '+pct(fraction)+' of '+b+'.';$('coverage').appendChild(line);
      }
    }
    function select(a,b) {first=a;second=b;zoom=1;pan={x:0,y:0};selected=null;hideTip();choices();draw();}
    $('first').addEventListener('change',e=>select(e.target.value,second));
    $('second').addEventListener('change',e=>select(first,e.target.value));
    $('in').addEventListener('click',()=>{zoom=Math.min(4,zoom*1.25);hideTip();transform();});
    $('out').addEventListener('click',()=>{zoom=Math.max(.6,zoom/1.25);hideTip();transform();});
    $('fit').addEventListener('click',()=>{zoom=1;pan={x:0,y:0};hideTip();transform();});
    svg.addEventListener('pointermove',e=>{
      if (drag) {pan={x:drag.px+e.clientX-drag.x,y:drag.py+e.clientY-drag.y};transform();return;}
      const region=hit(e);if (region) tip(e,lines(region));else hideTip();
    });
    svg.addEventListener('pointerleave',hideTip);
    svg.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY,px:pan.x,py:pan.y};svg.setPointerCapture(e.pointerId);hideTip();});
    svg.addEventListener('pointerup',e=>{
      const moved=drag&&(Math.abs(e.clientX-drag.x)+Math.abs(e.clientY-drag.y)>4);drag=null;
      if (!moved) {selected=hit(e);drawRegions(compare(config,first,second));}
    });
    svg.addEventListener('pointercancel',()=>{drag=null;});
    svg.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.6,Math.min(4,zoom*Math.exp(-e.deltaY*.001)));hideTip();transform();},{passive:false});
    const heatmap=document.querySelector('#strategy-overlap-jaccard .plotly-graph-div');
    if (heatmap && typeof heatmap.on==='function') heatmap.on('plotly_click',event=>{
      const point=event.points[0],a=config.strategies.find(s=>s.label===point.y),b=config.strategies.find(s=>s.label===point.x);
      if (a && b && a.key!==b.key) select(a.key,b.key);
    });
    let lastWidth=0;
    new ResizeObserver(()=>{const width=svg.clientWidth;if (width && width!==lastWidth) {lastWidth=width;hideTip();draw();}}).observe(root);
    choices();draw();
  }

  return {compare,intersectionArea,layout,init};
})();

document.addEventListener('DOMContentLoaded',()=>{
  const root=document.getElementById('strategy-set-overlap');
  if (root) StrategySetOverlap.init(root);
});
