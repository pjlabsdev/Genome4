// Grafiken des Überblicks. Die Daten stehen als JSON in #g-daten; alles läuft ohne Netz.
(function () {
  'use strict';
  var D = JSON.parse(document.getElementById('g-daten').textContent);
  var L = 'de';
  var picked = null;
  var NS = 'http://www.w3.org/2000/svg';
  var W = 760;

  function t(de, en) { return L === 'de' ? de : en; }
  var nfCache = {};
  function fmt(x, d) {
    d = d || 0;
    if (x === null || x === undefined || !isFinite(x)) return '-';
    var k = L + d;
    if (!nfCache[k]) nfCache[k] = new Intl.NumberFormat(L === 'de' ? 'de-DE' : 'en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
    return nfCache[k].format(x);
  }
  function pct(x, d) { return (x === null || x === undefined || !isFinite(x)) ? '-' : fmt(100 * x, d === undefined ? 1 : d) + ' %'; }
  function mb(x) { return fmt(x / 1e6, x < 1e7 ? 1 : 0) + ' Mb'; }
  function short(n) { return n.replace(/^chr/, ''); }

  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS, tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function txt(parent, x, y, s, cls, anchor) {
    var e = el('text', { x: x, y: y, 'class': cls || 'ax', 'text-anchor': anchor || 'start' }, parent);
    e.textContent = s;
    return e;
  }
  var tip = document.createElement('div');
  tip.className = 'tip';
  document.body.appendChild(tip);
  // Zeilen: erste fett, weitere als [Wert, Bezeichnung] oder Text.
  function showTip(ev, lines) {
    tip.replaceChildren();
    lines.forEach(function (l, i) {
      var d = document.createElement('div');
      if (i === 0) { var b = document.createElement('b'); b.textContent = l; d.appendChild(b); }
      else if (Array.isArray(l)) {
        var v = document.createElement('b'); v.textContent = l[0]; d.appendChild(v);
        var m = document.createElement('span'); m.className = 'm'; m.textContent = ' ' + l[1]; d.appendChild(m);
      } else { d.className = 'm'; d.textContent = l; }
      tip.appendChild(d);
    });
    tip.style.display = 'block';
    var r = tip.getBoundingClientRect();
    var x = ev.clientX + 14, y = ev.clientY + 14;
    if (x + r.width > innerWidth - 8) x = ev.clientX - r.width - 14;
    if (y + r.height > innerHeight - 8) y = ev.clientY - r.height - 14;
    tip.style.left = x + 'px';
    tip.style.top = y + 'px';
  }
  function hideTip() { tip.style.display = 'none'; }
  function niceMax(x) {
    if (!(x > 0)) return 1;
    var p = Math.pow(10, Math.floor(Math.log10(x)));
    var ms = [1, 2, 2.5, 5, 10];
    for (var i = 0; i < ms.length; i++) if (ms[i] * p >= x * 0.9999) return ms[i] * p;
    return 10 * p;
  }
  function frame(host, H) {
    host.replaceChildren();
    var svg = el('svg', { viewBox: '0 0 ' + W + ' ' + H, 'class': 'chart', role: 'img' });
    host.appendChild(svg);
    return svg;
  }
  function yAxis(svg, H, m, ymax, fy, label) {
    var ph = H - m.t - m.b;
    for (var i = 0; i <= 4; i++) {
      var y = m.t + ph * (1 - i / 4);
      el('line', { x1: m.l, x2: W - m.r, y1: y, y2: y, 'class': i === 0 ? 'base' : 'grid' }, svg);
      txt(svg, m.l - 6, y + 4, fy(ymax * i / 4), 'ax', 'end');
    }
    if (label) {
      var yc = m.t + ph / 2;
      txt(svg, 14, yc, label, 'ax', 'middle').setAttribute('transform', 'rotate(-90 14 ' + yc + ')');
    }
  }
  function roundTop(x, y, w, h, r) {
    if (h <= 0) return '';
    r = Math.min(r, h, w / 2);
    return 'M' + x + ',' + (y + h) + 'V' + (y + r) + 'Q' + x + ',' + y + ' ' + (x + r) + ',' + y + 'H' + (x + w - r) + 'Q' + (x + w) + ',' + y + ' ' + (x + w) + ',' + (y + r) + 'V' + (y + h) + 'Z';
  }
  function covAxis(ymax) { var d = ymax < 1 ? 2 : (ymax < 10 ? 1 : 0); return function (v) { return fmt(v, d) + 'x'; }; }
  function pctAxis(v) { return fmt(100 * v, v > 0 && v < 0.01 ? 2 : (v > 0 && v < 0.05 ? 1 : 0)) + ' %'; }

  // Säulen. items: {label, value, tip, cls, onclick}
  function bars(host, items, o) {
    var H = o.h || 220, m = { l: 58, r: 12, t: 16, b: 34 };
    var svg = frame(host, H);
    var ymax = o.ymax || niceMax(Math.max.apply(null, items.map(function (d) { return d.value; }).concat([0])));
    yAxis(svg, H, m, ymax, o.fy || pctAxis, o.ylabel);
    var pw = W - m.l - m.r, ph = H - m.t - m.b, bw = pw / items.length;
    items.forEach(function (d, i) {
      var h = Math.max(0, ph * Math.min(d.value, ymax) / ymax), x = m.l + bw * i, w = Math.max(1, bw - 2);
      var g = el('g', { 'class': 'hit' }, svg);
      el('rect', { x: x, y: m.t, width: bw, height: ph, 'class': 'hitbox' }, g);
      if (h > 0) el('path', { d: roundTop(x + 1, m.t + ph - h, w, h, 4), 'class': d.cls || 'bar' }, g);
      if (!o.labelEvery || i % o.labelEvery === 0) txt(svg, x + bw / 2, H - m.b + 16, d.label, 'ax', 'middle');
      g.addEventListener('pointermove', function (e) { showTip(e, d.tip); });
      g.addEventListener('pointerleave', hideTip);
      if (d.onclick) { g.style.cursor = 'pointer'; g.addEventListener('click', d.onclick); }
    });
    (o.refs || []).forEach(function (r) {
      var y = m.t + ph * (1 - r.y / ymax);
      el('line', { x1: m.l, x2: W - m.r, y1: y, y2: y, 'class': 'ref' }, svg);
      txt(svg, W - m.r - 4, y - 4, r.label, 'ax refl', 'end');
    });
    return svg;
  }

  // Fläche mit Fadenkreuz. xs aufsteigend, ys Werte.
  function area(host, xs, ys, o) {
    var H = o.h || 220, m = { l: 58, r: 12, t: 16, b: o.xlabel ? 38 : 26 };
    var svg = frame(host, H);
    var xmin = o.xmin !== undefined ? o.xmin : xs[0], xmax = o.xmax !== undefined ? o.xmax : xs[xs.length - 1];
    var top = 0;
    for (var i = 0; i < ys.length; i++) if (isFinite(ys[i]) && ys[i] > top) top = ys[i];
    var ymax = o.ymax || niceMax(top);
    yAxis(svg, H, m, ymax, o.fy || pctAxis, o.ylabel);
    var pw = W - m.l - m.r, ph = H - m.t - m.b;
    function X(x) { return m.l + pw * (x - xmin) / ((xmax - xmin) || 1); }
    function Y(y) { return m.t + ph * (1 - Math.max(0, Math.min(y, ymax)) / ymax); }
    (o.spans || []).forEach(function (s) {
      el('rect', { x: X(s[0]), y: m.t, width: Math.max(0.5, X(s[1]) - X(s[0])), height: ph, 'class': 'gap', opacity: 0.6 }, svg);
    });
    var d = 'M' + X(xs[0]).toFixed(1) + ',' + Y(0).toFixed(1);
    for (i = 0; i < xs.length; i++) d += 'L' + X(xs[i]).toFixed(1) + ',' + Y(isFinite(ys[i]) ? ys[i] : 0).toFixed(1);
    d += 'L' + X(xs[xs.length - 1]).toFixed(1) + ',' + Y(0).toFixed(1) + 'Z';
    el('path', { d: d, 'class': o.cls || 'area' }, svg);
    (o.refs || []).forEach(function (r) {
      if (r.x !== undefined) {
        if (r.x < xmin || r.x > xmax) return;
        el('line', { x1: X(r.x), x2: X(r.x), y1: m.t, y2: m.t + ph, 'class': 'ref' }, svg);
        txt(svg, X(r.x) + 4, m.t + 10, r.label, 'ax refl');
      } else {
        if (r.y > ymax) return;
        el('line', { x1: m.l, x2: W - m.r, y1: Y(r.y), y2: Y(r.y), 'class': 'ref' }, svg);
        txt(svg, W - m.r - 4, Y(r.y) - 4, r.label, 'ax refl', 'end');
      }
    });
    var step = niceMax((xmax - xmin) / 8);
    for (var v = Math.ceil(xmin / step) * step; v <= xmax + 1e-9; v += step) txt(svg, X(v), H - m.b + 15, (o.fx || fmt)(v), 'ax', 'middle');
    if (o.xlabel) txt(svg, m.l + pw / 2, H - 4, o.xlabel, 'ax', 'middle');
    var hair = el('line', { y1: m.t, y2: m.t + ph, 'class': 'hair', visibility: 'hidden' }, svg);
    var dot = el('circle', { r: 4, 'class': 'dot', visibility: 'hidden' }, svg);
    var hit = el('rect', { x: m.l, y: m.t, width: pw, height: ph, 'class': 'hitbox' }, svg);
    hit.addEventListener('pointermove', function (e) {
      var r = svg.getBoundingClientRect();
      var xv = xmin + ((e.clientX - r.left) * W / r.width - m.l) / pw * (xmax - xmin);
      var lo = 0, hi = xs.length - 1;
      while (hi - lo > 1) { var mid = (lo + hi) >> 1; if (xs[mid] < xv) lo = mid; else hi = mid; }
      var k = Math.abs(xs[lo] - xv) <= Math.abs(xs[hi] - xv) ? lo : hi;
      hair.setAttribute('x1', X(xs[k])); hair.setAttribute('x2', X(xs[k])); hair.setAttribute('visibility', 'visible');
      dot.setAttribute('cx', X(xs[k])); dot.setAttribute('cy', Y(isFinite(ys[k]) ? ys[k] : 0)); dot.setAttribute('visibility', 'visible');
      showTip(e, o.tip(k));
    });
    hit.addEventListener('pointerleave', function () { hair.setAttribute('visibility', 'hidden'); dot.setAttribute('visibility', 'hidden'); hideTip(); });
    return svg;
  }

  // Farbe für Abdeckung relativ zum Mittel der Autosomen: blau weniger, grau gleich, rot mehr.
  function css(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function hex(c) { c = c.replace('#', ''); return [parseInt(c.slice(0, 2), 16), parseInt(c.slice(2, 4), 16), parseInt(c.slice(4, 6), 16)]; }
  var colLo, colMid, colHi;
  function readColors() { colLo = hex(css('--lo')); colMid = hex(css('--mid')); colHi = hex(css('--hi')); }
  function relColor(rel) {
    // ±10 % um das Mittel bleibt grau: so viel schwankt die Abdeckung ohnehin.
    var t2 = Math.max(-1, Math.min(1, rel - 1));
    var a = t2 < 0 ? colLo : colHi, f = Math.max(0, (Math.abs(t2) - 0.1) / 0.9);
    return 'rgb(' + [0, 1, 2].map(function (i) { return Math.round(colMid[i] + (a[i] - colMid[i]) * f); }).join(',') + ')';
  }

  var chroms = D.chromosomen;
  var auto = D.autosomen_mittel || 1;
  function byName(n) { for (var i = 0; i < chroms.length; i++) if (chroms[i].name === n) return chroms[i]; return null; }

  // Karyogramm: jedes Chromosom eine Leiste in Fächern zu 1 Mio. Basen.
  function karyo(host) {
    var rows = chroms.filter(function (c) { return !/^(chr)?(M|MT)$/.test(c.name); });
    var rowH = 15, gap = 7, m = { l: 40, r: 70, t: 10, b: 10 };
    var H = m.t + rows.length * (rowH + gap) - gap + m.b;
    var svg = frame(host, H);
    var maxLen = Math.max.apply(null, rows.map(function (c) { return c.laenge; }));
    var pw = W - m.l - m.r, agg = D.zusammenfassen;
    rows.forEach(function (c, k) {
      var y = m.t + k * (rowH + gap);
      txt(svg, m.l - 8, y + rowH - 3, short(c.name), 'ax', 'end');
      var g = el('g', { 'class': 'chrom' + (c.name === picked ? ' sel' : '') }, svg);
      var wTot = pw * c.laenge / maxLen;
      var clipId = 'k' + k;
      var cp = el('clipPath', { id: clipId }, svg);
      el('rect', { x: m.l, y: y, width: wTot, height: rowH, rx: rowH / 2 }, cp);
      var inner = el('g', { 'clip-path': 'url(#' + clipId + ')' }, g);
      el('rect', { x: m.l, y: y, width: wTot, height: rowH, 'class': 'gap' }, inner);
      var n = Math.ceil(c.abdeckung.length / agg);
      for (var j = 0; j < n; j++) {
        var from = j * agg, sum = 0, len = 0;
        for (var q = from; q < Math.min(from + agg, c.abdeckung.length); q++) {
          var bl = Math.min(D.fach, c.laenge - q * D.fach);
          sum += c.abdeckung[q] * bl; len += bl;
        }
        var mean = len > 0 ? sum / len : 0;
        if (mean <= 0.05 * auto) continue;
        el('rect', { x: m.l + pw * from * D.fach / maxLen, y: y, width: pw * len / maxLen + 0.4, height: rowH, fill: relColor(mean / auto), 'class': 'cell' }, inner);
      }
      el('rect', { x: m.l, y: y, width: wTot, height: rowH, rx: rowH / 2, 'class': 'outline' }, g);
      txt(svg, m.l + wTot + 6, y + rowH - 3, fmt(c.mittel, 1) + 'x', 'ax');
      var hb = el('rect', { x: m.l, y: y - gap / 2, width: wTot, height: rowH + gap, 'class': 'hitbox' }, g);
      hb.addEventListener('pointermove', function (e) {
        var r = svg.getBoundingClientRect();
        var pos = ((e.clientX - r.left) * W / r.width - m.l) / pw * maxLen;
        var q2 = Math.max(0, Math.min(c.abdeckung.length - 1, Math.floor(pos / D.fach)));
        var a = q2 * D.fach;
        showTip(e, [c.name + ' · ' + mb(a) + '-' + mb(Math.min(c.laenge, a + D.fach)),
          [fmt(c.abdeckung[q2], 1) + 'x', t('Abdeckung', 'coverage')],
          [fmt(c.abdeckung[q2] / auto, 2), t('× Mittel der Autosomen', '× autosome mean')],
          t('Klicken für Details', 'Click for details')]);
      });
      hb.addEventListener('pointerleave', hideTip);
      hb.addEventListener('click', function () { pick(c.name, true); });
    });
    return svg;
  }

  function relBars(host) {
    var items = chroms.map(function (c) {
      var rel = c.mittel / auto;
      return {
        label: short(c.name), value: Math.min(rel, 2),
        cls: /^(chr)?[XY]$/.test(c.name) ? 'bar2' : (/^(chr)?(M|MT)$/.test(c.name) ? 'bar3' : 'bar'),
        tip: [c.name, [fmt(rel, 3), t('× Mittel der Autosomen', '× autosome mean')], [fmt(c.mittel, 2) + 'x', t('mittlere Abdeckung', 'mean coverage')], t('Klicken für Details', 'Click for details')],
        onclick: function () { pick(c.name, true); }
      };
    });
    bars(host, items, {
      ymax: 1.6, h: 220, ylabel: t('× Mittel der Autosomen', '× autosome mean'), fy: function (v) { return fmt(v, 2); },
      refs: [{ y: 1, label: t('2 Kopien', '2 copies') }, { y: 0.5, label: t('1 Kopie', '1 copy') }]
    });
  }

  function chips(host) {
    host.replaceChildren();
    chroms.forEach(function (c) {
      var b = document.createElement('button');
      b.textContent = short(c.name);
      if (c.name === picked) b.className = 'on';
      b.addEventListener('click', function () { pick(c.name, false); });
      host.appendChild(b);
    });
  }

  function tile(host, v, k, n) {
    var d = document.createElement('div'); d.className = 'tile';
    var a = document.createElement('div'); a.className = 'v'; a.textContent = v; d.appendChild(a);
    var b = document.createElement('div'); b.className = 'k'; b.textContent = k; d.appendChild(b);
    var c = document.createElement('div'); c.className = 'n'; c.textContent = n || ''; d.appendChild(c);
    host.appendChild(d);
  }
  function h3(host, s) { var h = document.createElement('h3'); h.textContent = s; host.appendChild(h); }
  function note(host, s) { var p = document.createElement('p'); p.className = 'note'; p.textContent = s; host.appendChild(p); }
  function box(host) { var d = document.createElement('div'); d.className = 'viz'; host.appendChild(d); return d; }

  function detail(host) {
    host.replaceChildren();
    var c = byName(picked);
    if (!c) return;
    var tiles = document.createElement('div'); tiles.className = 'tiles'; host.appendChild(tiles);
    tile(tiles, fmt(c.mittel, 1) + 'x', t('mittlere Abdeckung', 'mean coverage'), fmt(c.mittel / auto, 2) + ' ' + t('× Autosomen', '× autosomes'));
    tile(tiles, mb(c.laenge), t('Länge', 'Length'), fmt(c.laenge) + ' ' + t('Basen', 'bases'));
    tile(tiles, fmt(c.reads), 'Reads', t('Hauptzeilen mit Ort', 'primary mapped'));
    tile(tiles, pct(c.reads > 0 ? c.duplikate / c.reads : null, 2), t('Duplikate', 'Duplicates'), '');
    tile(tiles, pct(c.ab[0], 1), t('Stellen mit ≥ 1x', 'Positions ≥ 1x'), t('N-Strecken eingerechnet', 'N gaps included'));
    tile(tiles, pct(c.ab[2], 1), t('Stellen mit ≥ 20x', 'Positions ≥ 20x'), '');
    tile(tiles, c.mapq === null ? '-' : fmt(c.mapq, 1), t('MAPQ im Mittel', 'Mean MAPQ'), t('ohne Duplikate', 'without duplicates'));

    var xs = c.abdeckung.map(function (_, i) { return Math.min(c.laenge, (i + 0.5) * D.fach); });
    var spans = [];
    var s0 = -1;
    for (var i = 0; i <= c.abdeckung.length; i++) {
      var zero = i < c.abdeckung.length && c.abdeckung[i] <= 0.01;
      if (zero && s0 < 0) s0 = i;
      if (!zero && s0 >= 0) { spans.push([s0 * D.fach, Math.min(c.laenge, i * D.fach)]); s0 = -1; }
    }
    function tipAt(k) {
      var a = k * D.fach;
      return [c.name + ' · ' + mb(a) + '-' + mb(Math.min(c.laenge, a + D.fach)),
        [fmt(c.abdeckung[k], 1) + 'x', t('Abdeckung', 'coverage')],
        [fmt(c.abdeckung[k] / auto, 2), t('× Mittel der Autosomen', '× autosome mean')],
        [fmt(c.fach_reads[k]), 'Reads'],
        [pct(c.mapq0[k], 1), 'MAPQ 0'],
        [pct(c.dup[k], 2), t('Duplikate', 'duplicates')]];
    }
    var fx = function (v) { return fmt(v / 1e6, 0); };
    var covYmax = niceMax(Math.max(auto * 2, c.mittel * 1.6));
    h3(host, t('Abdeckung entlang von ', 'Coverage along ') + c.name + t(' (je 100 kb)', ' (per 100 kb)'));
    area(box(host), xs, c.abdeckung, {
      xmin: 0, xmax: c.laenge, h: 240, ymax: covYmax,
      fy: covAxis(covYmax), fx: fx, xlabel: t('Position in Mb', 'Position in Mb'),
      ylabel: t('Abdeckung', 'Coverage'), spans: spans,
      refs: [{ y: auto, label: t('Mittel Autosomen ', 'autosome mean ') + fmt(auto, 1) + 'x' }],
      tip: tipAt
    });
    note(host, t('Graue Flächen: keine Reads, meist N-Strecken der Referenz (Zentromer, Lücken). Werte über dem oberen Rand sind abgeschnitten; der genaue Wert steht beim Darüberfahren.',
      'Grey areas: no reads, mostly N gaps of the reference (centromere, gaps). Values above the top are clipped; hover for the exact value.'));
    h3(host, t('Anteil der Reads mit MAPQ 0 (Wiederholungen)', 'Share of reads with MAPQ 0 (repeats)'));
    area(box(host), xs, c.mapq0, {
      xmin: 0, xmax: c.laenge, h: 170, cls: 'area2', fx: fx, ymax: 1, spans: spans,
      ylabel: t('Anteil', 'Share'), tip: tipAt
    });
    note(host, t('Hohe Werte zeigen Bereiche, in denen kurze Reads keinen eindeutigen Ort haben (Segmentduplikationen, Satelliten). Varianten dort sind unsicher.',
      'High values mark regions where short reads have no unique location (segmental duplications, satellites). Variants there are uncertain.'));
    h3(host, t('Anteil der Duplikate', 'Share of duplicates'));
    var dmax = 0; c.dup.forEach(function (v) { if (v > dmax) dmax = v; });
    area(box(host), xs, c.dup, {
      xmin: 0, xmax: c.laenge, h: 150, cls: 'area3', fx: fx, ymax: niceMax(Math.max(dmax, 0.01)), spans: spans,
      ylabel: t('Anteil', 'Share'), tip: tipAt
    });
    if (c.verteilung) {
      h3(host, t('Abdeckung je Stelle auf ', 'Per-position coverage on ') + c.name);
      var hx = c.verteilung.map(function (_, i) { return i; });
      var total = c.verteilung_stellen || c.verteilung.reduce(function (a, b) { return a + b; }, 0) || 1;
      var share = c.verteilung.map(function (v) { return v / total; });
      area(box(host), hx.slice(1), share.slice(1), {
        h: 200, xmin: 1, xlabel: t('Abdeckung', 'Coverage'), ylabel: t('Anteil der Stellen', 'Share of positions'),
        refs: [{ x: auto, label: t('Mittel Autosomen', 'autosome mean') }],
        tip: function (k) { return [fmt(k + 1) + 'x', [pct(share[k + 1], 2), t('der Stellen', 'of positions')]]; }
      });
      note(host, t('Stellen ohne Abdeckung (0x): ', 'Positions without coverage (0x): ') + pct(share[0], 2) + t('; sie fehlen in der Grafik. Rechts abgeschnitten bei ', '; not drawn. Cut off at ') + fmt(c.verteilung.length - 1) + 'x.');
    }
  }

  function pick(name, scroll) {
    picked = name;
    chips(document.getElementById('v-chips'));
    detail(document.getElementById('v-detail'));
    karyo(document.getElementById('v-karyo'));
    if (scroll) document.getElementById('chromosomen').scrollIntoView({ behavior: 'smooth' });
  }

  function hist(host, h, o) {
    var xs = [], ys = [], tot = 0;
    for (var i = 0; i < h.length; i++) tot += h[i];
    tot = o.total || tot || 1;
    for (i = o.from || 0; i < h.length; i++) { xs.push(i); ys.push(h[i] / tot); }
    o.tip = o.tip || function (k) { return [fmt(xs[k]) + (o.unit || ''), [pct(ys[k], 2), o.of], [fmt(h[xs[k]]), o.count]]; };
    return area(host, xs, ys, o);
  }

  function render() {
    readColors();
    if (picked === null) picked = chroms.length ? chroms[0].name : null;
    karyo(document.getElementById('v-karyo'));
    relBars(document.getElementById('v-rel'));
    chips(document.getElementById('v-chips'));
    detail(document.getElementById('v-detail'));

    var mq = D.mapq;
    var mqTot = mq.reduce(function (a, b) { return a + b.n; }, 0) || 1;
    bars(document.getElementById('v-mapq'), mq.map(function (b, i) {
      return { label: b.label, value: b.n / mqTot, cls: i === mq.length - 1 ? 'bar2' : 'bar', tip: ['MAPQ ' + b.label, [pct(b.n / mqTot, 2), t('der Reads', 'of reads')], [fmt(b.n), 'Reads']] };
    }), { ylabel: t('Anteil der Reads', 'Share of reads') });

    if (D.insert) {
      hist(document.getElementById('v-insert'), D.insert.verteilung, {
        h: 230, xlabel: t('Insertgröße in Basen', 'Insert size in bases'), ylabel: t('Anteil der Paare', 'Share of pairs'),
        total: D.insert.gesamt,
        refs: [{ x: D.insert.median, label: 'Median ' + fmt(D.insert.median) }], unit: ' bp', of: t('der Paare', 'of pairs'), count: t('Paare', 'pairs')
      });
    }
    if (D.read_laenge) {
      hist(document.getElementById('v-readlen'), D.read_laenge, {
        h: 200, from: 1, xlabel: t('Länge in Basen', 'Length in bases'), ylabel: t('Anteil der Reads', 'Share of reads'), unit: ' bp', of: t('der Reads', 'of reads'), count: 'Reads'
      });
    }
    var dh = D.tiefe.verteilung, dt = 0;
    dh.forEach(function (v) { dt += v; });
    dt = dt || 1;
    hist(document.getElementById('v-depth'), dh, {
      h: 230, total: D.tiefe.stellen, xlabel: t('Abdeckung', 'Coverage'), ylabel: t('Anteil der Stellen', 'Share of positions'),
      refs: [{ x: D.tiefe.median, label: 'Median ' + fmt(D.tiefe.median) + 'x' }], unit: 'x', of: t('der Stellen', 'of positions'), count: t('Stellen', 'positions')
    });
    var cxs = [], cys = [], rest = D.tiefe.ab_ende;
    for (var i = dh.length - 1; i >= 0; i--) { rest += dh[i]; cys[i] = rest; }
    var all = D.tiefe.stellen || 1;
    for (i = 0; i < dh.length; i++) { cxs.push(i); cys[i] = cys[i] / all; }
    area(document.getElementById('v-depth-cum'), cxs, cys, {
      h: 210, ymax: 1, cls: 'area', xlabel: t('mindestens so viel Abdeckung', 'at least this coverage'), ylabel: t('Anteil der Stellen', 'Share of positions'),
      refs: [{ x: 20, label: '20x' }, { x: 30, label: '30x' }],
      tip: function (k) { return ['≥ ' + fmt(k) + 'x', [pct(cys[k], 2), t('der Stellen auf den Autosomen', 'of autosome positions')]]; }
    });

    document.querySelectorAll('tr.pick').forEach(function (tr) {
      tr.onclick = function () { pick(tr.dataset.c, true); };
    });
  }

  window.setLang = function (x) {
    L = x === 'en' ? 'en' : 'de';
    document.body.classList.toggle('en', L === 'en');
    document.documentElement.lang = L;
    document.querySelectorAll('.lang button').forEach(function (b) { b.classList.toggle('on', b.dataset.l === L); });
    try { localStorage.setItem('genome.lang', L); } catch (e) { }
    render();
  };
  var start = 'de';
  try { start = localStorage.getItem('genome.lang') || 'de'; } catch (e) { }
  if (window.matchMedia) {
    var mm = window.matchMedia('(prefers-color-scheme: dark)');
    if (mm.addEventListener) mm.addEventListener('change', render);
  }
  window.setLang(start);
})();
