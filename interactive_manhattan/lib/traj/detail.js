/* detail.js -- interactive per-SNP trajectory panel for the interactive
   Manhattan plots (05f_interactive_manhattan.R).

   Clicking a significant SNP in the ggiraph Manhattan calls
   SNPDetail.show("CHROM:POS"). This draws the SNP's allele-frequency slopegraph
   (Gen 1 -> Final) as an inline SVG: 10 O-replicate lines (blue) and 10
   B-replicate lines (red). Hovering a line -- or a numbered chip in the key --
   highlights that replicate and reads out its values, so individual replicates
   are traceable even where the 20 lines overlap. Data comes from the shared
   window.SNP_TRAJ (trajectories.js, built by 05h_trajectory_data.R).

   Dependency-free: vanilla JS + inline SVG. */
(function () {
  "use strict";

  var SVGNS = "http://www.w3.org/2000/svg";
  var COL = { O: "#1F6FB2", B: "#C0392B" };
  var FINAL_GEN = { O: "20", B: "56" };   // final generation per regime
  var W = 470, H = 330, M = { t: 14, r: 24, b: 44, l: 44 };
  var plotW = W - M.l - M.r, plotH = H - M.t - M.b;
  var xG1 = M.l + 8, xFin = M.l + plotW - 8;
  function yOf(f) { return M.t + (1 - f) * plotH; }
  function fmt(v) { return v == null ? "n/a" : v.toFixed(3); }

  var panel, svg, gLines, readout, chipEls = {}, lineRefs = {}, pinned = null, curKey = null;

  function elem(tag, attrs, ns) {
    var e = ns ? document.createElementNS(ns, tag) : document.createElement(tag);
    if (attrs) for (var k in attrs) e.setAttribute(k, attrs[k]);
    return e;
  }
  function svgEl(tag, attrs) { return elem(tag, attrs, SVGNS); }

  function buildPanel() {
    panel = elem("div", { id: "snp-detail" });
    panel.hidden = true;
    panel.innerHTML =
      '<div class="sd-head">' +
        '<span class="sd-pos"></span>' +
        '<span class="sd-gene"></span>' +
        '<span class="sd-stats"></span>' +
      '</div>' +
      '<div class="sd-hint">Hover a line or a numbered chip to trace one replicate.</div>' +
      '<div class="sd-body">' +
        '<div class="sd-plot"></div>' +
        '<div class="sd-key">' +
          '<div class="sd-key-row" data-reg="O"><span class="sd-key-label">O</span><span class="sd-key-chips"></span></div>' +
          '<div class="sd-key-row" data-reg="B"><span class="sd-key-label">B</span><span class="sd-key-chips"></span></div>' +
          '<div class="sd-readout"></div>' +
        '</div>' +
      '</div>';
    document.body.appendChild(panel);
    readout = panel.querySelector(".sd-readout");
    buildKey();
    buildSvgFrame();
  }

  function buildKey() {
    ["O", "B"].forEach(function (reg) {
      var box = panel.querySelector('.sd-key-row[data-reg="' + reg + '"] .sd-key-chips');
      chipEls[reg] = [];
      for (var r = 0; r < 10; r++) {
        var chip = elem("div", { class: "sd-chip", title: reg + " replicate " + (r + 1) });
        chip.style.background = COL[reg];
        chip.textContent = String(r + 1);
        (function (reg, r, chip) {
          chip.addEventListener("mouseenter", function () { onEnter(reg, r); });
          chip.addEventListener("mouseleave", onLeave);
          chip.addEventListener("click", function () { onClick(reg, r); });
        })(reg, r, chip);
        box.appendChild(chip);
        chipEls[reg].push(chip);
      }
    });
  }

  function buildSvgFrame() {
    svg = svgEl("svg", { viewBox: "0 0 " + W + " " + H, role: "img" });
    // y gridlines + labels at 0, .25, .5, .75, 1
    [0, 0.25, 0.5, 0.75, 1].forEach(function (f) {
      var yy = yOf(f);
      svg.appendChild(svgEl("line", { class: f === 0 ? "sd-axis" : "sd-grid", x1: M.l, y1: yy, x2: xFin + 8, y2: yy }));
      var t = svgEl("text", { class: "sd-ylab", x: M.l - 6, y: yy + 4, "text-anchor": "end" });
      t.textContent = f.toFixed(2);
      svg.appendChild(t);
    });
    // y-axis line
    svg.appendChild(svgEl("line", { class: "sd-axis", x1: M.l, y1: M.t, x2: M.l, y2: M.t + plotH }));
    // y-axis title
    var yt = svgEl("text", { class: "sd-ylab", x: 12, y: M.t + plotH / 2, transform: "rotate(-90 12 " + (M.t + plotH / 2) + ")", "text-anchor": "middle" });
    yt.textContent = "ALT allele frequency";
    svg.appendChild(yt);
    // x ticks
    var xg = svgEl("text", { class: "sd-xtick", x: xG1, y: H - M.b + 18, "text-anchor": "middle" });
    xg.textContent = "Gen 1";
    svg.appendChild(xg);
    var xf = svgEl("text", { class: "sd-xtick", x: xFin, y: H - M.b + 18, "text-anchor": "middle" });
    xf.textContent = "Final";
    svg.appendChild(xf);
    var xf2 = svgEl("text", { class: "sd-ylab", x: xFin, y: H - M.b + 32, "text-anchor": "middle" });
    xf2.textContent = "O:20 · B:56";
    svg.appendChild(xf2);
    // groups: visible lines below, transparent hit lines on top
    gLines = svgEl("g", {});
    var gHit = svgEl("g", {});
    gHit.setAttribute("data-role", "hit");
    svg.appendChild(gLines);
    svg.appendChild(gHit);
    svg._gHit = gHit;
    panel.querySelector(".sd-plot").appendChild(svg);
  }

  function clearData() {
    while (gLines.firstChild) gLines.removeChild(gLines.firstChild);
    while (svg._gHit.firstChild) svg._gHit.removeChild(svg._gHit.firstChild);
    lineRefs = { O: [], B: [] };
  }

  function drawRegime(reg, pairs) {
    for (var r = 0; r < 10; r++) {
      var p = pairs[r];
      if (!p || p[0] == null || p[1] == null) { lineRefs[reg].push(null); continue; }
      var x1 = xG1, y1 = yOf(p[0]), x2 = xFin, y2 = yOf(p[1]);
      var line = svgEl("line", { class: "sd-line", x1: x1, y1: y1, x2: x2, y2: y2, stroke: COL[reg] });
      var d1 = svgEl("circle", { class: "sd-dot", cx: x1, cy: y1, r: 2.4, fill: COL[reg] });
      var d2 = svgEl("circle", { class: "sd-dot", cx: x2, cy: y2, r: 2.4, fill: COL[reg] });
      gLines.appendChild(line); gLines.appendChild(d1); gLines.appendChild(d2);
      var hit = svgEl("line", { class: "sd-hit", x1: x1, y1: y1, x2: x2, y2: y2 });
      (function (reg, r) {
        hit.addEventListener("mouseenter", function () { onEnter(reg, r); });
        hit.addEventListener("mouseleave", onLeave);
        hit.addEventListener("click", function () { onClick(reg, r); });
      })(reg, r);
      svg._gHit.appendChild(hit);
      lineRefs[reg].push({ line: line, d1: d1, d2: d2, v: p });
    }
  }

  function applyHL(reg, r) {
    var ref = lineRefs[reg] && lineRefs[reg][r];
    if (!ref) return;
    panel.classList.add("sd-hl");
    ["O", "B"].forEach(function (g) {
      (lineRefs[g] || []).forEach(function (o) {
        if (!o) return;
        o.line.classList.remove("sd-on"); o.d1.classList.remove("sd-on"); o.d2.classList.remove("sd-on");
      });
      chipEls[g].forEach(function (c) { c.classList.remove("sd-active"); });
    });
    ref.line.classList.add("sd-on"); ref.d1.classList.add("sd-on"); ref.d2.classList.add("sd-on");
    // raise the highlighted line + its dots above the other lines
    gLines.appendChild(ref.line); gLines.appendChild(ref.d1); gLines.appendChild(ref.d2);
    chipEls[reg][r].classList.add("sd-active");
    var d = ref.v[1] - ref.v[0], sign = d >= 0 ? "+" : "−";
    readout.innerHTML =
      '<span class="sd-swatch" style="background:' + COL[reg] + '"></span>' +
      '<b>' + reg + " replicate " + (r + 1) + "</b> &nbsp; " +
      "gen 1 " + fmt(ref.v[0]) + " → gen " + FINAL_GEN[reg] + " " + fmt(ref.v[1]) +
      ' &nbsp; <span class="sd-muted">(Δ ' + sign + fmt(Math.abs(d)) + ")</span>";
  }

  function clearHL() {
    panel.classList.remove("sd-hl");
    ["O", "B"].forEach(function (g) {
      (lineRefs[g] || []).forEach(function (o) {
        if (!o) return;
        o.line.classList.remove("sd-on"); o.d1.classList.remove("sd-on"); o.d2.classList.remove("sd-on");
      });
      chipEls[g].forEach(function (c) { c.classList.remove("sd-active"); });
    });
    readout.innerHTML = '<span class="sd-muted">Hover a line or a chip.</span>';
  }

  function onEnter(reg, r) { if (!pinned) applyHL(reg, r); }
  function onLeave() { if (pinned) applyHL(pinned.reg, pinned.r); else clearHL(); }
  function onClick(reg, r) {
    if (pinned && pinned.reg === reg && pinned.r === r) { pinned = null; clearHL(); }
    else { pinned = { reg: reg, r: r }; applyHL(reg, r); }
  }

  function show(key) {
    if (!panel) buildPanel();
    var data = window.SNP_TRAJ || {};
    var rec = data[key];
    curKey = key; pinned = null;
    var parts = key.split(":");
    panel.querySelector(".sd-pos").textContent = parts[0] + ":" + Number(parts[1]).toLocaleString("en-US");
    panel.hidden = false;
    if (!rec) {
      panel.querySelector(".sd-gene").textContent = "";
      panel.querySelector(".sd-stats").textContent = "";
      clearData(); clearHL();
      readout.innerHTML = '<span class="sd-muted">No trajectory data for this SNP.</span>';
      panel.scrollIntoView({ behavior: "smooth", block: "nearest" });
      return;
    }
    panel.querySelector(".sd-gene").textContent = rec.g || "";
    var stats = [];
    if (rec.c != null) stats.push("CMH −log₁₀p = " + rec.c);
    if (rec.a != null) stats.push("adapted = " + rec.a);
    if (rec.sig) stats.push("sig: " + rec.sig);
    panel.querySelector(".sd-stats").textContent = stats.join("  ·  ");
    clearData();
    drawRegime("O", rec.O || []);
    drawRegime("B", rec.B || []);
    clearHL();
    panel.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  // showEl: called from the ggiraph point's inline onclick as showEl(this);
  // the SNP key lives in the element's data-id attribute.
  function showEl(el) { if (el && el.getAttribute) show(el.getAttribute("data-id")); }

  window.SNPDetail = { show: show, showEl: showEl };
})();
