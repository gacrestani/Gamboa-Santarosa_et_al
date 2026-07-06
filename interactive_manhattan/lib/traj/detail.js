/* detail.js -- interactive per-SNP trajectory panel for the interactive
   Manhattan plots (05f_interactive_manhattan.R).

   Clicking a significant SNP in the ggiraph Manhattan calls
   SNPDetail.show("CHROM:POS"). This draws the SNP's allele-frequency slopegraph
   (Gen 1 -> Final) as an inline SVG, grouped by the four Shahrestani treatments
   (5 replicates each): selected OBO/OB (blues, gen 1 -> 20) and control
   nBO/nB (reds, gen 1 -> 56). Hovering a line -- or a numbered chip in the key
   -- highlights that replicate and reads out its values, so individual
   replicates are traceable even where the 20 lines overlap. Data comes from the
   shared window.SNP_TRAJ (trajectories.js, built by 05h_trajectory_data.R).

   Dependency-free: vanilla JS + inline SVG. */
(function () {
  "use strict";

  var SVGNS = "http://www.w3.org/2000/svg";
  var COL = { OBO: "#1F6FB2", OB: "#6DB3E8", nBO: "#C0392B", nB: "#E8896B" };
  var TREATS = ["OBO", "OB", "nBO", "nB"];   // draw/iteration order
  var NREP = 5;
  var GROUPS = [
    { label: "selected (O)", treats: ["OBO", "OB"] },
    { label: "control (B)",  treats: ["nBO", "nB"] },
  ];
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
        '<div class="sd-key"></div>' +
      '</div>';
    document.body.appendChild(panel);
    buildKey();
    buildSvgFrame();
    // Clicking anywhere that is not a line hit-area or a chip clears the pin
    // (those handlers stopPropagation, so this only fires for "outside" clicks).
    document.addEventListener("click", function () {
      if (pinned) { pinned = null; clearHL(); }
    });
  }

  function buildKey() {
    var key = panel.querySelector(".sd-key");
    GROUPS.forEach(function (grp) {
      var gh = elem("div", { class: "sd-key-group" });
      gh.textContent = grp.label;
      key.appendChild(gh);
      grp.treats.forEach(function (tr) {
        var row = elem("div", { class: "sd-key-row" });
        row.setAttribute("data-reg", tr);
        var lab = elem("span", { class: "sd-key-label" });
        lab.textContent = tr;
        var chips = elem("span", { class: "sd-key-chips" });
        row.appendChild(lab); row.appendChild(chips);
        key.appendChild(row);
        chipEls[tr] = [];
        for (var r = 0; r < NREP; r++) {
          var chip = elem("div", { class: "sd-chip", title: tr + " replicate " + (r + 1) });
          chip.style.background = COL[tr];
          chip.textContent = String(r + 1);
          (function (tr, r, chip) {
            chip.addEventListener("mouseenter", function () { onEnter(tr, r); });
            chip.addEventListener("mouseleave", onLeave);
            chip.addEventListener("click", function (e) { e.stopPropagation(); onClick(tr, r); });
          })(tr, r, chip);
          chips.appendChild(chip);
          chipEls[tr].push(chip);
        }
      });
    });
    readout = elem("div", { class: "sd-readout" });
    key.appendChild(readout);
  }

  function buildSvgFrame() {
    svg = svgEl("svg", { viewBox: "0 0 " + W + " " + H, role: "img" });
    [0, 0.25, 0.5, 0.75, 1].forEach(function (f) {
      var yy = yOf(f);
      svg.appendChild(svgEl("line", { class: f === 0 ? "sd-axis" : "sd-grid", x1: M.l, y1: yy, x2: xFin + 8, y2: yy }));
      var t = svgEl("text", { class: "sd-ylab", x: M.l - 6, y: yy + 4, "text-anchor": "end" });
      t.textContent = f.toFixed(2);
      svg.appendChild(t);
    });
    svg.appendChild(svgEl("line", { class: "sd-axis", x1: M.l, y1: M.t, x2: M.l, y2: M.t + plotH }));
    var yt = svgEl("text", { class: "sd-ylab", x: 12, y: M.t + plotH / 2, transform: "rotate(-90 12 " + (M.t + plotH / 2) + ")", "text-anchor": "middle" });
    yt.textContent = "ALT allele frequency";
    svg.appendChild(yt);
    var xg = svgEl("text", { class: "sd-xtick", x: xG1, y: H - M.b + 18, "text-anchor": "middle" });
    xg.textContent = "Gen 1";
    svg.appendChild(xg);
    var xf = svgEl("text", { class: "sd-xtick", x: xFin, y: H - M.b + 18, "text-anchor": "middle" });
    xf.textContent = "Final";
    svg.appendChild(xf);
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
    lineRefs = {};
    TREATS.forEach(function (tr) { lineRefs[tr] = []; });
  }

  function drawTreat(tr, pairs) {
    for (var r = 0; r < NREP; r++) {
      var p = pairs[r];
      if (!p || p[0] == null || p[1] == null) { lineRefs[tr].push(null); continue; }
      var x1 = xG1, y1 = yOf(p[0]), x2 = xFin, y2 = yOf(p[1]);
      var line = svgEl("line", { class: "sd-line", x1: x1, y1: y1, x2: x2, y2: y2, stroke: COL[tr] });
      var d1 = svgEl("circle", { class: "sd-dot", cx: x1, cy: y1, r: 2.4, fill: COL[tr] });
      var d2 = svgEl("circle", { class: "sd-dot", cx: x2, cy: y2, r: 2.4, fill: COL[tr] });
      gLines.appendChild(line); gLines.appendChild(d1); gLines.appendChild(d2);
      var hit = svgEl("line", { class: "sd-hit", x1: x1, y1: y1, x2: x2, y2: y2 });
      (function (tr, r) {
        hit.addEventListener("mouseenter", function () { onEnter(tr, r); });
        hit.addEventListener("mouseleave", onLeave);
        hit.addEventListener("click", function (e) { e.stopPropagation(); onClick(tr, r); });
      })(tr, r);
      svg._gHit.appendChild(hit);
      lineRefs[tr].push({ line: line, d1: d1, d2: d2, v: p });
    }
  }

  function clearOn() {
    TREATS.forEach(function (tr) {
      (lineRefs[tr] || []).forEach(function (o) {
        if (!o) return;
        o.line.classList.remove("sd-on"); o.d1.classList.remove("sd-on"); o.d2.classList.remove("sd-on");
      });
      (chipEls[tr] || []).forEach(function (c) { c.classList.remove("sd-active"); });
    });
  }

  function applyHL(tr, r) {
    var ref = lineRefs[tr] && lineRefs[tr][r];
    if (!ref) return;
    panel.classList.add("sd-hl");
    clearOn();
    ref.line.classList.add("sd-on"); ref.d1.classList.add("sd-on"); ref.d2.classList.add("sd-on");
    gLines.appendChild(ref.line); gLines.appendChild(ref.d1); gLines.appendChild(ref.d2);
    chipEls[tr][r].classList.add("sd-active");
    var d = ref.v[1] - ref.v[0], sign = d >= 0 ? "+" : "−";
    readout.innerHTML =
      '<span class="sd-swatch" style="background:' + COL[tr] + '"></span>' +
      "<b>" + tr + " rep. " + (r + 1) + "</b> &nbsp; " +
      fmt(ref.v[0]) + " → " + fmt(ref.v[1]) +
      ' &nbsp; <span class="sd-muted">(Δ ' + sign + fmt(Math.abs(d)) + ")</span>";
  }

  function clearHL() {
    panel.classList.remove("sd-hl");
    clearOn();
    readout.innerHTML = '<span class="sd-muted">Hover a line or a chip.</span>';
  }

  function onEnter(tr, r) { if (!pinned) applyHL(tr, r); }
  function onLeave() { if (pinned) applyHL(pinned.tr, pinned.r); else clearHL(); }
  function onClick(tr, r) {
    if (pinned && pinned.tr === tr && pinned.r === r) { pinned = null; clearHL(); }
    else { pinned = { tr: tr, r: r }; applyHL(tr, r); }
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
    if (rec.oo != null) stats.push("O vs O −log₁₀p = " + rec.oo);
    if (rec.ob != null) stats.push("O vs B −log₁₀p = " + rec.ob);
    if (rec.bb != null) stats.push("B vs B −log₁₀p = " + rec.bb);
    if (rec.ad != null) stats.push("adapted (O vs O) = " + rec.ad);
    if (rec.sig) stats.push("sig: " + rec.sig);
    panel.querySelector(".sd-stats").innerHTML = stats.join("<br>");
    clearData();
    TREATS.forEach(function (tr) { drawTreat(tr, rec[tr] || []); });
    clearHL();
    panel.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  // showEl: called from the ggiraph point's inline onclick as showEl(this);
  // the SNP key lives in the element's data-id attribute.
  function showEl(el) { if (el && el.getAttribute) show(el.getAttribute("data-id")); }

  window.SNPDetail = { show: show, showEl: showEl };
})();
