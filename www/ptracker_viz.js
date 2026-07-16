// Data-driven visualization for ProteoformTracker's real Shiny app.
// Two entry points: PT.renderSection1(payload) and PT.renderSection2(payload),
// called from server-injected <script> blocks after each "Run analysis" click.
// All fragment-level tier coloring (common/partial/unique) and propensity
// scores are computed server-side in R and arrive ready-to-draw here --
// this file only does cheap client-side work (exon-tier comparison across a
// handful of proteoforms, propensity-threshold filtering, pixel math for
// zoom/pan/click), so there is no repeat of the earlier prototype's
// per-position-recompute performance bug.
window.PT = (function () {
  const TIER_COLOR = { common: "#eda100", partial: "#8952e0", unique: "#2a78d6", neutral: "#8a8a86" };
  const PX0 = 90, PX1 = 620, ROW_H = 84;
  const FILTER_LABELS = ["All fragments", "Elevated (score>1)", "High (score>=4)", "Very high (score>=9)"];

  function passesFilter(score, level) {
    if (level === 0) return true;
    if (level === 1) return score > 1;
    if (level === 2) return score >= 4;
    return score >= 9;
  }
  function xScale(pos, vs, ve) { return PX0 + (pos - vs) / (ve - vs) * (PX1 - PX0); }
  function posFromPx(px, vs, ve) { return vs + (px - PX0) / (PX1 - PX0) * (ve - vs); }
  function clampView(vs, ve, total) {
    let w = ve - vs;
    if (w < 10) w = 10;
    if (w > total - 1) w = total - 1;
    if (vs < 1) { vs = 1; ve = 1 + w; }
    if (ve > total) { ve = total; vs = total - w; }
    return [vs, ve];
  }

  // Minimum axis-units-per-pixel resolution (inverse of zoom) below which
  // individual ticks are far enough apart to hover distinctly -- below this
  // zoom level ticks would overlap and hovering would be ambiguous/noisy,
  // so the hover tooltip stays off until the user has zoomed in enough.
  const MIN_PX_PER_UNIT_FOR_HOVER = 2.2;
  // Wider than the hover threshold -- a hoverable tick just needs to be
  // separable from its neighbor, but a legible monospace letter needs
  // meaningfully more width per residue.
  const MIN_PX_PER_UNIT_FOR_AA = 7;

  function zoomAt(state, factor, centerUnit) {
    const total = state.axisLength;
    let w = (state.ve - state.vs) * factor;
    const minW = Math.max(6, total * 0.008);
    if (w < minW) w = minW;
    if (w > total - 1) w = total - 1;
    const oldW = state.ve - state.vs;
    const frac = oldW > 0 ? (centerUnit - state.vs) / oldW : 0.5;
    let vs = centerUnit - frac * w;
    let ve = vs + w;
    const clamped = clampView(vs, ve, total);
    state.vs = clamped[0]; state.ve = clamped[1];
  }

  function zoomControlsHtml(idPrefix) {
    return `<div class="pt-zoom-controls">` +
      `<span style="font-size:11.5px;color:#666;">Zoom:</span>` +
      `<button type="button" id="${idPrefix}-zoom-out" title="Zoom out">&minus;</button>` +
      `<button type="button" id="${idPrefix}-zoom-in" title="Zoom in">+</button>` +
      `<button type="button" class="pt-zoom-reset" id="${idPrefix}-zoom-reset" title="Reset zoom">Reset</button>` +
      `<span class="pt-zoom-range" id="${idPrefix}-zoom-range"></span>` +
      `</div>`;
  }

  function updateZoomRangeLabel(idPrefix, state) {
    const el = document.getElementById(idPrefix + "-zoom-range");
    if (!el) return;
    const w = Math.round(state.ve - state.vs);
    const pct = Math.round((w / state.axisLength) * 100);
    el.textContent = `${Math.round(state.vs)}-${Math.round(state.ve)} of ${Math.round(state.axisLength)} (${pct}%)`;
  }

  function wireZoomButtons(idPrefix, state, redraw) {
    const mid = (state.vs + state.ve) / 2;
    const outBtn = document.getElementById(idPrefix + "-zoom-out");
    const inBtn = document.getElementById(idPrefix + "-zoom-in");
    const resetBtn = document.getElementById(idPrefix + "-zoom-reset");
    if (outBtn) outBtn.onclick = () => { zoomAt(state, 1.5, (state.vs + state.ve) / 2); updateZoomRangeLabel(idPrefix, state); redraw(); };
    if (inBtn) inBtn.onclick = () => { zoomAt(state, 1 / 1.5, (state.vs + state.ve) / 2); updateZoomRangeLabel(idPrefix, state); redraw(); };
    if (resetBtn) resetBtn.onclick = () => { state.vs = 1; state.ve = state.axisLength; updateZoomRangeLabel(idPrefix, state); redraw(); };
    updateZoomRangeLabel(idPrefix, state);
  }

  function ensureTooltip() {
    let el = document.getElementById("pt-hover-tooltip");
    if (!el) {
      el = document.createElement("div");
      el.id = "pt-hover-tooltip";
      document.body.appendChild(el);
    }
    return el;
  }
  function showTooltip(x, y, html) {
    const el = ensureTooltip();
    el.innerHTML = html;
    el.style.display = "block";
    el.style.left = (x + 14) + "px";
    el.style.top = (y + 14) + "px";
  }
  function hideTooltip() {
    const el = document.getElementById("pt-hover-tooltip");
    if (el) el.style.display = "none";
  }

  // Shared pan/zoom wiring for any of this app's SVG panels: mouse-wheel
  // (and pinch-trackpad, which browsers report as wheel+ctrlKey) zooms
  // around the cursor; drag pans; the cursor itself is left at its default
  // arrow throughout (no grab/grabbing hand) per explicit user preference.
  // opts.onHover(px, py, pxPerUnit)/opts.onLeave() let callers layer
  // hover-specific behavior (e.g. the ion tooltip) on top.
  function attachPanZoom(svgEl, state, redraw, opts) {
    opts = opts || {};
    let dragging = false, dragStartX = 0, dragStart = [1, 1];
    svgEl.onmousedown = e => { dragging = true; dragStartX = e.clientX; dragStart = [state.vs, state.ve]; if (opts.onLeave) opts.onLeave(); };
    window.addEventListener("mousemove", e => {
      if (!dragging) return;
      const rect = svgEl.getBoundingClientRect();
      const scaleX = 640 / rect.width;
      const dx = (e.clientX - dragStartX) * scaleX / (PX1 - PX0) * (dragStart[1] - dragStart[0]);
      state.vs = dragStart[0] - dx; state.ve = dragStart[1] - dx;
      [state.vs, state.ve] = clampView(state.vs, state.ve, state.axisLength);
      if (opts.onZoomChange) opts.onZoomChange();
      redraw();
    });
    window.addEventListener("mouseup", () => { dragging = false; });
    svgEl.addEventListener("wheel", e => {
      e.preventDefault();
      const rect = svgEl.getBoundingClientRect();
      const vb = svgEl.viewBox.baseVal;
      const px = (e.clientX - rect.left) * vb.width / rect.width;
      const centerUnit = posFromPx(px, state.vs, state.ve);
      const factor = e.deltaY > 0 ? 1.15 : 1 / 1.15;
      zoomAt(state, factor, centerUnit);
      if (opts.onZoomChange) opts.onZoomChange();
      redraw();
      if (opts.onHover) {
        const py = (e.clientY - rect.top) * vb.height / rect.height;
        const pxPerUnit = (PX1 - PX0) / (state.ve - state.vs);
        opts.onHover(px, py, pxPerUnit);
      }
    }, { passive: false });
    svgEl.addEventListener("mousemove", e => {
      if (dragging || !opts.onHover) return;
      const rect = svgEl.getBoundingClientRect();
      const vb = svgEl.viewBox.baseVal;
      const px = (e.clientX - rect.left) * vb.width / rect.width;
      const py = (e.clientY - rect.top) * vb.height / rect.height;
      const pxPerUnit = (PX1 - PX0) / (state.ve - state.vs);
      opts.onHover(px, py, pxPerUnit, e.clientX, e.clientY);
    });
    svgEl.addEventListener("mouseleave", () => { if (opts.onLeave) opts.onLeave(); });
  }

  function legendHtml(showNeutral) {
    let h = '<div style="display:flex;gap:16px;flex-wrap:wrap;padding:4px 0;font-size:12px;color:#555;">';
    h += '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;border-radius:2px;display:inline-block;"></span>common</span>';
    h += '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#8952e0;border-radius:2px;display:inline-block;"></span>partial</span>';
    h += '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#2a78d6;border-radius:2px;display:inline-block;"></span>unique</span>';
    if (showNeutral) h += '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#8a8a86;border-radius:2px;display:inline-block;"></span>neutral (nothing to compare)</span>';
    h += '<span style="color:#888;">tick opacity = fragmentation propensity</span>';
    h += "</div>";
    return h;
  }

  function filterButtonsHtml(idPrefix, active) {
    let h = '<div style="display:flex;gap:6px;flex-wrap:wrap;align-items:center;padding:4px 0;">';
    h += '<span style="font-size:12px;color:#666;">Fragment filter:</span>';
    FILTER_LABELS.forEach((lbl, i) => {
      const on = i === active;
      h += `<button data-filter="${i}" class="pt-filter-btn" style="font-size:11.5px;padding:4px 9px;border-radius:6px;border:1px solid ${on ? '#0ca30c' : '#ccc'};background:${on ? 'rgba(12,163,12,.15)' : '#f5f5f5'};cursor:pointer;font-weight:${on ? '600' : '400'};" id="${idPrefix}-f${i}">${lbl}</button>`;
    });
    h += '<span id="' + idPrefix + '-fcount" style="font-size:11.5px;color:#888;margin-left:6px;"></span>';
    h += "</div>";
    return h;
  }

  function drawMS1(svgEl, entries) {
    if (!entries.length) { svgEl.innerHTML = ""; return; }
    let allMz = [];
    entries.forEach(e => e.env.forEach(p => allMz.push(p.mz)));
    const mn = Math.min(...allMz) * 0.95, mx = Math.max(...allMz) * 1.05;
    const mzx = mz => 40 + (mz - mn) / (mx - mn) * 560;
    const baseY = 86;
    let els = `<line x1="40" y1="${baseY}" x2="600" y2="${baseY}" stroke="#ccc" stroke-width="1"/>`;
    [0, 0.5, 1].forEach(v => {
      const y = baseY - v * 70;
      els += `<line x1="36" y1="${y}" x2="600" y2="${y}" stroke="#ddd" stroke-width="0.5" stroke-dasharray="2,2"/><text x="32" y="${y + 3}" font-size="9" fill="#666" text-anchor="end">${v.toFixed(1)}</text>`;
    });
    els += `<text x="10" y="51" font-size="9" fill="#666" text-anchor="middle" transform="rotate(-90, 10, 51)">Relative intensity</text>`;
    entries.forEach((e, i) => {
      e.env.forEach(p => {
        const x = mzx(p.mz), h = p.rel * 70;
        els += `<line x1="${x.toFixed(1)}" y1="${baseY}" x2="${x.toFixed(1)}" y2="${(baseY - h).toFixed(1)}" stroke="${e.color}" stroke-width="2" stroke-opacity="0.85"/>`;
      });
      els += `<text x="44" y="${12 + i * 12}" font-size="10" fill="${e.color}" font-weight="600">${e.label}: ${e.mass.toFixed(1)} Da</text>`;
    });
    const ticks = [Math.round(mn), Math.round((mn + mx) / 2), Math.round(mx)];
    const anchors = ["start", "middle", "end"];
    ticks.forEach((t, i) => { els += `<text x="${mzx(t).toFixed(1)}" y="98" font-size="9" fill="#666" text-anchor="${anchors[i]}">${t} m/z</text>`; });
    svgEl.setAttribute("viewBox", "0 0 640 130");
    svgEl.innerHTML = els;
  }

  function exonTier(block, others) {
    if (others.length === 0) return "neutral";
    let mc = 0;
    others.forEach(o => { if (o.exon_blocks.some(b => b.start === block.start && b.end === block.end)) mc++; });
    if (mc === others.length) return "common";
    if (mc === 0) return "unique";
    return "partial";
  }

  function renderLadderRows(svgEl, proteoforms, axisLength, state) {
    const height = Math.max(1, proteoforms.length) * ROW_H + 6;
    svgEl.setAttribute("viewBox", `0 0 640 ${height}`);
    let els = "";
    let visibleCount = 0, totalCount = 0;
    const ticks = [];
    // Amino-acid letters only pay off once each residue has enough width to
    // read a monospace character without overlapping its neighbors -- below
    // that they'd just be an illegible smear, so they stay off until the
    // user has zoomed in past MIN_PX_PER_UNIT_FOR_AA.
    const pxPerUnit = (PX1 - PX0) / (state.ve - state.vs);
    const showAA = pxPerUnit >= MIN_PX_PER_UNIT_FOR_AA;
    proteoforms.forEach((pf, ri) => {
      const top = ri * ROW_H + 4;
      const others = proteoforms.filter((_, j) => j !== ri);
      els += `<text x="4" y="${top + 9}" font-size="10" font-weight="600" fill="${pf.color}">${pf.label}</text>`;
      els += `<text x="4" y="${top + 20}" font-size="9" fill="#666">${pf.len} aa, ${pf.mass.toFixed(1)} Da</text>`;
      els += `<line x1="${PX0}" y1="${top + 24}" x2="${PX1}" y2="${top + 24}" stroke="#ddd" stroke-width="1"/>`;

      pf.ptms.forEach(ptm => {
        const x = xScale(ptm.axis_pos, state.vs, state.ve);
        if (x < PX0 - 2 || x > PX1 + 2) return;
        els += `<line x1="${x.toFixed(1)}" y1="${top + 24}" x2="${x.toFixed(1)}" y2="${top + 13}" stroke="#e0433d" stroke-width="1.5"/>`;
        els += `<circle cx="${x.toFixed(1)}" cy="${top + 11}" r="3" fill="#e0433d"/>`;
        els += `<text x="${x.toFixed(1)}" y="${top + 9}" font-size="7" text-anchor="middle" fill="#555">${ptm.name}</text>`;
      });

      pf.exon_blocks.forEach(block => {
        if (block.end < state.vs || block.start > state.ve) return;
        const x0 = Math.max(xScale(block.start, state.vs, state.ve), PX0);
        const x1 = Math.min(xScale(block.end, state.vs, state.ve), PX1);
        if (x1 <= x0) return;
        const t = exonTier(block, others);
        els += `<rect x="${x0.toFixed(1)}" y="${top + 28}" width="${(x1 - x0).toFixed(1)}" height="8" fill="${TIER_COLOR[t]}" rx="1.5"/>`;
      });

      if (showAA && pf.sequence) {
        // axis_pos[p-1] is the p-th peptide-bond position (same convention
        // the b/y ticks below use); a residue's own letter sits at its
        // preceding bond's position, so this reuses that array directly
        // instead of needing a separate per-residue axis mapping.
        for (let p = 1; p <= pf.sequence.length; p++) {
          const axisP = p <= pf.axis_pos.length ? pf.axis_pos[p - 1] : pf.axis_pos[pf.axis_pos.length - 1] + 1;
          if (axisP < 0 || axisP < state.vs || axisP > state.ve) continue;
          const x = xScale(axisP, state.vs, state.ve);
          if (x < PX0 - 2 || x > PX1 + 2) continue;
          els += `<text x="${x.toFixed(1)}" y="${top + 45}" font-size="10" font-family="monospace" text-anchor="middle" fill="#333">${pf.sequence[p - 1]}</text>`;
        }
      }

      for (let p = 1; p <= pf.b_mass.length; p++) {
        const axisP = pf.axis_pos[p - 1];
        if (axisP < 0 || axisP < state.vs || axisP > state.ve) continue;
        const score = pf.propensity[p - 1];
        totalCount++;
        if (!passesFilter(score, state.filter)) continue;
        visibleCount++;
        const x = xScale(axisP, state.vs, state.ve);
        if (x < PX0 - 2 || x > PX1 + 2) continue;
        const op = state.filter === 0 ? (0.22 + 0.78 * Math.max(0, Math.min(1, (score - 1) / 9.5))) : 0.95;
        const w = state.filter === 0 ? 1.2 : state.filter === 1 ? 1.6 : state.filter === 2 ? 2.2 : 2.8;
        els += `<line data-row="${ri}" data-p="${p}" data-ion="b" x1="${x.toFixed(2)}" y1="${top + 50}" x2="${x.toFixed(2)}" y2="${top + 62}" stroke-width="${w}" stroke="${TIER_COLOR[pf.tier_b[p - 1]]}" stroke-opacity="${op.toFixed(2)}"/>`;
        els += `<line data-row="${ri}" data-p="${p}" data-ion="y" x1="${x.toFixed(2)}" y1="${top + 66}" x2="${x.toFixed(2)}" y2="${top + 78}" stroke-width="${w}" stroke="${TIER_COLOR[pf.tier_y[p - 1]]}" stroke-opacity="${op.toFixed(2)}"/>`;
        ticks.push({ ri, p, axisP, xPx: x });
      }
      els += `<text x="4" y="${top + 56}" font-size="9" fill="#999">b</text>`;
      els += `<text x="4" y="${top + 74}" font-size="9" fill="#999">y</text>`;
    });
    svgEl.innerHTML = els;
    state.ticks = ticks;
    return { visibleCount, totalCount };
  }

  // Finds the b/y tick nearest the cursor, restricted to the row/ion band
  // the cursor is currently over, within a small pixel radius -- feeds the
  // hover tooltip below. Returns null when nothing is close enough or the
  // cursor isn't over a b/y band at all.
  function findNearestTick(state, proteoforms, px, py) {
    const ri = Math.floor(py / ROW_H);
    if (ri < 0 || ri >= proteoforms.length) return null;
    const top = ri * ROW_H + 4;
    const wr = py - top;
    if (wr < 46 || wr > 78) return null;
    const ion = wr < 62 ? "b" : "y";
    let best = null, bestD = 5;
    (state.ticks || []).forEach(t => {
      if (t.ri !== ri) return;
      const d = Math.abs(t.xPx - px);
      if (d < bestD) { bestD = d; best = t; }
    });
    if (!best) return null;
    return { ri, ion, p: best.p };
  }

  function tickTooltipHtml(pf, ion, localP) {
    const mass = ion === "b" ? pf.b_mass[localP - 1] : pf.y_mass[localP - 1];
    const tier = ion === "b" ? pf.tier_b[localP - 1] : pf.tier_y[localP - 1];
    const label = ion === "b" ? "b" + localP : "y" + (pf.len - localP);
    return `<strong>${pf.label} ${label}</strong> (residue ${localP})<br/>mass ${mass.toFixed(1)} Da &middot; propensity ${pf.propensity[localP - 1]} &middot; tier: <strong>${tier}</strong>`;
  }

  function wireInteraction(svgEl, state, proteoforms, redraw, hintEl) {
    function onHover(px, py, pxPerUnit, clientX, clientY) {
      if (pxPerUnit < MIN_PX_PER_UNIT_FOR_HOVER) {
        hideTooltip();
        if (hintEl) hintEl.textContent = "Zoom in further to hover individual b/y ions.";
        return;
      }
      const hit = findNearestTick(state, proteoforms, px, py);
      if (!hit || clientX === undefined) { hideTooltip(); if (hintEl) hintEl.textContent = ""; return; }
      const pf = proteoforms[hit.ri];
      showTooltip(clientX, clientY, tickTooltipHtml(pf, hit.ion, hit.p));
      if (hintEl) hintEl.textContent = "";
    }
    function onLeave() { hideTooltip(); if (hintEl) hintEl.textContent = ""; }
    attachPanZoom(svgEl, state, redraw, { onHover, onLeave, onZoomChange: () => hideTooltip() });
  }

  function renderSection1(payload) {
    const ms1El = document.getElementById("s1-ms1");
    const ladderEl = document.getElementById("s1-ladder");
    const legendEl = document.getElementById("s1-legend");
    const filterEl = document.getElementById("s1-filters");
    const zoomEl = document.getElementById("s1-zoom");
    const infoEl = document.getElementById("s1-info");
    if (!ms1El || !ladderEl) return;

    const proteoforms = payload.proteoforms;
    const state = { vs: 1, ve: payload.axis_length, filter: 0, axisLength: payload.axis_length };

    drawMS1(ms1El, proteoforms.map(p => ({ label: p.label, color: p.color, mass: p.mass, env: p.env || [] })));
    legendEl.innerHTML = legendHtml(proteoforms.length < 2);
    filterEl.innerHTML = filterButtonsHtml("s1", state.filter);
    if (zoomEl) zoomEl.innerHTML = zoomControlsHtml("s1");

    function redraw() {
      const counts = renderLadderRows(ladderEl, proteoforms, state.axisLength, state);
      const cEl = document.getElementById("s1-fcount");
      if (cEl) cEl.textContent = counts.visibleCount + " of " + counts.totalCount + " bonds shown";
      updateZoomRangeLabel("s1", state);
    }
    redraw();
    wireInteraction(ladderEl, state, proteoforms, redraw, infoEl);
    if (zoomEl) wireZoomButtons("s1", state, redraw);
    function rewireFilterButtons() {
      FILTER_LABELS.forEach((_, i) => {
        const btn = document.getElementById("s1-f" + i);
        if (btn) btn.onclick = () => { state.filter = i; filterEl.innerHTML = filterButtonsHtml("s1", i); redraw(); rewireFilterButtons(); };
      });
    }
    rewireFilterButtons();
  }

  function renderSection2(payload) {
    const ms1El = document.getElementById("s2-ms1");
    const ladderEl = document.getElementById("s2-ladder");
    const legendEl = document.getElementById("s2-legend");
    const filterEl = document.getElementById("s2-filters");
    const zoomEl = document.getElementById("s2-zoom");
    const infoEl = document.getElementById("s2-info");
    if (!ms1El || !ladderEl) return;

    const t = payload.target;
    const targetColor = "#1f8a70";
    const entries = [{ label: "TARGET: " + t.id, color: targetColor, mass: t.mass, env: t.env || [] }];
    payload.confounders.forEach(c => entries.push({ label: c.id + " (confounder)", color: "#c0392b", mass: c.mass, env: c.env.map(e => ({ mz: e.mz, rel: e.rel })) }));
    drawMS1(ms1El, entries);

    const proteoform = {
      label: t.id, color: targetColor, mass: t.mass, len: t.len, sequence: t.sequence,
      ptms: t.ptms.map(p => Object.assign({}, p, { axis_pos: p.pos })),
      exon_blocks: t.exon_blocks.map(b => ({ start: b.start, end: b.end })),
      b_mass: t.b_mass, y_mass: t.y_mass, propensity: t.propensity,
      axis_pos: t.b_mass.map((_, i) => i + 1),
      tier_b: t.tier_b, tier_y: t.tier_y
    };
    const state = { vs: 1, ve: t.len, filter: 0, axisLength: t.len };
    legendEl.innerHTML = legendHtml(payload.confounders.length === 0) +
      `<div style="font-size:12px;color:#666;">Window: +/-${payload.window_da} Da at best charge state ${payload.best_charge_state}. ${payload.confounders.length} real confounder(s) found.</div>`;
    filterEl.innerHTML = filterButtonsHtml("s2", state.filter);
    if (zoomEl) zoomEl.innerHTML = zoomControlsHtml("s2");

    function redraw() {
      const counts = renderLadderRows(ladderEl, [proteoform], state.axisLength, state);
      const cEl = document.getElementById("s2-fcount");
      if (cEl) cEl.textContent = counts.visibleCount + " of " + counts.totalCount + " bonds shown";
      updateZoomRangeLabel("s2", state);
    }
    redraw();
    wireInteraction(ladderEl, state, [proteoform], redraw, infoEl);
    if (zoomEl) wireZoomButtons("s2", state, redraw);
    function rewireFilterButtons() {
      FILTER_LABELS.forEach((_, i) => {
        const btn = document.getElementById("s2-f" + i);
        if (btn) btn.onclick = () => { state.filter = i; filterEl.innerHTML = filterButtonsHtml("s2", i); redraw(); rewireFilterButtons(); };
      });
    }
    rewireFilterButtons();
  }

  // Multi-row exon-presence comparison for Option 2 (novel FASTA sequence
  // vs. one or more user-picked known isoforms), rendered BEFORE any
  // MS1/MS2 analysis -- this is purely about "does this exon exist here",
  // independent of translation/ORF choice. Each row draws a thin connecting
  // line across its FULL width first (the intron backbone), then draws its
  // own exon blocks on top at their own true genomic-derived position --
  // a stretch a row doesn't cover simply shows the bare line through it,
  // which is the "placed gap" alignment. Two tracks' exons that overlap
  // but don't share identical boundaries render at genuinely overlapping
  // (not merged, not force-identical) positions, same as a real
  // genome-browser multi-track view. Tier color (common/partial/unique)
  // comes from the server-computed shared_count (real genomic overlap with
  // the other compared rows), not axis-position matching.
  function renderExonAlignment(svgEl, payload, state) {
    if (!payload) { svgEl.innerHTML = '<text x="10" y="20" font-size="11" fill="#999">No known transcript exon data available for this comparison.</text>'; svgEl.setAttribute("viewBox", "0 0 640 40"); return; }
    const vs = state ? state.vs : 0, ve = state ? state.ve : payload.axis_length;
    const rowH = 32, topPad = 24;
    const n = payload.n_tracks || payload.tracks.length;
    const height = payload.tracks.length * rowH + 24 + (topPad - 10);
    const xs = pos => PX0 + (pos - vs) / (ve - vs) * (PX1 - PX0);
    const rowColors = ["#1f8a70", "#c0392b", "#7a5cff", "#d9730d", "#2a78d6", "#5b6470"];

    // 5'/3' end labels -- the axis always runs 5' (left) to 3' (right)
    // regardless of the gene's genomic strand (minus-strand genes are
    // already flipped onto this convention before reaching this payload).
    let els = `<text x="${PX0}" y="10" font-size="10" font-weight="600" fill="#333" text-anchor="middle">5'</text>` +
      `<text x="${PX1}" y="10" font-size="10" font-weight="600" fill="#333" text-anchor="middle">3'</text>`;
    payload.tracks.forEach((t, ti) => {
      const top = topPad + ti * rowH;
      const color = rowColors[ti % rowColors.length];
      els += `<text x="4" y="${top + 10}" font-size="10" font-weight="600" fill="${color}">${t.label}</text>`;
      els += `<line x1="${PX0}" y1="${top + 16}" x2="${PX1}" y2="${top + 16}" stroke="#bbb" stroke-width="1.5"/>`;
      t.blocks.forEach(b => {
        if (b.end < vs || b.start > ve) return;
        const tierColor = b.shared_count === n ? TIER_COLOR.common : b.shared_count === 1 ? TIER_COLOR.unique : TIER_COLOR.partial;
        const x0 = Math.max(xs(b.start), PX0), x1 = Math.min(xs(b.end), PX1);
        if (x1 <= x0) return;
        if (b.coding === undefined) {
          // coding status unknown for this track (e.g. novel sequence
          // before any TransDecoder candidate is selected) -- render solid,
          // same as before this distinction existed.
          els += `<rect x="${x0.toFixed(1)}" y="${top + 10}" width="${Math.max(1, x1 - x0).toFixed(1)}" height="12" fill="${tierColor}" rx="1.5"/>`;
        } else {
          // Light/translucent base layer spans the whole block (reads as
          // non-coding/UTR); a full-opacity overlay is drawn only over the
          // real coding (CDS) sub-range(s), so a "hybrid" exon that's part
          // UTR, part CDS shows both within the same block.
          els += `<rect x="${x0.toFixed(1)}" y="${top + 10}" width="${Math.max(1, x1 - x0).toFixed(1)}" height="12" fill="${tierColor}" fill-opacity="0.3" rx="1.5"/>`;
          b.coding.forEach(c => {
            if (c.end < vs || c.start > ve) return;
            const cx0 = Math.max(xs(c.start), PX0), cx1 = Math.min(xs(c.end), PX1);
            if (cx1 <= cx0) return;
            els += `<rect x="${cx0.toFixed(1)}" y="${top + 10}" width="${Math.max(1, cx1 - cx0).toFixed(1)}" height="12" fill="${tierColor}"/>`;
          });
        }
      });
    });

    els += `<text x="${PX0}" y="${height - 4}" font-size="9" fill="#888">${payload.seqname}, ${payload.strand} strand -- exon widths shown to scale relative to each other; introns shown as connecting lines, not to scale</text>`;

    svgEl.setAttribute("viewBox", `0 0 640 ${height}`);
    svgEl.innerHTML = els;
  }

  function exonAlignmentLegendHtml() {
    return '<div style="display:flex;flex-direction:column;gap:4px;padding:4px 0;font-size:12px;color:#555;">' +
      '<div style="display:flex;gap:16px;flex-wrap:wrap;">' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;border-radius:2px;display:inline-block;"></span>present in all compared</span>' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#8952e0;border-radius:2px;display:inline-block;"></span>present in some, not all</span>' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#2a78d6;border-radius:2px;display:inline-block;"></span>present in only one</span>' +
      "</div>" +
      '<div style="display:flex;gap:16px;flex-wrap:wrap;">' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;border-radius:2px;display:inline-block;"></span>coding (CDS)</span>' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;opacity:0.3;border-radius:2px;display:inline-block;"></span>non-coding (UTR)</span>' +
      "</div>" +
      "</div>";
  }

  return {
    renderSection1, renderSection2, renderExonAlignment, exonAlignmentLegendHtml,
    zoomControlsHtml, wireZoomButtons, updateZoomRangeLabel, attachPanZoom
  };
})();

if (window.Shiny) {
  // Persisted across re-renders (isoform/ORF selection changes) so the
  // user's zoom/pan doesn't reset on every tweak -- only when the
  // underlying axis actually changes length does the view snap back to
  // full-width, since an old vs/ve pair could otherwise fall outside a
  // shorter new axis.
  let exonAlignState = null;
  Shiny.addCustomMessageHandler("pt_render_exon_alignment", function (msg) {
    const svgEl = document.getElementById("fasta-exon-align");
    const legendEl = document.getElementById("fasta-exon-align-legend");
    const zoomEl = document.getElementById("fasta-exon-zoom");
    if (!svgEl) return;
    const payload = msg.has_data ? msg.payload : null;
    if (legendEl) legendEl.innerHTML = payload ? PT.exonAlignmentLegendHtml() : "";
    if (!payload) { PT.renderExonAlignment(svgEl, null); if (zoomEl) zoomEl.innerHTML = ""; exonAlignState = null; return; }

    if (!exonAlignState || exonAlignState.axisLength !== payload.axis_length) {
      exonAlignState = { vs: 0, ve: payload.axis_length, axisLength: payload.axis_length };
    }
    const state = exonAlignState;

    function redraw() {
      PT.renderExonAlignment(svgEl, payload, state);
      PT.updateZoomRangeLabel("fasta-exon", state);
    }
    redraw();
    if (zoomEl) {
      zoomEl.innerHTML = PT.zoomControlsHtml("fasta-exon");
      PT.wireZoomButtons("fasta-exon", state, redraw);
    }
    PT.attachPanZoom(svgEl, state, redraw, {});
  });
}
