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
  // Distinct qualitative colors for confounder MS1 curves (renderSection2),
  // cycled if there are more confounders than colors. Deliberately avoids
  // the target's own green (#1f8a70) and the TIER_COLOR hues (used
  // elsewhere for a different meaning -- fragment-tier common/partial/
  // unique) so nothing reads as double-booked across the two charts.
  // Before this, every confounder shared one uniform red -- harmless when
  // confounders were mass-window-only and rarely landed on the SAME m/z as
  // the target, but the m/z-collision search (search_confounding_
  // proteins_combined()) specifically curates confounders that DO share a
  // peak with the target, so a uniform color made every confounder's curve
  // indistinguishable from every other's at exactly the positions that
  // matter most.
  const CONFOUNDER_PALETTE = ["#c0392b", "#e67e22", "#c2185b", "#6f42c1", "#2471a3", "#795548", "#607d8b", "#b8860b"];
  // Colors for exon-alignment "highlight" boxes (e.g. Option 3's rMATS
  // differential exon(s)) -- drawn as a dashed outline ON TOP of a block's
  // normal common/partial/unique tier fill, cycling if there's more than
  // one highlight (MXE has two: its 1st and 2nd mutually exclusive exons).
  const HIGHLIGHT_COLORS = ["#111111", "#0b5fff", "#c0392b", "#0ca35c"];
  // Tier cutpoints, calibrated separately per scoring mode -- the two
  // models are fit fully independently and must never share thresholds,
  // even though both now output a genuine 0-1 probability (see
  // R/fragmentation_propensity.R vs R/fragmentation_propensity_rf.R).
  // "glm": a single joint logistic regression (data/fragmentation_
  // propensity_glm.rds), replacing an earlier version assembled from
  // several separately-calibrated pieces multiplied together -- that
  // assembly was never validated as a whole and turned out to score much
  // worse (cross-study holdout AUROC ~0.64) than one joint fit on the same
  // features (~0.78-0.82). Thresholds selected via the same fold-
  // enrichment sweep as before: ~3.1x/~4.7x over the ~6.3% pooled baseline
  // match rate. "rf": length-free ranking mode, fold-enrichment ~2.5x/~3.5x
  // (scripts/build_propensity_rf_model.R). A third "Very high" tier existed
  // at both scales but was dropped -- confirmed directly it almost never
  // fires for realistic-length proteoforms in either mode, so it was
  // mostly just an extra click that led to an empty ladder.
  const TIER_THRESHOLDS = {
    glm: { elevated: 0.10, high: 0.20, opNorm: [0, 0.4],
           labels: ["All fragments", "Elevated (score>0.10)", "High (score>=0.20)"] },
    rf: { elevated: 0.08, high: 0.15, opNorm: [0, 0.3],
          labels: ["All fragments", "Elevated (score>0.08)", "High (score>=0.15)"] }
  };
  function tierThresholds(mode) { return TIER_THRESHOLDS[mode] || TIER_THRESHOLDS.glm; }
  // isotope-panel computation gate, mirroring R's GLM_ELEVATED_THRESHOLD /
  // RF_ELEVATED_THRESHOLD (R/fragmentation_propensity.R /
  // R/fragmentation_propensity_rf.R) -- both modes now gate on their own
  // Elevated cutoff rather than GLM's old baseline-relative "1".
  const GLM_ISOTOPE_GATE_JS = 0.10;
  const RF_ISOTOPE_GATE_JS = 0.08;

  function passesFilter(score, level, mode) {
    const t = tierThresholds(mode);
    if (level === 0) return true;
    if (level === 1) return score > t.elevated;
    return score >= t.high;
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

  // minWidth overrides the default "at least 0.8% of the full span, floor
  // 6 units" rule -- the residue axis this was written for has a natural
  // minimum (1 residue) way looser than what the MS1 chart needs (real
  // isotope peaks can be a small fraction of a single m/z unit apart at
  // high charge states), so callers with a tighter natural resolution pass
  // their own floor instead of inheriting one tuned for a different domain.
  function zoomAt(state, factor, centerUnit, minWidth) {
    const total = state.axisLength;
    let w = (state.ve - state.vs) * factor;
    const minW = minWidth != null ? minWidth : Math.max(6, total * 0.008);
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

  // formatFn(state) => label string, for domains where "residue N-M of
  // total" doesn't read naturally (e.g. the MS1 chart's m/z domain wants
  // real m/z values with decimals, not an integer offset-from-minimum).
  function updateZoomRangeLabel(idPrefix, state, formatFn) {
    const el = document.getElementById(idPrefix + "-zoom-range");
    if (!el) return;
    if (formatFn) { el.textContent = formatFn(state); return; }
    const w = Math.round(state.ve - state.vs);
    const pct = Math.round((w / state.axisLength) * 100);
    el.textContent = `${Math.round(state.vs)}-${Math.round(state.ve)} of ${Math.round(state.axisLength)} (${pct}%)`;
  }

  function wireZoomButtons(idPrefix, state, redraw, minZoomWidth, formatFn) {
    const mid = (state.vs + state.ve) / 2;
    const outBtn = document.getElementById(idPrefix + "-zoom-out");
    const inBtn = document.getElementById(idPrefix + "-zoom-in");
    const resetBtn = document.getElementById(idPrefix + "-zoom-reset");
    if (outBtn) outBtn.onclick = () => { zoomAt(state, 1.5, (state.vs + state.ve) / 2, minZoomWidth); updateZoomRangeLabel(idPrefix, state, formatFn); redraw(); };
    if (inBtn) inBtn.onclick = () => { zoomAt(state, 1 / 1.5, (state.vs + state.ve) / 2, minZoomWidth); updateZoomRangeLabel(idPrefix, state, formatFn); redraw(); };
    if (resetBtn) resetBtn.onclick = () => { state.vs = 1; state.ve = state.axisLength; updateZoomRangeLabel(idPrefix, state, formatFn); redraw(); };
    updateZoomRangeLabel(idPrefix, state, formatFn);
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
  //
  // opts.signal (an AbortSignal) lets a caller that re-wires the SAME svgEl
  // repeatedly (e.g. the MS2 fragment-isotope chart, re-rendered on every
  // bond click, unlike the MS1/ladder charts which only ever wire once per
  // "Run analysis") cancel every listener from the PREVIOUS wiring in one
  // shot before attaching new ones -- without it, each click would pile up
  // another window-level mousemove/mouseup listener (never removed) plus
  // more wheel/mousemove/mouseleave listeners on svgEl itself, each closing
  // over that click's now-stale `state`/`entries`, silently multiplying
  // duplicate work and hover/drag glitches with every subsequent click.
  // Omitted entirely, addEventListener behaves exactly as before (a
  // `signal: undefined` option is equivalent to no signal at all).
  function attachPanZoom(svgEl, state, redraw, opts) {
    opts = opts || {};
    const signal = opts.signal;
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
    }, { signal });
    window.addEventListener("mouseup", () => { dragging = false; }, { signal });
    svgEl.addEventListener("wheel", e => {
      e.preventDefault();
      const rect = svgEl.getBoundingClientRect();
      const vb = svgEl.viewBox.baseVal;
      const px = (e.clientX - rect.left) * vb.width / rect.width;
      const centerUnit = posFromPx(px, state.vs, state.ve);
      const factor = e.deltaY > 0 ? 1.15 : 1 / 1.15;
      zoomAt(state, factor, centerUnit, opts.minZoomWidth);
      if (opts.onZoomChange) opts.onZoomChange();
      redraw();
      if (opts.onHover) {
        const py = (e.clientY - rect.top) * vb.height / rect.height;
        const pxPerUnit = (PX1 - PX0) / (state.ve - state.vs);
        opts.onHover(px, py, pxPerUnit);
      }
    }, { passive: false, signal });
    svgEl.addEventListener("mousemove", e => {
      if (dragging || !opts.onHover) return;
      const rect = svgEl.getBoundingClientRect();
      const vb = svgEl.viewBox.baseVal;
      const px = (e.clientX - rect.left) * vb.width / rect.width;
      const py = (e.clientY - rect.top) * vb.height / rect.height;
      const pxPerUnit = (PX1 - PX0) / (state.ve - state.vs);
      opts.onHover(px, py, pxPerUnit, e.clientX, e.clientY);
    }, { signal });
    svgEl.addEventListener("mouseleave", () => { if (opts.onLeave) opts.onLeave(); }, { signal });
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

  function filterButtonsHtml(idPrefix, active, mode) {
    let h = '<div style="display:flex;gap:6px;flex-wrap:wrap;align-items:center;padding:4px 0;">';
    h += '<span style="font-size:12px;color:#666;">Fragment filter:</span>';
    tierThresholds(mode).labels.forEach((lbl, i) => {
      const on = i === active;
      h += `<button data-filter="${i}" class="pt-filter-btn" style="font-size:11.5px;padding:4px 9px;border-radius:6px;border:1px solid ${on ? '#0ca30c' : '#ccc'};background:${on ? 'rgba(12,163,12,.15)' : '#f5f5f5'};cursor:pointer;font-weight:${on ? '600' : '400'};" id="${idPrefix}-f${i}">${lbl}</button>`;
    });
    h += '<span id="' + idPrefix + '-fcount" style="font-size:11.5px;color:#888;margin-left:6px;"></span>';
    h += "</div>";
    return h;
  }

  // Aggregate MS1 stat tiles (clean/overlapping/total charge-state peaks),
  // rendered beside the MS1 chart itself. Per-proteoform MS2 tallies are no
  // longer a separate panel here at all -- they're drawn INLINE inside the
  // ladder SVG next to each row's own title (rowBadgesSvg(), above
  // renderLadderRows()), recomputed live against whichever fragment filter
  // is currently active, so there's nothing left to keep in sync with the
  // ladder's layout. Pure presentation of numbers already computed server-
  // side (envelope_crowding_check()/search_confounding_proteins_combined()
  // -- R/viz_json.R, server.R) -- no client-side counting logic here.
  function statTile(iconHtml, value, label) {
    return `<div class="pt-stat-tile">${iconHtml}<div><div class="pt-stat-value">${value}</div><div class="pt-stat-label">${label}</div></div></div>`;
  }
  function glyphIcon(glyph, color) { return `<span class="pt-stat-icon" style="color:${color};">${glyph}</span>`; }

  function renderStatsStrip(stripEl, ms1Stats, ms1Label) {
    if (!stripEl) return;
    if (!ms1Stats) { stripEl.innerHTML = ""; return; }
    let h = `<div class="pt-stat-group-title">${ms1Label}</div>`;
    h += statTile(glyphIcon("&#10003;", "#1baf7a"), ms1Stats.clean_peaks, "Clean peaks (no overlap)");
    h += statTile(glyphIcon("&#9888;", "#e34948"), ms1Stats.crowded_peaks, "Overlapping peaks");
    h += statTile(glyphIcon("&Sigma;", "#666"), ms1Stats.total_peaks, "Total charge-state peaks");
    stripEl.innerHTML = h;
  }

  // Isotope peaks (or, for an unresolved charge state, the sample points of
  // its smooth envelope curve) are only ever a small fraction of a single
  // m/z unit apart at realistic charge states, while the chart's default
  // view spans the WHOLE charge-state envelope (hundreds to thousands of
  // m/z) -- so at the default zoomed-all-the-way-out view, one charge
  // state's entire isotope comb can be under a single pixel wide,
  // indistinguishable from the single-stick-per-charge-state rendering this
  // replaced. Computes: the full m/z domain (mzMin + axisLength) for the
  // reset view, and a data-driven minimum zoom width (a small multiple of
  // the tightest real gap between any two plotted points anywhere in this
  // chart) so "zoom in as far as the controls allow" actually lands on
  // "see individual isotope peaks clearly" for THIS data, not a fixed
  // guess that's too loose for a highly-charged small protein or
  // needlessly tight for a large one.
  function computeMS1Domain(entries) {
    let allMz = [];
    entries.forEach(e => (e.env || []).forEach(cs => cs.points.forEach(p => allMz.push(p.mz))));
    if (!allMz.length) return null;
    allMz.sort((a, b) => a - b);
    const mn = allMz[0] * 0.95, mx = allMz[allMz.length - 1] * 1.05;
    let minGap = Infinity;
    for (let i = 1; i < allMz.length; i++) {
      const g = allMz[i] - allMz[i - 1];
      if (g > 1e-9 && g < minGap) minGap = g;
    }
    if (!isFinite(minGap)) minGap = (mx - mn) * 0.01;
    return { mzMin: mn, axisLength: mx - mn, minZoomWidth: Math.max(minGap * 10, 1e-6) };
  }

  // entries: [{label, color, mass, env}], where env is an array of charge
  // states (predict_ms1_peaks()/ms1_peaks_json(), R/isotope_envelope.R +
  // R/viz_json.R): [{z, resolved, points: [{mz, rel}]}]. `resolved` is
  // decided per charge state against the instrument's own resolving power
  // (fwhm_mz() vs. real isotope peak spacing at that m/z) -- resolved
  // charge states draw as a real discrete isotope comb (one stick per
  // actual isotope peak); unresolved ones draw as the smooth envelope
  // shape those unresolved peaks would actually blur into on a real
  // spectrum, using `points` as samples along that curve instead of
  // discrete peak positions. See predict_ms1_peaks()'s doc comment for why.
  //
  // state = {vs, ve, axisLength, mzMin} (computeMS1Domain() above), same
  // shape/semantics as the ladder chart's zoom state so the SAME pan/zoom
  // machinery (attachPanZoom/zoomAt/clampView, all written against an
  // abstract [1, axisLength] "unit" domain) works unchanged here too --
  // state.vs/ve are an OFFSET from state.mzMin, not an absolute m/z value.
  function drawMS1(svgEl, entries, state) {
    if (!entries.length || !state) { svgEl.innerHTML = ""; return; }
    const mzx = mz => PX0 + (mz - state.mzMin - state.vs) / (state.ve - state.vs) * (PX1 - PX0);
    const loMz = state.mzMin + state.vs, hiMz = state.mzMin + state.ve;
    const baseY = 86;
    let els = `<line x1="${PX0}" y1="${baseY}" x2="${PX1}" y2="${baseY}" stroke="#ccc" stroke-width="1"/>`;
    [0, 0.5, 1].forEach(v => {
      const y = baseY - v * 70;
      els += `<line x1="${PX0 - 4}" y1="${y}" x2="${PX1}" y2="${y}" stroke="#ddd" stroke-width="0.5" stroke-dasharray="2,2"/><text x="${PX0 - 8}" y="${y + 3}" font-size="9" fill="#666" text-anchor="end">${v.toFixed(1)}</text>`;
    });
    els += `<text x="20" y="51" font-size="9" fill="#666" text-anchor="middle" transform="rotate(-90, 20, 51)">Relative intensity</text>`;
    // Draw order is independent of `entries`' own order (which drives the
    // top-left label stack and the data-entry index used elsewhere for
    // hover/tooltip lookups) -- entries flagged e.emphasize draw LAST, on
    // top of everything else. Without this, a proteoform whose peaks are
    // specifically curated to overlap another's (e.g. the section-2 target
    // vs. its m/z-collision confounders -- by construction, many of the
    // target's own peaks sit at nearly the confounders' exact m/z) gets
    // silently painted over by whichever entry happens to be drawn last in
    // array order, however visually important that entry actually is.
    const drawOrder = entries.map((_, i) => i).sort((a, b) => (entries[a].emphasize ? 1 : 0) - (entries[b].emphasize ? 1 : 0));
    drawOrder.forEach(i => {
      const e = entries[i];
      const strokeW = e.emphasize ? 2.4 : 1.5;
      const strokeOp = e.emphasize ? 1 : 0.85;
      (e.env || []).forEach(cs => {
        if (cs.resolved) {
          cs.points.forEach(p => {
            if (p.mz < loMz || p.mz > hiMz) return;
            const x = mzx(p.mz), h = p.rel * 70;
            // class + data-* attributes: hover target for the cross-
            // proteoform overlap tooltip (wireMS1PeakHover(), below) --
            // scoped to resolved discrete peaks only, not the unresolved
            // envelope curves, since "this specific peak" isn't a
            // well-defined thing to hover on a smooth curve the same way.
            els += `<line class="pt-ms1-peak" data-entry="${i}" data-mz="${p.mz}" data-rel="${p.rel}" data-z="${cs.z}" x1="${x.toFixed(1)}" y1="${baseY}" x2="${x.toFixed(1)}" y2="${(baseY - h).toFixed(1)}" stroke="${e.color}" stroke-width="${strokeW}" stroke-opacity="${strokeOp}"/>`;
          });
        } else {
          // Keep one point just outside each edge (if available) so the
          // curve drawn across a zoomed-in window still meets the chart
          // edges correctly instead of visibly starting/ending mid-air.
          const all = cs.points;
          let lo = all.findIndex(p => p.mz >= loMz);
          let hi = all.length - 1 - [...all].reverse().findIndex(p => p.mz <= hiMz);
          if (lo === -1 || hi < 0 || lo > hi) return;
          const pts = all.slice(Math.max(0, lo - 1), Math.min(all.length, hi + 2));
          if (pts.length < 2) return;
          let path = `M ${mzx(pts[0].mz).toFixed(1)} ${baseY}`;
          pts.forEach(p => { path += ` L ${mzx(p.mz).toFixed(1)} ${(baseY - p.rel * 70).toFixed(1)}`; });
          path += ` L ${mzx(pts[pts.length - 1].mz).toFixed(1)} ${baseY} Z`;
          els += `<path d="${path}" fill="${e.color}" fill-opacity="${e.emphasize ? 0.5 : 0.35}" stroke="${e.color}" stroke-width="${strokeW}" stroke-opacity="${strokeOp}"/>`;
        }
      });
    });
    const decimals = hiMz - loMz < 5 ? 3 : hiMz - loMz < 50 ? 1 : 0;
    const ticks = [loMz, (loMz + hiMz) / 2, hiMz];
    const anchors = ["start", "middle", "end"];
    ticks.forEach((t, i) => { els += `<text x="${mzx(t).toFixed(1)}" y="98" font-size="9" fill="#666" text-anchor="${anchors[i]}">${t.toFixed(decimals)} m/z</text>`; });
    // Entry labels (id + mass) sit BELOW the x-axis/m-z ticks, not overlaid
    // on the plot area itself (previously a top-left stack at font-size 10,
    // which both crowded the actual curves and, for a confounder search
    // with a couple dozen entries, ran off the bottom of the chart). Small
    // font, wrapped across a few columns so a large confounder list stays a
    // reasonable height instead of one tall single column.
    const labelFontSize = 8, labelLineH = 10, labelCols = entries.length > 8 ? 3 : entries.length > 3 ? 2 : 1;
    const labelColW = 640 / labelCols;
    const labelStartY = 112;
    const labelRows = Math.ceil(entries.length / labelCols);
    entries.forEach((e, i) => {
      const col = i % labelCols, row = Math.floor(i / labelCols);
      const x = 4 + col * labelColW, y = labelStartY + row * labelLineH;
      els += `<text x="${x}" y="${y}" font-size="${labelFontSize}" fill="${e.color}" font-weight="600">${e.label}: ${e.mass.toFixed(1)} Da</text>`;
    });
    const totalHeight = labelStartY - labelLineH + labelRows * labelLineH + 6;
    svgEl.setAttribute("viewBox", `0 0 640 ${totalHeight}`);
    svgEl.innerHTML = els;
  }

  // Real m/z values with decimals (finer once zoomed in enough that whole
  // numbers alone would be meaningless) -- e.g. "838.0-839.5 m/z (12%)" --
  // rather than the ladder chart's "residue N-M of total" wording, which
  // doesn't mean anything in the m/z domain.
  function ms1ZoomRangeText(state) {
    const lo = state.mzMin + state.vs, hi = state.mzMin + state.ve;
    const w = hi - lo;
    const decimals = w < 5 ? 3 : w < 50 ? 1 : 0;
    const pct = Math.round((w / state.axisLength) * 100);
    return `${lo.toFixed(decimals)}-${hi.toFixed(decimals)} m/z (${pct}%)`;
  }

  // Orbitrap resolving power / m/z-domain FWHM -- exact port of
  // resolving_power()/fwhm_mz() in R/resolving_power.R, so the hover
  // tooltip's "how many peak-widths apart" readout uses the SAME instrument
  // model the rest of the app already scores resolvability with, computed
  // client-side (no server round-trip per hover).
  function fwhmMzJs(mz, rRef, mzRef) {
    const rp = rRef * Math.sqrt(mzRef / mz);
    return mz / rp;
  }

  // Below this relative intensity, a point is numerical tail/noise, not a
  // real detectable peak -- comparing a hovered peak against a candidate
  // this faint (however close in m/z) produced a misleading "same detected
  // peak" verdict driven purely by matching against negligible signal, not
  // a real peak collision. Confirmed directly: at IL32 z=20, the two
  // proteoforms' actual APEX peaks (rel=1.0 each) are 0.70 m/z apart --
  // 24x the resolving-power FWHM, clearly resolved -- but an unfiltered
  // nearest-point search instead matched bare's apex against a mod-form
  // point at rel=0.0007 sitting 0.0008 m/z away, reporting "same peak"
  // for a comparison that was never physically meaningful.
  const MIN_PEAK_REL_FOR_COMPARISON = 0.01;

  // For a hovered peak (entries[entryIdx], at mz), the closest MEANINGFUL
  // peak (rel >= MIN_PEAK_REL_FOR_COMPARISON) belonging to a DIFFERENT
  // proteoform shown in the same chart -- scoped to other RESOLVED peaks
  // only (an unresolved charge state's envelope curve has no single
  // well-defined "peak" to compare a distance against the same way).
  // Searches every OTHER entry's points -- resolved discrete isotope peaks
  // AND unresolved envelope-curve samples alike, no distinction -- for the
  // nearest one to `mz`. Deliberately NOT scoped to resolved peaks only:
  // "is this proteoform distinguishable from that one here" is a pairwise
  // mass-vs-resolving-power question, and it has the same well-defined
  // answer regardless of whether either proteoform's OWN isotope fine
  // structure happens to be individually resolved into sticks or blurred
  // into one curve -- an 80 kDa protein's charge state not resolving its
  // own isotope comb doesn't change how far its apex sits from a checked
  // neighbor's apex. Treating curve sample points as ordinary points here
  // is exactly what lets ms1PeakTooltipHtml()'s resolving-power-ratio
  // wording ("clearly resolved" / "partially separated" / "effectively the
  // same detected peak") apply identically to both cases, instead of the
  // unresolved case getting a separate, differently-worded "overlap"
  // verdict with no quantitative ratio.
  function findNearestOtherPeak(entries, entryIdx, mz) {
    let best = null, bestDist = Infinity;
    entries.forEach((e, j) => {
      if (j === entryIdx) return;
      (e.env || []).forEach(cs => {
        cs.points.forEach(p => {
          if (p.rel < MIN_PEAK_REL_FOR_COMPARISON) return;
          const d = Math.abs(p.mz - mz);
          if (d < bestDist) { bestDist = d; best = { label: e.label, mass: e.mass, z: cs.z, mz: p.mz, rel: p.rel, deltaMz: d }; }
        });
      });
    });
    return best;
  }

  // Hover tooltip content for one MS1 peak: its own m/z/charge/intensity,
  // plus -- the point of this feature -- how it relates to the nearest
  // MEANINGFUL peak from any other proteoform shown in the same chart,
  // expressed in units of the instrument's own resolving-power FWHM at
  // that m/z (not a raw Da/m/z number alone, which doesn't by itself say
  // whether two peaks are actually distinguishable): under 1x means the
  // instrument would show these as literally the same detected peak; a few
  // x means partially separated; many x means cleanly resolved.
  function ms1PeakTooltipHtml(entries, entryIdx, mz, rel, z, rRef, mzRef) {
    const e = entries[entryIdx];
    let html = `<strong>${e.label}</strong><br/>z=${z}, m/z=${mz.toFixed(3)}, rel=${rel.toFixed(2)}`;
    const nearest = findNearestOtherPeak(entries, entryIdx, mz);
    if (nearest) {
      const fwhm = fwhmMzJs(mz, rRef, mzRef);
      const ratio = nearest.deltaMz / fwhm;
      let verdict;
      if (ratio < 1) verdict = "within the instrument's resolving width -- effectively the SAME detected peak";
      else if (ratio < 3) verdict = `partially separated (${ratio.toFixed(1)}&times; the resolving width apart)`;
      else verdict = `clearly resolved (${ratio.toFixed(1)}&times; the resolving width apart)`;
      html += `<br/><br/>Nearest other MEANINGFUL peak (&ge;${(MIN_PEAK_REL_FOR_COMPARISON * 100).toFixed(0)}% rel):<br/>` +
              `<strong>${nearest.label}</strong> (z=${nearest.z}, rel=${nearest.rel.toFixed(2)})` +
              `<br/>&Delta;m/z = ${nearest.deltaMz.toFixed(3)}<br/>${verdict}`;
    } else {
      html += `<br/><br/>No other proteoform has a peak &ge;${(MIN_PEAK_REL_FOR_COMPARISON * 100).toFixed(0)}% relative intensity nearby.`;
    }
    return html;
  }

  // Wires the MS1 chart's own pan/zoom, independent of the ladder chart's
  // (they're two different SVGs with two different domains -- m/z vs.
  // residue position -- so they need their own state object each, even
  // though both reuse the same generic attachPanZoom/zoomAt machinery), and
  // the peak-overlap hover tooltip (rRef/mzRef needed for the latter --
  // see ms1PeakTooltipHtml()). Returns the redraw function so the caller
  // can invoke it once for the initial paint.
  //
  // Hover target-finding searches the underlying DATA (nearest point to the
  // cursor's m/z, within a small pixel tolerance) rather than relying on
  // exactly which SVG element the cursor lands on. Deliberate: when two
  // proteoforms' combs are offset by close to an integer multiple of their
  // own isotope spacing (confirmed directly for IL32 bare vs +Methyl: 0.70
  // m/z apart, grid spacing 0.05 m/z, ratio = 14.02 -- almost exactly 14),
  // huge swaths of their two point sets sit within a fraction of a pixel of
  // each other, and whichever proteoform was drawn LAST paints on top and
  // silently swallows all mouse hits there -- reproduced directly: 30 of
  // bare's 44 points at z=20 had a mod-form point within 0.01 m/z. Hit-
  // testing the data instead of the DOM makes this exactly as easy to hover
  // regardless of paint order.
  // signal (an AbortSignal, optional): see attachPanZoom()'s doc comment --
  // pass one whenever the SAME ms1El/zoomEl pair may be re-wired more than
  // once (e.g. the MS2 fragment-isotope chart, re-rendered per bond click),
  // so the previous wiring's listeners are cleaned up first. The MS1/
  // confounder charts wire once per render and can omit it.
  function wireMS1Zoom(idPrefix, ms1El, zoomEl, entries, rRef, mzRef, signal) {
    const domain = computeMS1Domain(entries);
    if (!domain) { ms1El.innerHTML = ""; return () => {}; }
    const state = { vs: 1, ve: domain.axisLength, axisLength: domain.axisLength, mzMin: domain.mzMin };
    function redraw() { drawMS1(ms1El, entries, state); }
    if (zoomEl) {
      zoomEl.innerHTML = zoomControlsHtml(idPrefix);
      wireZoomButtons(idPrefix, state, redraw, domain.minZoomWidth, ms1ZoomRangeText);
    }
    attachPanZoom(ms1El, state, redraw, {
      minZoomWidth: domain.minZoomWidth,
      onZoomChange: () => updateZoomRangeLabel(idPrefix, state, ms1ZoomRangeText),
      signal
    });
    ms1El.addEventListener("mousemove", ev => {
      const rect = ms1El.getBoundingClientRect();
      const vb = ms1El.viewBox.baseVal;
      const px = (ev.clientX - rect.left) * vb.width / rect.width;
      const py = (ev.clientY - rect.top) * vb.height / rect.height;
      if (py < 8 || py > 92) { hideTooltip(); return; } // outside the plot area (axis-label rows)
      const mzCursor = state.mzMin + posFromPx(px, state.vs, state.ve);
      const pxPerUnit = (PX1 - PX0) / (state.ve - state.vs);
      const maxMzDist = 6 / pxPerUnit; // ~6px hover tolerance
      // Searches every point (resolved discrete isotope peaks AND
      // unresolved envelope-curve samples alike -- no `cs.resolved` gate)
      // for the nearest one to the cursor. Unified deliberately: whether a
      // GIVEN charge state's own isotope fine structure happens to resolve
      // into sticks or blur into one curve doesn't change the pairwise
      // "is this proteoform distinguishable from that one here" question,
      // so both cases get the exact same resolving-power-ratio verdict via
      // ms1PeakTooltipHtml() below, instead of a separately-worded
      // "overlap yes/no" verdict for unresolved charge states with no
      // quantitative ratio -- confirmed as a real point of user confusion
      // when the wording differed between a small, mostly-resolved
      // proteoform (IL32) and a large, entirely-unresolved one (CD44) in
      // the SAME comparison.
      let best = null, bestDist = maxMzDist;
      entries.forEach((e, i) => {
        (e.env || []).forEach(cs => {
          cs.points.forEach(p => {
            // Same floor as findNearestOtherPeak(), applied here too --
            // without it, the PRIMARY hover target itself could resolve to
            // a negligible-intensity point instead of the real peak the
            // user is visually pointing at. Confirmed directly: two
            // proteoforms offset by close to an integer multiple of their
            // own isotope spacing (IL32 bare vs +Methyl: ratio 14.02) have
            // huge stretches where a real, tall peak in one and a rel~0
            // "ghost" point in the other sit within a thousandth of an m/z
            // unit of each other -- sub-pixel mouse-precision noise alone
            // was enough to flip which one won, silently swapping which
            // PROTEOFORM the tooltip claimed was even being hovered.
            if (p.rel < MIN_PEAK_REL_FOR_COMPARISON) return;
            const d = Math.abs(p.mz - mzCursor);
            if (d < bestDist) { bestDist = d; best = { entryIdx: i, mz: p.mz, rel: p.rel, z: cs.z }; }
          });
        });
      });
      if (!best) { hideTooltip(); return; }
      showTooltip(ev.clientX, ev.clientY, ms1PeakTooltipHtml(entries, best.entryIdx, best.mz, best.rel, best.z, rRef, mzRef));
    }, { signal });
    ms1El.addEventListener("mouseleave", hideTooltip, { signal });
    return redraw;
  }

  function exonTier(block, others) {
    if (others.length === 0) return "neutral";
    let mc = 0;
    others.forEach(o => { if (o.exon_blocks.some(b => b.start === block.start && b.end === block.end)) mc++; });
    if (mc === others.length) return "common";
    if (mc === 0) return "unique";
    return "partial";
  }

  // Compact unique/partial/common badges (dot + count, no text label --
  // the shared common/partial/unique legend, s1-legend/s2-legend, already
  // explains the colors once) drawn INLINE inside the ladder SVG, to the
  // right of a row's own title (transcript name + "N aa, M Da"), instead of
  // a separate HTML column that had to be height-synced against the SVG's
  // own responsive size. Right-aligned in 3 fixed-width slots ending near
  // PX1 -- clear of the title text (which never reaches this far right) and
  // clear of the exon/tick rows below (a different y range entirely); the
  // only element that could ever share this (x, y) region is a PTM marker
  // whose residue happens to map into the chart's rightmost ~90px at the
  // CURRENT zoom, an occasional minor visual crowd rather than a
  // structural conflict.
  function rowBadgesSvg(counts, top) {
    const badgeY = top + 16;
    const slotW = 34, blockEnd = PX1, blockStart = blockEnd - 3 * slotW;
    const tiers = [["unique", counts.unique], ["partial", counts.partial], ["common", counts.common]];
    let s = "";
    tiers.forEach(([tier, n], i) => {
      const sx = blockStart + i * slotW;
      s += `<circle cx="${sx}" cy="${badgeY - 3}" r="3.2" fill="${TIER_COLOR[tier]}"/>`;
      s += `<text x="${sx + 6}" y="${badgeY}" font-size="9" font-weight="600" fill="#333">${n}</text>`;
    });
    return s;
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

      // MS2 unique/partial/common tally for THIS row, recomputed here (not
      // read from the server-computed, filter-independent ms2_stats) so it
      // reflects only the bonds the CURRENT fragment filter is showing --
      // every redraw() call (pan/zoom, and critically every filter-button
      // click) recomputes it fresh. Both the b- and y-ion at a passing bond
      // count, matching tally_ms2_tiers()'s server-side definition, just
      // scoped to the visible subset instead of the full ladder.
      const rowCounts = { unique: 0, partial: 0, common: 0 };
      for (let p = 1; p <= pf.b_mass.length; p++) {
        if (!passesFilter(pf.propensity[p - 1], state.filter, state.scoringMode)) continue;
        if (rowCounts[pf.tier_b[p - 1]] !== undefined) rowCounts[pf.tier_b[p - 1]]++;
        if (rowCounts[pf.tier_y[p - 1]] !== undefined) rowCounts[pf.tier_y[p - 1]]++;
      }

      els += `<text x="4" y="${top + 9}" font-size="10" font-weight="600" fill="${pf.color}">${pf.label}</text>`;
      els += `<text x="4" y="${top + 20}" font-size="9" fill="#666">${pf.len} aa, ${pf.mass.toFixed(1)} Da</text>`;
      els += rowBadgesSvg(rowCounts, top);
      els += `<line x1="${PX0}" y1="${top + 24}" x2="${PX1}" y2="${top + 24}" stroke="#ddd" stroke-width="1"/>`;

      pf.ptms.forEach(ptm => {
        const x = xScale(ptm.axis_pos, state.vs, state.ve);
        if (x < PX0 - 2 || x > PX1 + 2) return;
        els += `<line x1="${x.toFixed(1)}" y1="${top + 24}" x2="${x.toFixed(1)}" y2="${top + 13}" stroke="#e0433d" stroke-width="1.5"/>`;
        els += `<circle cx="${x.toFixed(1)}" cy="${top + 11}" r="3" fill="#e0433d"/>`;
        els += `<text x="${x.toFixed(1)}" y="${top + 9}" font-size="7" text-anchor="middle" fill="#555">${ptm.name}</text>`;
      });

      els += `<text x="4" y="${top + 35}" font-size="9" fill="#999">exon</text>`;
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
        if (!passesFilter(score, state.filter, state.scoringMode)) continue;
        visibleCount++;
        const x = xScale(axisP, state.vs, state.ve);
        if (x < PX0 - 2 || x > PX1 + 2) continue;
        const [normLo, normHi] = tierThresholds(state.scoringMode).opNorm;
        const op = state.filter === 0 ? (0.22 + 0.78 * Math.max(0, Math.min(1, (score - normLo) / (normHi - normLo)))) : 0.95;
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

  // Keyed by the ladder <svg> element itself (WeakMap, no manual cleanup
  // needed) rather than an idPrefix string -- same reasoning as
  // wireMs1ChartOnce() above: ladderEl is STATIC ui.R markup that outlives
  // any single "Run analysis" click, so re-wiring it on every click without
  // aborting the previous listeners would stack duplicates.
  const ladderInteractionAborts = new WeakMap();
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
    const prev = ladderInteractionAborts.get(svgEl);
    if (prev) prev.abort();
    const ctrl = new AbortController();
    ladderInteractionAborts.set(svgEl, ctrl);
    attachPanZoom(svgEl, state, redraw, { onHover, onLeave, onZoomChange: () => hideTooltip(), signal: ctrl.signal });
  }

  // Cross-proteoform isotope-peak comparison for ONE specific fragment ion
  // (b or y at one backbone position), shown on click of its ladder tick.
  // Unlike the MS1 chart -- always the SAME molecule's charge states,
  // comparable across the whole checked set on one shared m/z axis -- a
  // ladder has hundreds of DIFFERENT fragment ions at wildly different
  // masses; plotting all of them at once would be meaningless, so this is
  // deliberately click-to-select rather than always-on. Matches "the same
  // fragment ion" across proteoforms by shared axis position + ion type
  // (build_section1_payload()'s axis_pos -- the same coordinate the
  // ladder's own tick x-position already uses), not by local bond number,
  // so it stays correct even when the compared proteoforms have different
  // lengths (different N-terminal exon structure, or a middle-down
  // peptide's own local numbering).
  //
  // WHICH bonds qualify (propensity above this mode's isotope gate) is decided entirely from data
  // already in `proteoforms` (propensity/tier_b/tier_y/axis_pos -- cheap
  // numbers, needed for the ladder chart anyway). The actual isotope
  // PATTERN for a qualifying fragment is NOT sent up front: a ~750-residue
  // proteoform can have 100+ qualifying bonds, and eagerly computing +
  // embedding all of their isotope patterns was confirmed directly to cost
  // ~2MB of JSON and ~16s of IsoSpecPy calls PER proteoform -- exactly what
  // made "Run analysis" crawl (and briefly hang the page) on a 3-way,
  // ~80kDa comparison. Instead, a click sends a small request
  // (Shiny.setInputValue('frag_ms1_request', ...)) for just the 1-3
  // fragments actually being looked at, and the server computes only those
  // on demand (server.R's input$frag_ms1_request observer).
  const fragSections = {}; // sectionPrefix -> {proteoforms, rRef, mzRef, activeRequestId, listenerAbort, zoomAbort}
  let fragRequestSeq = 0;
  let fragResponseHandlerRegistered = false;

  // Per-idPrefix AbortController for the ALWAYS-ON MS1 charge-envelope
  // chart's own wireMS1Zoom() wiring (s1-ms1/s2-ms1), so repeat "Run
  // analysis" clicks -- renderSection1()/renderSection2() run once per
  // click, but ms1El itself is STATIC ui.R markup that persists across
  // every click, unlike the fragment chart's per-click re-wire -- abort the
  // previous wiring's listeners before attaching new ones. Without this,
  // Shiny re-executing viz_script's <script> tag (which it does on every
  // DOM rebind of that output, not strictly once per "Run analysis") stacks
  // duplicate mousemove listeners closing over progressively stale
  // entries/state objects; confirmed directly as the cause of the MS1
  // chart's hover tooltip going completely silent on an UNRESOLVED curve
  // (~81 kDa, 3-proteoform comparison) while the SAME hover logic, wired
  // with this same cleanup pattern, worked correctly on the MS2 fragment
  // chart right next to it.
  const ms1ChartAborts = {};
  function wireMs1ChartOnce(idPrefix, ms1El, zoomEl, entries, rRef, mzRef) {
    if (ms1ChartAborts[idPrefix]) ms1ChartAborts[idPrefix].abort();
    const ctrl = new AbortController();
    ms1ChartAborts[idPrefix] = ctrl;
    wireMS1Zoom(idPrefix, ms1El, zoomEl, entries, rRef, mzRef, ctrl.signal)();
  }

  function registerFragResponseHandler() {
    if (fragResponseHandlerRegistered || !window.Shiny) return;
    fragResponseHandlerRegistered = true;
    Shiny.addCustomMessageHandler("pt_frag_ms1_response", data => {
      const sec = fragSections[data.section];
      if (!sec || sec.activeRequestId !== data.requestId) return; // stale (superseded by a later click)
      const fragEl = document.getElementById(data.section + "-frag-ms1");
      const fragZoomEl = document.getElementById(data.section + "-frag-ms1-zoom");
      const labelEl = document.getElementById(data.section + "-frag-label");
      if (!fragEl) return;
      const byId = Object.fromEntries(sec.proteoforms.map(pf => [pf.id, pf]));
      const entries = data.entries.map(e => {
        const pf = byId[e.id] || {};
        const score = (sec.lastScores || {})[e.id + ":" + e.p + ":" + e.ion];
        const scoreText = score == null ? "" : ` (score ${score})`;
        return { label: `${pf.label || e.id} ${e.ion}${e.p}${scoreText}`, color: pf.color || "#333", mass: e.mass, env: e.env };
      });
      if (labelEl) {
        labelEl.textContent = sec.lastLabelPrefix + (entries.length === 0
          ? " -- no isotope data returned."
          : ` ${entries.length} of ${sec.lastQualifiedTotal} checked proteoform(s) shown below.`);
      }
      if (sec.zoomAbort) sec.zoomAbort.abort();
      sec.zoomAbort = new AbortController();
      wireMS1Zoom(data.section + "-frag-ms1", fragEl, fragZoomEl, entries, sec.rRef, sec.mzRef, sec.zoomAbort.signal)();
    });
  }

  function wireFragmentClick(sectionPrefix, ladderEl, proteoforms, rRef, mzRef, state) {
    const fragEl = document.getElementById(sectionPrefix + "-frag-ms1");
    const fragZoomEl = document.getElementById(sectionPrefix + "-frag-ms1-zoom");
    const labelEl = document.getElementById(sectionPrefix + "-frag-label");
    if (!fragEl) return;
    registerFragResponseHandler();

    const sec = fragSections[sectionPrefix] || {};
    sec.proteoforms = proteoforms; sec.rRef = rRef; sec.mzRef = mzRef;
    fragSections[sectionPrefix] = sec;

    // Re-wired on every "Run analysis" re-render (unlike the MS1/ladder
    // charts' one-time-per-render wiring, this same ladderEl persists
    // across MULTIPLE renders) -- abort the previous click listener first,
    // same AbortController pattern as attachPanZoom()/wireMS1Zoom(), or
    // repeat clicks would fire once per past render, each with stale data.
    if (sec.listenerAbort) sec.listenerAbort.abort();
    sec.listenerAbort = new AbortController();

    ladderEl.addEventListener("click", ev => {
      const rect = ladderEl.getBoundingClientRect();
      const vb = ladderEl.viewBox.baseVal;
      const px = (ev.clientX - rect.left) * vb.width / rect.width;
      const py = (ev.clientY - rect.top) * vb.height / rect.height;
      const hit = findNearestTick(state, proteoforms, px, py);
      if (!hit) return;
      const clickedPf = proteoforms[hit.ri];
      const axisPos = clickedPf.axis_pos[hit.p - 1];

      // Independent of whether the CLICKED proteoform's own bond here
      // qualifies (propensity above this mode's Elevated threshold) --
      // another checked proteoform can have a qualifying fragment at this
      // same aligned axis position even when the one actually clicked
      // doesn't. Only cheap per-bond arrays are scanned here -- no isotope
      // data exists client-side yet. Threshold is scoring-mode-aware since
      // "glm" and "rf" scores are on entirely different numeric scales.
      const isotopeGate = state.scoringMode === "rf" ? RF_ISOTOPE_GATE_JS : GLM_ISOTOPE_GATE_JS;
      const qualified = [];
      proteoforms.forEach(pf => {
        const localP = pf.axis_pos.indexOf(axisPos) + 1; // 1-based to match p
        if (localP < 1) return;
        const score = pf.propensity[localP - 1];
        if (!(score > isotopeGate)) return;
        qualified.push({ id: pf.id, p: localP, score });
      });

      if (qualified.length === 0) {
        if (sec.zoomAbort) sec.zoomAbort.abort();
        fragEl.innerHTML = "";
        if (fragZoomEl) fragZoomEl.innerHTML = "";
        if (labelEl) {
          labelEl.textContent = `${hit.ion}${hit.p} (${clickedPf.label}): fragmentation-propensity score at or below ${isotopeGate} for every checked proteoform at this bond -- no isotope-peak detail computed (bond unlikely to actually fragment).`;
        }
        return;
      }

      sec.lastLabelPrefix = `${hit.ion}-ion at axis position ${axisPos}:`;
      sec.lastQualifiedTotal = qualified.length;
      sec.lastScores = {};
      qualified.forEach(q => { sec.lastScores[q.id + ":" + q.p + ":" + hit.ion] = q.score; });
      if (labelEl) labelEl.textContent = `${sec.lastLabelPrefix} computing isotope peaks for ${qualified.length} of ${proteoforms.length} checked proteoform(s)...`;

      const requestId = ++fragRequestSeq;
      sec.activeRequestId = requestId;
      Shiny.setInputValue("frag_ms1_request", {
        section: sectionPrefix, requestId,
        requests: qualified.map(q => ({ id: q.id, p: q.p, ion: hit.ion }))
      }, { priority: "event" });
    }, { signal: sec.listenerAbort.signal });
  }

  function renderSection1(payload) {
    const ms1El = document.getElementById("s1-ms1");
    const ms1ZoomEl = document.getElementById("s1-ms1-zoom");
    const ladderEl = document.getElementById("s1-ladder");
    const legendEl = document.getElementById("s1-legend");
    const filterEl = document.getElementById("s1-filters");
    const zoomEl = document.getElementById("s1-zoom");
    const infoEl = document.getElementById("s1-info");
    if (!ms1El || !ladderEl) return;
    expandCollapsible("s1-collapse-body");

    const proteoforms = payload.proteoforms;
    const state = { vs: 1, ve: payload.axis_length, filter: 0, axisLength: payload.axis_length, scoringMode: payload.scoring_mode || "glm" };

    const ms1Entries = proteoforms.map(p => ({ label: p.label, color: p.color, mass: p.mass, env: p.env || [] }));
    wireMs1ChartOnce("s1-ms1", ms1El, ms1ZoomEl, ms1Entries, payload.r_ref, payload.mz_ref);
    legendEl.innerHTML = legendHtml(proteoforms.length < 2);
    filterEl.innerHTML = filterButtonsHtml("s1", state.filter, state.scoringMode);
    if (zoomEl) zoomEl.innerHTML = zoomControlsHtml("s1");
    renderStatsStrip(document.getElementById("s1-stats-strip"), payload.ms1_stats, "MS1 charge-envelope peaks");

    function redraw() {
      const counts = renderLadderRows(ladderEl, proteoforms, state.axisLength, state);
      const cEl = document.getElementById("s1-fcount");
      if (cEl) cEl.textContent = counts.visibleCount + " of " + counts.totalCount + " bonds shown";
      updateZoomRangeLabel("s1", state);
    }
    redraw();
    wireInteraction(ladderEl, state, proteoforms, redraw, infoEl);
    wireFragmentClick("s1", ladderEl, proteoforms, payload.r_ref, payload.mz_ref, state);
    if (zoomEl) wireZoomButtons("s1", state, redraw);
    function rewireFilterButtons() {
      tierThresholds(state.scoringMode).labels.forEach((_, i) => {
        const btn = document.getElementById("s1-f" + i);
        if (btn) btn.onclick = () => { state.filter = i; filterEl.innerHTML = filterButtonsHtml("s1", i, state.scoringMode); redraw(); rewireFilterButtons(); };
      });
    }
    rewireFilterButtons();
  }

  function renderSection2(payload) {
    const ms1El = document.getElementById("s2-ms1");
    const ms1ZoomEl = document.getElementById("s2-ms1-zoom");
    const ladderEl = document.getElementById("s2-ladder");
    const legendEl = document.getElementById("s2-legend");
    const filterEl = document.getElementById("s2-filters");
    const zoomEl = document.getElementById("s2-zoom");
    const infoEl = document.getElementById("s2-info");
    if (!ms1El || !ladderEl) return;
    expandCollapsible("s2-collapse-body");
    // The MS1/MS2 sub-panels start collapsed too (ui.R) -- only pop open
    // once "Compare selected confounders" (or a later scoring-mode
    // rescore, which also calls this) actually produces something to show,
    // same reasoning as the outer s2-collapse-body above.
    expandCollapsible("s2-ms1-collapse-body");
    expandCollapsible("s2-ms2-collapse-body");

    const t = payload.target;
    const targetColor = "#1f8a70";
    const entries = [{ label: "TARGET: " + t.id, color: targetColor, mass: t.mass, env: t.env || [], emphasize: true }];
    payload.confounders.forEach((c, i) => entries.push({ label: c.id, color: CONFOUNDER_PALETTE[i % CONFOUNDER_PALETTE.length], mass: c.mass, env: c.env || [] }));
    wireMs1ChartOnce("s2-ms1", ms1El, ms1ZoomEl, entries, payload.r_ref, payload.mz_ref);

    // Target + every confounder, each on its OWN residue axis (0..len) --
    // unlike section 1's checked proteoforms (usually one gene, so they can
    // share a real exon axis), a confounder is by definition an unrelated
    // protein from a different gene, so there's no shared coordinate system
    // to align them on. Same fallback section 1 itself uses when no exon
    // table is available (build_section1_payload(), R/viz_json.R). The MS2
    // ladder used to only draw the target's own row here even though the
    // MS1 chart above already overlays every confounder -- this makes both
    // charts show the same set of proteins.
    const proteoforms = [{
      id: t.id, label: "TARGET: " + t.id, color: targetColor, mass: t.mass, len: t.len, sequence: t.sequence,
      ptms: t.ptms.map(p => Object.assign({}, p, { axis_pos: p.pos })),
      exon_blocks: t.exon_blocks.map(b => ({ start: b.start, end: b.end })),
      b_mass: t.b_mass, y_mass: t.y_mass, propensity: t.propensity,
      axis_pos: t.b_mass.map((_, i) => i + 1),
      tier_b: t.tier_b, tier_y: t.tier_y, ms2_stats: t.ms2_stats
    }];
    payload.confounders.forEach((c, i) => {
      if (!c.b_mass || !c.b_mass.length) return; // no ladder computed for this candidate (shouldn't normally happen)
      proteoforms.push({
        id: c.id, label: c.id, color: CONFOUNDER_PALETTE[i % CONFOUNDER_PALETTE.length], mass: c.mass,
        len: c.len, sequence: c.sequence, ptms: [], exon_blocks: [],
        b_mass: c.b_mass, y_mass: c.y_mass, propensity: c.propensity,
        axis_pos: c.b_mass.map((_, j) => j + 1),
        tier_b: c.tier_b, tier_y: c.tier_y, ms2_stats: c.ms2_stats
      });
    });
    const axisLength = Math.max(...proteoforms.map(p => p.len || 1));
    const state = { vs: 1, ve: axisLength, filter: 0, axisLength: axisLength, scoringMode: payload.scoring_mode || "glm" };
    legendEl.innerHTML = legendHtml(payload.confounders.length === 0) +
      `<div style="font-size:12px;color:#666;">Window: +/-${payload.window_da} Da at best charge state ${payload.best_charge_state}. ${payload.confounders.length} real confounder(s) found.</div>`;
    renderStatsStrip(document.getElementById("s2-stats-strip"), payload.ms1_stats, "Target's MS1 peaks");
    filterEl.innerHTML = filterButtonsHtml("s2", state.filter, state.scoringMode);
    if (zoomEl) zoomEl.innerHTML = zoomControlsHtml("s2");

    function redraw() {
      const counts = renderLadderRows(ladderEl, proteoforms, state.axisLength, state);
      const cEl = document.getElementById("s2-fcount");
      if (cEl) cEl.textContent = counts.visibleCount + " of " + counts.totalCount + " bonds shown";
      updateZoomRangeLabel("s2", state);
    }
    redraw();
    wireInteraction(ladderEl, state, proteoforms, redraw, infoEl);
    wireFragmentClick("s2", ladderEl, proteoforms, payload.r_ref, payload.mz_ref, state);
    if (zoomEl) wireZoomButtons("s2", state, redraw);
    function rewireFilterButtons() {
      tierThresholds(state.scoringMode).labels.forEach((_, i) => {
        const btn = document.getElementById("s2-f" + i);
        if (btn) btn.onclick = () => { state.filter = i; filterEl.innerHTML = filterButtonsHtml("s2", i, state.scoringMode); redraw(); rewireFilterButtons(); };
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
    // Tracks which highlights actually landed on a real block below (drawn
    // as a box on top of that block) -- an event arm can have ZERO
    // matching transcripts (e.g. one of MXE's two mutually exclusive
    // exons), in which case NO track ever has a block at that exact
    // position, so its box could never be drawn there. Those get a
    // fallback marker afterwards instead of silently going unhighlighted.
    const highlightMatched = (payload.highlights || []).map(() => false);
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
        // Box any block that IS one of the event's own differential
        // exon(s) (e.g. Option 3's rMATS-flagged exon(s)) -- drawn as a
        // dashed outline on top of the normal tier fill above, exact-match
        // on axis coordinates since both are derived from the same
        // genomic->axis mapping. A track that doesn't have this exon
        // simply has no block here at all, so no box is drawn for it --
        // that absence is itself the "exon skipped" signal.
        (payload.highlights || []).forEach((hl, hi) => {
          if (Math.abs(b.start - hl.start) > 0.5 || Math.abs(b.end - hl.end) > 0.5) return;
          const hcolor = HIGHLIGHT_COLORS[hi % HIGHLIGHT_COLORS.length];
          els += `<rect x="${x0.toFixed(1)}" y="${top + 10}" width="${Math.max(1, x1 - x0).toFixed(1)}" height="12" fill="none" stroke="${hcolor}" stroke-width="2" stroke-dasharray="3,2" rx="1.5"/>`;
          highlightMatched[hi] = true;
        });
      });
    });

    // Fallback for any highlight that never landed on a real block above
    // (its arm had zero matching transcripts, so no track has an exon
    // there) -- a full-height dashed vertical band plus a small label, so
    // the event's differential exon(s) are always shown somewhere even
    // when no candidate transcript actually contains them.
    (payload.highlights || []).forEach((hl, hi) => {
      if (highlightMatched[hi]) return;
      const hcolor = HIGHLIGHT_COLORS[hi % HIGHLIGHT_COLORS.length];
      const x0 = Math.max(xs(hl.start), PX0), x1 = Math.min(xs(hl.end), PX1);
      if (x1 < PX0 || x0 > PX1) return;
      const cx0 = Math.max(x0, PX0), cx1 = Math.min(Math.max(x1, x0 + 1), PX1);
      els += `<rect x="${cx0.toFixed(1)}" y="${(topPad - 4).toFixed(1)}" width="${Math.max(1, cx1 - cx0).toFixed(1)}" height="${(height - topPad + 4).toFixed(1)}" fill="none" stroke="${hcolor}" stroke-width="2" stroke-dasharray="3,2"/>`;
      els += `<text x="${((cx0 + cx1) / 2).toFixed(1)}" y="${(topPad - 8).toFixed(1)}" font-size="8" fill="${hcolor}" text-anchor="middle">${hl.label || "differential exon"} (no matching transcript)</text>`;
    });

    els += `<text x="${PX0}" y="${height - 4}" font-size="9" fill="#888">${payload.seqname}, ${payload.strand} strand -- exon widths shown to scale relative to each other; introns shown as connecting lines, not to scale</text>`;

    svgEl.setAttribute("viewBox", `0 0 640 ${height}`);
    svgEl.innerHTML = els;
  }

  function exonAlignmentLegendHtml(highlights) {
    let h = '<div style="display:flex;flex-direction:column;gap:4px;padding:4px 0;font-size:12px;color:#555;">' +
      '<div style="display:flex;gap:16px;flex-wrap:wrap;">' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;border-radius:2px;display:inline-block;"></span>present in all compared</span>' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#8952e0;border-radius:2px;display:inline-block;"></span>present in some, not all</span>' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#2a78d6;border-radius:2px;display:inline-block;"></span>present in only one</span>' +
      "</div>" +
      '<div style="display:flex;gap:16px;flex-wrap:wrap;">' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;border-radius:2px;display:inline-block;"></span>coding (CDS)</span>' +
      '<span style="display:flex;align-items:center;gap:6px;"><span style="width:10px;height:10px;background:#eda100;opacity:0.3;border-radius:2px;display:inline-block;"></span>non-coding (UTR)</span>' +
      "</div>";
    if (highlights && highlights.length) {
      h += '<div style="display:flex;gap:16px;flex-wrap:wrap;">' +
        highlights.map((hl, hi) => {
          const color = HIGHLIGHT_COLORS[hi % HIGHLIGHT_COLORS.length];
          const label = hl.label || "differential exon";
          return `<span style="display:flex;align-items:center;gap:6px;"><span style="width:14px;height:10px;border:2px dashed ${color};border-radius:2px;display:inline-block;"></span>${label}</span>`;
        }).join("") +
        "</div>";
    }
    return h + "</div>";
  }

  // Formats a bare sequence into numbered blocks for the PTM-lookup
  // popover: e.g. residues 1-25 as
  //   1     6     11    16    21
  //   MDKFW WHAAW GLCLV PLSLA ABCDE
  // -- a ruler line (the residue number at the START of each 5-residue
  // group) directly above the sequence line, monospace-aligned so each
  // number's first digit sits above its group's first letter. groupSize=5
  // (not the more common 10) is deliberate -- finer-grained position
  // anchors are the whole point of this feature (quickly counting "which
  // group is residue 133 in", not just "which line").
  function formatSequenceWithRuler(seq, groupSize, groupsPerLine) {
    groupSize = groupSize || 5;
    groupsPerLine = groupsPerLine || 10;
    const lineWidth = groupSize * groupsPerLine;
    const lines = [];
    for (let start = 0; start < seq.length; start += lineWidth) {
      const lineSeq = seq.slice(start, start + lineWidth);
      const groups = [];
      let ruler = "";
      for (let g = 0; g < lineSeq.length; g += groupSize) {
        const groupSeq = lineSeq.slice(g, g + groupSize);
        groups.push(groupSeq);
        const pos = start + g + 1; // 1-based
        ruler += String(pos).padEnd(groupSeq.length + 1, " ");
      }
      lines.push(ruler.replace(/\s+$/, "") + "\n" + groups.join(" "));
    }
    return lines.join("\n\n");
  }

  // Sequence-lookup popover for the isoform catalog table: hover a
  // transcript id (.pt-seq-hover, data-seq carries the raw sequence -- see
  // server.R's isoform_catalog_ui row_for()) to see its full sequence with
  // residue-position anchors, for filling in the PTM spec box right next to
  // it (e.g. "133_Thr_Phospho") without having to count residues by eye.
  //
  // A SEPARATE element from #pt-hover-tooltip (used for the b/y ion quick
  // glance) rather than reusing it: that one is pointer-events:none and
  // 260px wide, both wrong here -- this needs to stay open while the mouse
  // moves INTO it (so its content can actually be read/selected/copied) and
  // needs real width for a monospace block. A short close-delay on
  // mouseleave (both the trigger's and the popover's own) is what makes
  // "move mouse from trigger into the popover" work without a flicker-close
  // in between.
  function ensureSeqPopover() {
    let el = document.getElementById("pt-seq-popover");
    if (!el) {
      el = document.createElement("div");
      el.id = "pt-seq-popover";
      document.body.appendChild(el);
    }
    return el;
  }

  function wireSeqHoverDelegation() {
    // Lazily created/wired on first real hover, not here -- this whole
    // function runs at SCRIPT LOAD time, while the <script> tag (in
    // tags$head()) is still being parsed, before document.body exists yet.
    // ensureSeqPopover()'s document.body.appendChild() would throw if
    // called this early; the mouseover/mouseout listeners below are safe to
    // attach to `document` immediately (it always exists), since they only
    // actually fire later, well after body exists.
    let popover = null;
    let wired = false;
    let closeTimer = null;
    function cancelClose() { if (closeTimer) { clearTimeout(closeTimer); closeTimer = null; } }
    function scheduleClose() { cancelClose(); closeTimer = setTimeout(() => { if (popover) popover.style.display = "none"; }, 250); }

    function showFor(trigger) {
      const seq = trigger.getAttribute("data-seq");
      if (!seq) return;
      if (!popover) {
        popover = ensureSeqPopover();
        if (!wired) {
          popover.addEventListener("mouseover", cancelClose);
          popover.addEventListener("mouseout", e => { if (!popover.contains(e.relatedTarget)) scheduleClose(); });
          wired = true;
        }
      }
      const id = trigger.textContent;
      popover.innerHTML = `<div class="pt-seq-popover-header">${id} (${seq.length} aa)</div><pre>${formatSequenceWithRuler(seq)}</pre>`;
      const r = trigger.getBoundingClientRect();
      popover.style.display = "block";
      // Measure after display:block (offsetWidth/Height are 0 while
      // display:none) so the viewport-edge clamp below uses real values.
      const pw = popover.offsetWidth, ph = popover.offsetHeight;
      let left = r.left, top = r.bottom + 4;
      if (left + pw > window.innerWidth - 8) left = Math.max(8, window.innerWidth - pw - 8);
      if (top + ph > window.innerHeight - 8) top = Math.max(8, r.top - ph - 4);
      popover.style.left = left + "px";
      popover.style.top = top + "px";
    }

    document.addEventListener("mouseover", e => {
      const trigger = e.target.closest && e.target.closest(".pt-seq-hover");
      if (trigger) { cancelClose(); showFor(trigger); }
    });
    document.addEventListener("mouseout", e => {
      const trigger = e.target.closest && e.target.closest(".pt-seq-hover");
      if (trigger && !trigger.contains(e.relatedTarget)) scheduleClose();
    });
  }
  wireSeqHoverDelegation();

  // Collapse/expand toggle for the isoform catalog table (ui.R's
  // #isoform-catalog-toggle / #isoform-catalog-collapse-body) -- a long
  // isoform list otherwise forces scrolling all the way past a list the
  // user is already done picking from just to reach the results below.
  // Click delegation on `document` (not the header element directly) for
  // the same reason as wireSeqHoverDelegation() above: this whole IIFE
  // runs while the <script> tag (in tags$head()) is still being parsed,
  // before document.body -- and hence the header itself -- exists yet.
  function updateIsoformCatalogSummary() {
    const body = document.getElementById("isoform-catalog-collapse-body");
    const summaryEl = document.getElementById("isoform-catalog-summary");
    if (!body || !summaryEl) return;
    const boxes = Array.from(body.querySelectorAll('.pt-isorow input[type="checkbox"]'));
    const checked = boxes.filter(b => b.checked).length;
    summaryEl.textContent = boxes.length ? `(${checked} of ${boxes.length} selected)` : "";
  }

  function wireIsoformCatalogCollapse() {
    document.addEventListener("click", e => {
      const header = e.target.closest && e.target.closest("#isoform-catalog-toggle");
      if (!header) return;
      const body = document.getElementById("isoform-catalog-collapse-body");
      const triangle = document.getElementById("isoform-catalog-triangle");
      if (!body || !triangle) return;
      const collapsed = body.style.display === "none";
      body.style.display = collapsed ? "" : "none";
      triangle.classList.toggle("pt-collapsed", !collapsed);
    });
    // Keep the "(N of M selected)" summary live so it's still trustworthy
    // while collapsed -- delegated the same way (the checkboxes themselves
    // are inside a renderUI() block that gets replaced wholesale on every
    // new gene load, so listeners bound directly to them would go stale).
    document.addEventListener("change", e => {
      if (e.target.matches && e.target.matches('.pt-isorow input[type="checkbox"]') &&
          e.target.closest("#isoform-catalog-collapse-body")) {
        updateIsoformCatalogSummary();
      }
    });
    // A genuinely NEW isoform_catalog_ui render (new gene loaded) resets
    // the panel back to expanded -- collapsing is a per-session convenience
    // for a list the user has already finished with, not something that
    // should carry over and hide a freshly loaded, unreviewed list.
    if (window.jQuery) {
      jQuery(document).on("shiny:value", function (event) {
        if (event.name !== "isoform_catalog_ui") return;
        const body = document.getElementById("isoform-catalog-collapse-body");
        const triangle = document.getElementById("isoform-catalog-triangle");
        if (body) body.style.display = "";
        if (triangle) triangle.classList.remove("pt-collapsed");
        // Shiny finishes swapping the new HTML in shortly after this event
        // fires; defer one tick so the checkbox count reflects the NEW list.
        setTimeout(updateIsoformCatalogSummary, 0);
      });
    }
  }
  wireIsoformCatalogCollapse();

  // Generic collapse/expand for the results sections (ui.R's
  // #s1-collapse-body / #s2-collapse-body, matched via a shared
  // data-collapse-target attribute rather than one function per section --
  // unlike the isoform catalog above, there are two of these and more could
  // follow). Both start collapsed (ui.R sets display:none + the triangle's
  // pt-collapsed class inline) since there is nothing to show before the
  // user has run an analysis; expandCollapsible() is called from inside
  // renderSection1()/renderSection2() themselves, so a section only ever
  // pops open at the exact moment it actually has a figure to show --
  // renderSection2() in particular only runs at all when a real confounder
  // target exists, so section 2 correctly stays collapsed (with its "(run
  // analysis to populate)" hint intact) when there's genuinely nothing
  // there, rather than opening onto an empty chart.
  function expandCollapsible(bodyId) {
    const body = document.getElementById(bodyId);
    if (!body) return;
    body.style.display = "";
    const header = document.querySelector('.pt-collapsible-header[data-collapse-target="' + bodyId + '"]');
    if (!header) return;
    const triangle = header.querySelector(".pt-collapse-triangle");
    if (triangle) triangle.classList.remove("pt-collapsed");
    const hint = header.querySelector('[id$="-collapse-hint"]');
    if (hint) hint.textContent = "";
  }

  // Inverse of expandCollapsible() -- puts a results section back to its
  // pre-"Run analysis" collapsed state (used by the Reset button, below):
  // display:none, triangle rotated closed, hint text restored.
  function collapseSection(bodyId, hintText) {
    const body = document.getElementById(bodyId);
    if (body) body.style.display = "none";
    const header = document.querySelector('.pt-collapsible-header[data-collapse-target="' + bodyId + '"]');
    if (!header) return;
    const triangle = header.querySelector(".pt-collapse-triangle");
    if (triangle) triangle.classList.add("pt-collapsed");
    const hint = header.querySelector('[id$="-collapse-hint"]');
    if (hint) hint.textContent = hintText || "";
  }

  function wireResultsCollapse() {
    document.addEventListener("click", e => {
      const header = e.target.closest && e.target.closest(".pt-collapsible-header[data-collapse-target]");
      if (!header) return;
      const body = document.getElementById(header.getAttribute("data-collapse-target"));
      const triangle = header.querySelector(".pt-collapse-triangle");
      if (!body || !triangle) return;
      const collapsed = body.style.display === "none";
      body.style.display = collapsed ? "" : "none";
      triangle.classList.toggle("pt-collapsed", !collapsed);
    });
  }
  wireResultsCollapse();

  // Tears down every pan/zoom/click listener wired to the s1/s2 MS1 charts,
  // ladders, and on-demand fragment-isotope charts -- used by the Reset
  // button (pt_reset_analysis, below). Clearing a chart's innerHTML alone
  // leaves its wheel/drag listeners (attached to the SVG ELEMENT itself,
  // which Reset doesn't remove, only empties) still live: each one's
  // closure still holds the PREVIOUS render's entries/state, so a stray
  // scroll or drag over the now-visually-empty chart called its redraw()
  // and silently repainted the stale pre-reset result right back in.
  // Reproduced directly: Reset, then scroll the mouse wheel over where the
  // MS1 chart used to be -- the old curves reappeared. All the abort
  // machinery already exists per-chart for RE-render (wireMs1ChartOnce/
  // wireInteraction/wireFragmentClick each abort their OWN previous
  // controller before attaching new listeners); this just triggers that
  // same abort immediately, instead of waiting for a render that Reset
  // means may never come.
  function resetInteractions() {
    ["s1-ms1", "s2-ms1"].forEach(idPrefix => {
      const ctrl = ms1ChartAborts[idPrefix];
      if (ctrl) ctrl.abort();
      delete ms1ChartAborts[idPrefix];
    });
    ["s1-ladder", "s2-ladder"].forEach(id => {
      const el = document.getElementById(id);
      if (!el) return;
      const ctrl = ladderInteractionAborts.get(el);
      if (ctrl) ctrl.abort();
      ladderInteractionAborts.delete(el);
    });
    ["s1", "s2"].forEach(sectionPrefix => {
      const sec = fragSections[sectionPrefix];
      if (!sec) return;
      if (sec.listenerAbort) sec.listenerAbort.abort();
      if (sec.zoomAbort) sec.zoomAbort.abort();
      delete fragSections[sectionPrefix];
    });
  }

  return {
    renderSection1, renderSection2, renderExonAlignment, exonAlignmentLegendHtml,
    zoomControlsHtml, wireZoomButtons, updateZoomRangeLabel, attachPanZoom,
    expandCollapsible, collapseSection, resetInteractions
  };
})();

if (window.Shiny) {
  // Persisted across re-renders (isoform/ORF selection changes) so the
  // user's zoom/pan doesn't reset on every tweak -- only when the
  // underlying axis actually changes length does the view snap back to
  // full-width, since an old vs/ve pair could otherwise fall outside a
  // shorter new axis. Keyed by msg.instance so Option 2's (FASTA) and
  // Option 3's (rMATS) independent exon-alignment previews -- which both
  // use this one handler, just pointed at their own DOM element ids --
  // don't clobber each other's zoom/pan state.
  let exonAlignStates = {};
  // AbortController per instance ("fasta"/"rmats"), so a later render (or
  // the Reset button, below) can tear down the PREVIOUS render's pan/zoom
  // listeners on the exon-alignment SVG before attaching new ones -- same
  // stale-closure bug PT.resetInteractions() already fixes for the MS1/
  // ladder charts (see that function's own doc comment), just never
  // extended to this SVG since attachPanZoom() was called here with no
  // opts.signal at all. Reproduced directly: Reset the rMATS panel, then
  // scroll the mouse wheel over where the exon alignment used to be -- the
  // old alignment reappeared, because the wheel listener (still attached to
  // the SVG element, which Reset only empties, not removes) closed over the
  // previous render's now-stale `payload`/`state`/`redraw`.
  let exonAlignAborts = {};
  Shiny.addCustomMessageHandler("pt_render_exon_alignment", function (msg) {
    const instance = msg.instance || "fasta";
    const svgId = msg.svg_id || "fasta-exon-align";
    const legendId = msg.legend_id || "fasta-exon-align-legend";
    const zoomId = msg.zoom_id || "fasta-exon-zoom";
    const idPrefix = msg.id_prefix || "fasta-exon";
    const svgEl = document.getElementById(svgId);
    const legendEl = document.getElementById(legendId);
    const zoomEl = document.getElementById(zoomId);
    if (!svgEl) return;
    // Abort the PREVIOUS render's pan/zoom listeners before doing anything
    // else (including the no-data early return below) -- see
    // exonAlignAborts' own doc comment.
    if (exonAlignAborts[instance]) { exonAlignAborts[instance].abort(); delete exonAlignAborts[instance]; }
    const payload = msg.has_data ? msg.payload : null;
    if (legendEl) legendEl.innerHTML = payload ? PT.exonAlignmentLegendHtml(payload.highlights) : "";
    if (!payload) { PT.renderExonAlignment(svgEl, null); if (zoomEl) zoomEl.innerHTML = ""; exonAlignStates[instance] = null; return; }

    // msg.force_reset_zoom (set by Option 3's rMATS preview, not Option 2's)
    // always snaps back to full view -- each "Find matching transcripts"
    // click is a brand new gene/event, unlike Option 2 where the SAME
    // sequence is being incrementally refined (isoform/ORF re-selection),
    // so preserving zoom there is useful but here it isn't. Without this,
    // a DIFFERENT rMATS run whose axis_length happens to coincide with a
    // previous one (e.g. the same gene, a different event row) would keep
    // showing whatever zoomed/panned sub-range was left over from before,
    // reading as if most of the gene's exons had vanished.
    if (!exonAlignStates[instance] || msg.force_reset_zoom || exonAlignStates[instance].axisLength !== payload.axis_length) {
      exonAlignStates[instance] = { vs: 0, ve: payload.axis_length, axisLength: payload.axis_length };
    }
    const state = exonAlignStates[instance];

    function redraw() {
      PT.renderExonAlignment(svgEl, payload, state);
      PT.updateZoomRangeLabel(idPrefix, state);
    }
    redraw();
    if (zoomEl) {
      zoomEl.innerHTML = PT.zoomControlsHtml(idPrefix);
      PT.wireZoomButtons(idPrefix, state, redraw);
    }
    const controller = new AbortController();
    exonAlignAborts[instance] = controller;
    PT.attachPanZoom(svgEl, state, redraw, { signal: controller.signal });
  });

  // Toggles a button's disabled state -- used to grey out "Run analysis"
  // in middle-down mode until at least one peptide candidate is checked, so
  // a click can't run the analysis on nothing.
  Shiny.addCustomMessageHandler("pt_set_button_enabled", function (msg) {
    const btn = document.getElementById(msg.id);
    if (btn) btn.disabled = !msg.enabled;
  });

  // Greys out/locks every element matched by msg.selector (used for the 3
  // "Input selection" radio inputs, which share name="input_mode" rather
  // than individually-targetable ids, hence a selector-based handler
  // instead of pt_set_button_enabled's single-id one) -- also toggles a
  // class on each input's enclosing <label> so the Bootstrap radio-inline
  // text itself visibly greys out, not just the (small, easy-to-miss)
  // radio circle.
  Shiny.addCustomMessageHandler("pt_set_selector_enabled", function (msg) {
    document.querySelectorAll(msg.selector).forEach(el => {
      el.disabled = !msg.enabled;
      const label = el.closest("label");
      if (label) label.classList.toggle("pt-disabled-label", !msg.enabled);
    });
  });

  // Clears rendered SVG/legend/filter/zoom content left over from a
  // previous "Run analysis"/FASTA-alignment pass -- needed for the reset
  // button: nulling the underlying reactive values on the server stops
  // FUTURE renders, but content a script tag already injected into these
  // elements stays in the DOM until something explicitly clears it.
  // Handler MUST declare exactly one parameter (even though this message
  // carries no real payload) -- Shiny 1.14's addCustomMessageHandler
  // enforces handler.length === 1 and silently throws (aborting the rest
  // of this script block, so THIS handler was simply never registered) for
  // any handler with a different arity. Confirmed directly: the raw
  // WebSocket frame for this message ({"custom":{"pt_reset_analysis":...}})
  // WAS arriving correctly; only the zero-arg handler's registration failed.
  // Live re-render triggered by a scoring-mode switch (no "Run analysis"
  // needed -- see server.R's observeEvent(input$scoring_mode, ...)):
  // tier_b/tier_y/propensity/ms2_stats depend on scoring_mode, but MS1
  // envelopes and exon/PTM/axis data don't, so the server reuses the last
  // full payload and only patches the scoring-dependent fields before
  // resending. A full PT.renderSection1/2() re-render (same call "Run
  // analysis" itself uses) is simplest and correctly re-wires listeners via
  // the existing abort-controller cleanup in wireInteraction/
  // wireFragmentClick -- the one visible tradeoff is that zoom/pan/filter
  // view resets to default on every switch, same as "Run analysis" already
  // does.
  Shiny.addCustomMessageHandler("pt_rescore_section", function (msg) {
    if (msg.section === "s1") PT.renderSection1(msg.payload);
    else if (msg.section === "s2") PT.renderSection2(msg.payload);
  });

  Shiny.addCustomMessageHandler("pt_reset_analysis", function (msg) {
    // Tear down stale pan/zoom/click listeners FIRST -- otherwise a wheel
    // event arriving between the innerHTML clears below and this handler
    // returning could still hit a live (about-to-be-orphaned) listener.
    // See resetInteractions()'s own doc comment for the bug this fixes.
    PT.resetInteractions();
    // Same fix, for the FASTA/rMATS exon-alignment SVGs' own pan/zoom
    // listeners -- see exonAlignAborts' own doc comment. onmousedown is a
    // plain property assignment (attachPanZoom sets svgEl.onmousedown =
    // ..., not addEventListener), so it isn't covered by the abort signal
    // and has to be nulled out directly, or a stray drag (not just a
    // scroll) could still revive the stale alignment the same way.
    ["fasta", "rmats"].forEach(instance => {
      if (exonAlignAborts[instance]) { exonAlignAborts[instance].abort(); delete exonAlignAborts[instance]; }
      exonAlignStates[instance] = null;
    });
    ["fasta-exon-align", "rmats-exon-align"].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.onmousedown = null;
    });
    PT.collapseSection("s1-collapse-body", "(run analysis to populate)");
    PT.collapseSection("s2-collapse-body", "(run analysis to populate)");
    PT.collapseSection("s2-ms1-collapse-body", "");
    PT.collapseSection("s2-ms2-collapse-body", "");
    ["s1-ms1", "s1-ladder", "s2-ms1", "s2-ladder", "s1-frag-ms1", "s2-frag-ms1",
     "fasta-exon-align", "rmats-exon-align"].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.innerHTML = "";
    });
    ["s1-legend", "s1-filters", "s1-zoom", "s1-stats-strip", "s1-frag-ms1-zoom",
     "s2-legend", "s2-filters", "s2-zoom", "s2-stats-strip", "s2-frag-ms1-zoom",
     "fasta-exon-zoom", "fasta-exon-align-legend", "rmats-exon-zoom", "rmats-exon-align-legend"].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.innerHTML = "";
    });
    ["s1-info", "s2-info", "s1-frag-label", "s2-frag-label"].forEach(id => {
      const el = document.getElementById(id);
      if (el) el.textContent = "";
    });
    const tip = document.getElementById("pt-hover-tooltip");
    if (tip) tip.style.display = "none";
  });
}
