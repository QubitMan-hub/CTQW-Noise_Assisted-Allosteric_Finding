"use strict";

/* ---------- small DOM helpers (all text goes through textContent) ---------- */
const $ = (sel) => document.querySelector(sel);
const SVGNS = "http://www.w3.org/2000/svg";

function el(tag, attrs = {}, ...kids) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") node.className = v;
    else if (k === "text") node.textContent = v;
    else node.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat()) if (kid !== null && kid !== undefined && kid !== false)
    node.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return node;
}
function svg(tag, attrs = {}, parent) {
  const node = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "text") node.textContent = v;
    else node.setAttribute(k, v);
  }
  if (parent) parent.appendChild(node);
  return node;
}
const fmt = (v, d = 3) => (v === null || v === undefined || Number.isNaN(v) ? "n/a" : Number(v).toFixed(d));
const fmtInt = (v) => Number(v).toLocaleString("en-US");
const fmtGamma = (g) => String(+Number(g).toPrecision(3));
const pct = (v) => (v === null || v === undefined ? "n/a" : `${(100 * v).toFixed(1)}%`);

/* ---------- tooltip ---------- */
const tip = $("#tooltip");
function showTip(clientX, clientY, head, rows) {
  tip.replaceChildren();
  if (head) tip.append(el("div", { class: "tt-head", text: head }));
  for (const r of rows) {
    const row = el("div", { class: "tt-row" });
    if (r.key) row.append(el("span", { class: "tt-key" + (r.key === "dashed" ? " dashed" : "") }));
    row.append(el("strong", { text: r.value }));
    if (r.label) row.append(el("span", { text: r.label }));
    tip.append(row);
  }
  tip.hidden = false;
  const pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
  let x = clientX + pad, y = clientY - h - pad;
  if (x + w > window.innerWidth - 8) x = clientX - w - pad;
  if (y < 8) y = clientY + pad;
  tip.style.left = `${Math.max(8, x)}px`;
  tip.style.top = `${y}px`;
}
const hideTip = () => { tip.hidden = true; };

/* ---------- input: file or PDB id, exactly one ---------- */
const form = $("#run-form");
const fileInput = $("#file-input");
const dropzone = $("#dropzone");
const pdbInput = $("#pdb-id");
let chosenFile = null;

function setFile(file) {
  chosenFile = file || null;
  dropzone.classList.toggle("has-file", !!chosenFile);
  dropzone.querySelector(".dz-empty").hidden = !!chosenFile;
  dropzone.querySelector(".dz-chosen").hidden = !chosenFile;
  if (chosenFile) {
    $("#dz-file").textContent = chosenFile.name;
    $("#dz-size").textContent = `${(chosenFile.size / 1024).toFixed(0)} KB  ·  click to replace`;
    pdbInput.value = "";
  } else {
    fileInput.value = "";
  }
  clearError();
}
dropzone.addEventListener("click", () => fileInput.click());
fileInput.addEventListener("change", () => setFile(fileInput.files[0]));
["dragenter", "dragover"].forEach((ev) => dropzone.addEventListener(ev, (e) => {
  e.preventDefault(); dropzone.classList.add("drag");
}));
["dragleave", "drop"].forEach((ev) => dropzone.addEventListener(ev, (e) => {
  e.preventDefault(); dropzone.classList.remove("drag");
}));
dropzone.addEventListener("drop", (e) => { if (e.dataTransfer.files[0]) setFile(e.dataTransfer.files[0]); });
pdbInput.addEventListener("input", () => {
  pdbInput.value = pdbInput.value.replace(/[^0-9a-z]/gi, "").toUpperCase();
  if (pdbInput.value && chosenFile) setFile(null);
  clearError();
});
document.querySelectorAll('input[name="labels"]').forEach((r) => r.addEventListener("change", () => {
  $("#custom-labels").hidden = form.labels.value !== "custom";
}));

function showError(msg) { const e = $("#error"); e.textContent = msg; e.hidden = false; }
function clearError() { $("#error").hidden = true; }

/* ---------- run ---------- */
let timer = null;
form.addEventListener("submit", async (e) => {
  e.preventDefault();
  clearError();
  const pdb = pdbInput.value.trim();
  if (!chosenFile && !pdb) return showError("Choose a structure file or type a PDB id.");
  if (chosenFile && pdb) return showError("Give either a structure file or a PDB id, not both.");

  const data = new FormData(form);
  data.delete("structure");
  if (chosenFile) data.append("structure", chosenFile);
  else data.set("pdb_id", pdb);

  const btn = $("#run-btn");
  btn.disabled = true; btn.textContent = "Running";
  $("#results").replaceChildren();
  $("#loading").hidden = false;
  $("#loading").scrollIntoView({ behavior: "smooth", block: "center" });
  const t0 = Date.now();
  const tick = () => {
    const s = Math.floor((Date.now() - t0) / 1000);
    $("#elapsed").textContent = `${String(Math.floor(s / 60)).padStart(2, "0")}:${String(s % 60).padStart(2, "0")}`;
  };
  tick(); timer = setInterval(tick, 1000);

  try {
    const res = await fetch("/api/run", { method: "POST", body: data });
    let body;
    try { body = await res.json(); } catch { body = { error: `The server answered ${res.status} without a result.` }; }
    if (!res.ok || body.error) showError(body.error || "The run failed.");
    else renderResults(body);
  } catch {
    showError("Could not reach the local server. Is web/app.py still running?");
  } finally {
    clearInterval(timer);
    $("#loading").hidden = true;
    btn.disabled = false; btn.textContent = "Run";
  }
});

/* ---------- results ---------- */
let current = null;
const charts = [];
function renderResults(r) {
  current = r;
  charts.length = 0;
  const root = $("#results");
  root.replaceChildren(...[summaryCard(r), r.allosteric ? allostericCard(r) : null, transportCard(r),
    mapCard(r), rankingCard(r), filesCard(r)].filter(Boolean));
  [...root.children].forEach((c, i) => { c.classList.add("reveal"); c.style.animationDelay = `${i * 40}ms`; });
  redraw();
  root.firstElementChild.scrollIntoView({ behavior: "smooth", block: "start" });
}
function redraw() { charts.forEach((c) => drawChart(c)); drawMap(); }

function verdictText(stats, what) {
  if (stats.verdict === "hump") return `${what}: hump at γ = ${fmtGamma(stats.peak_gamma)}`;
  if (stats.verdict === "quantum") return `${what}: best fully quantum`;
  return `${what}: best most dephased`;
}

function summaryCard(r) {
  const s = r.summary, t = r.transport.stats, a = r.allosteric;
  const pill = (st, what) => el("span", { class: "pill " + (st.verdict === "hump" ? "solid" : "outline"), text: verdictText(st, what) });
  const stat = (label, value) => el("div", { class: "stat" },
    el("div", { class: "stat-label", text: label }), el("div", { class: "stat-value", text: value }));
  const notes = r.notes.length ? el("ul", { class: "notes" }, r.notes.map((n) => el("li", { text: n }))) : null;
  return el("section", { class: "card" },
    el("div", { class: "summary-top" },
      el("div", {},
        el("h2", { class: "run-name", text: r.name }),
        el("div", { class: "run-meta", text: `chains ${s.chains.join(", ")}  ·  ran in ${r.elapsed_s} s` })),
      el("div", { class: "pills" }, a ? pill(a.stats, "allosteric AUC") : null, pill(t, "transport"))),
    el("div", { class: "stats" },
      stat("Residues", fmtInt(s.residues)), stat("Contacts", fmtInt(s.contacts)),
      a ? stat("Best AUC", fmt(a.stats.peak_value)) : stat("Transport peak γ", fmtGamma(t.peak_gamma)),
      a ? stat("Best AUC at γ", fmtGamma(a.stats.peak_gamma)) : stat("Gain over ends", t.verdict === "hump" ? pct(t.gain_rel) : "none")),
    el("dl", { class: "facts" },
      el("dt", { text: "Labels" }), el("dd", { text: r.label_note }),
      el("dt", { text: "Transport" }),
      el("dd", {}, "from ", el("span", { class: s.walk_from === "active site" ? "" : "mono", text: s.walk_from }),
        ` to the ${s.distal_count} residues at least ${s.distal_min_hops} contacts away`),
      el("dt", { text: "Site energies" }),
      el("dd", {}, el("span", { class: "mono", text: `${s.site_energy}, scale ${s.scale}` })),
    ),
    notes);
}

function download(href, label, ext, dark = false) {
  return el("a", { class: "btn" + (dark ? " dark" : ""), href, download: "" },
    el("span", { text: label }), ext ? el("span", { class: "ext", text: ext }) : null);
}

function statsRow(st, unit) {
  const items = [["quantum (γ = 0)", fmt(st.quantum_value)], [`best (γ = ${fmtGamma(st.peak_gamma)})`, fmt(st.peak_value)],
    ["most dephased", fmt(st.dephased_value)]];
  if (st.verdict === "hump") {
    items.push(["gain over better end", `+${fmt(st.gain_abs, unit === "auc" ? 3 : 4)} (${pct(st.gain_rel)})`]);
    items.push(["noise helps for γ", `${fmtGamma(st.helps_range[0])} – ${fmtGamma(st.helps_range[1])}`]);
  }
  return el("div", { class: "metric-row" }, items.map(([k, v]) => el("span", {}, k, el("strong", { text: v }))));
}

function controlNote(ctrl, st) {
  if (!ctrl) return null;
  const humps = ctrl.seed_verdicts ? `${ctrl.seed_verdicts.filter((v) => v === "hump").length} of ${ctrl.seed_verdicts.length} seeds show a hump; ` : "";
  return el("p", { class: "note" }, "Random-energy control: ",
    el("span", { class: "mono", text: verdictText(ctrl.stats, "mean curve").toLowerCase() }), `. ${humps}`,
    `the real curve's gain exceeds ${ctrl.exceeds_seeds} of ${ctrl.n_seeds} seed gains`,
    ctrl.z !== null ? ` (z = ${fmt(ctrl.z, 2)}).` : ".",
    st.verdict === "hump" && ctrl.exceeds_seeds < ctrl.n_seeds
      ? " A hump that random energies reproduce is not specific to the hydropathy model." : "");
}

function allostericCard(r) {
  const a = r.allosteric, st = a.stats;
  const series = [{ values: a.auc, label: `${r.summary.site_energy} site energies`, style: "main" }];
  if (a.control) series.push({ values: a.control.mean, sd: a.control.sd, label: "random energies (mean ± sd)", style: "control" });
  const chart = { id: "chart-allo", gammas: r.gammas, series, yLabel: "ROC AUC, allosteric residues", peakGamma: st.peak_gamma,
    baselines: Object.entries(a.baselines).map(([k, v]) => ({ value: v, label: k })), chance: 0.5, auc: true, digits: 3 };
  charts.push(chart);
  const miss = [...a.missing.active, ...a.missing.allosteric];
  return el("section", { class: "card" },
    el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Allosteric recovery versus noise" }),
        el("p", { class: "card-sub", text: `Walk from the ${a.n_active} active-site residues; how well the signal each residue receives ranks the ${a.n_allosteric} known allosteric residues among ${a.n_candidates} candidates (ROC AUC; 0.5 is chance).` }))),
    legend(series, a.baselines),
    el("div", { class: "chart-wrap", id: chart.id }),
    statsRow(st, "auc"),
    el("div", { class: "metric-row" }, Object.entries(a.baselines).map(([k, v]) => el("span", {}, `baseline: ${k}`, el("strong", { text: fmt(v) })))),
    controlNote(a.control, st),
    miss.length || a.name_mismatch.length ? el("p", { class: "note", text:
      `${miss.length ? `Not in the network: ${miss.join(", ")}. ` : ""}${a.name_mismatch.length ? `Name mismatch: ${a.name_mismatch.join("; ")}.` : ""}` }) : null,
    el("div", { class: "btn-row" },
      download(r.files.allosteric_png, "Figure for paper", "PNG 300 dpi", true),
      download(r.files.allosteric_svg, "Figure", "SVG"),
      download(r.files.allosteric_csv, "Data", "CSV")),
    dataTable(r.gammas, series, 3),
    r.parameters.labels_used && r.parameters.labels_used.citation
      ? el("p", { class: "cite", text: `Labels: ${r.parameters.labels_used.citation}` }) : null);
}

function transportCard(r) {
  const t = r.transport, st = t.stats, s = r.summary;
  const series = [{ values: t.distal_mean, label: `${s.site_energy} site energies`, style: "main" }];
  if (t.control) series.push({ values: t.control.mean, sd: t.control.sd, label: "random energies (mean ± sd)", style: "control" });
  const chart = { id: "chart-transport", gammas: r.gammas, series, yLabel: "signal reaching distal residues", peakGamma: st.peak_gamma,
    baselines: [], chance: null, auc: false, digits: 4 };
  charts.push(chart);
  return el("section", { class: "card" },
    el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Transport versus noise" }),
        el("p", { class: "card-sub" }, "Mean signal reaching the far part of the network from ",
          el("span", { class: s.walk_from === "active site" ? "" : "mono", text: s.walk_from }), ` across ${r.gammas.length} dephasing rates`))),
    series.length > 1 ? legend(series, {}) : null,
    el("div", { class: "chart-wrap", id: chart.id }),
    statsRow(st, "transport"),
    controlNote(t.control, st),
    el("div", { class: "btn-row" },
      download(r.files.hump_png, "Figure for paper", "PNG 300 dpi", !r.allosteric),
      download(r.files.hump_svg, "Figure", "SVG"),
      download(r.files.hump_csv, "Data", "CSV")),
    dataTable(r.gammas, series, 4));
}

function legend(series, baselines) {
  return el("div", { class: "legend" },
    series.map((s) => el("span", { class: "legend-item" }, el("span", { class: "key-line" + (s.style === "control" ? " dashed" : "") }), el("span", { text: s.label }))),
    Object.keys(baselines).length ? el("span", { class: "legend-item" }, el("span", { class: "key-line dotted" }), el("span", { text: "baselines" })) : null);
}

function dataTable(gammas, series, digits) {
  const head = el("tr", {}, el("th", { text: "γ" }), series.map((s) => el("th", { class: "num", text: s.style === "control" ? "control mean" : "value" })),
    series.some((s) => s.sd) ? el("th", { class: "num", text: "control sd" }) : null);
  const rows = gammas.map((g, i) => el("tr", {}, el("td", { text: fmtGamma(g) }),
    series.map((s) => el("td", { class: "num", text: fmt(s.values[i], digits) })),
    series.filter((s) => s.sd).map((s) => el("td", { class: "num", text: fmt(s.sd[i], digits) }))));
  return el("details", { class: "fold" }, el("summary", { text: "Data table" }),
    el("div", { class: "table-scroll" }, el("table", { class: "data-table" }, el("thead", {}, head), el("tbody", {}, rows))));
}

/* line chart; x on a log(1 + γ/0.05) axis so γ = 0 sits at the left edge */
const gx = (g) => Math.log10(1 + g / 0.05);
function niceAxis(lo, hi) {                // round step (1, 2 or 5 x 10^k), about 4-6 ticks
  const raw = (hi - lo) / 4, p = 10 ** Math.floor(Math.log10(raw)), n = raw / p;
  const step = (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * p;
  const decimals = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return { step, min: Math.floor(lo / step + 1e-9) * step, max: Math.ceil(hi / step - 1e-9) * step, decimals };
}

function drawChart(c) {
  const host = document.getElementById(c.id);
  if (!host) return;
  host.replaceChildren();
  const W = Math.max(300, host.clientWidth), H = Math.round(Math.min(340, Math.max(220, W * 0.46)));
  const m = { l: 64, r: c.baselines.length ? 112 : 16, t: 24, b: 44 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const gammas = c.gammas, gmax = gammas[gammas.length - 1];
  const vals = c.series.flatMap((s) => s.values.map((v, i) => [v - (s.sd ? s.sd[i] : 0), v + (s.sd ? s.sd[i] : 0)]).flat())
    .concat(c.baselines.map((b) => b.value)).concat(c.chance !== null ? [c.chance] : []).filter((v) => v !== null);
  let lo = c.auc ? Math.min(...vals) - 0.03 : 0, hi = Math.max(...vals) * (c.auc ? 1 : 1.08) + (c.auc ? 0.03 : 0);
  if (c.auc) { lo = Math.max(0, lo); hi = Math.min(1, hi); }
  const axis = niceAxis(lo, hi);
  const X = (g) => m.l + (gx(g) / gx(gmax)) * iw;
  const Y = (v) => m.t + ih - ((v - axis.min) / (axis.max - axis.min)) * ih;
  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img", "aria-label": c.yLabel + " versus dephasing rate" }, host);

  for (let v = axis.min; v <= axis.max + axis.step / 2; v += axis.step) {
    const y = Y(v);
    svg("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: Math.abs(v - axis.min) < 1e-12 ? "#d4d4d4" : "#f0f0f0", "stroke-width": 1, "shape-rendering": "crispEdges" }, root);
    svg("text", { x: m.l - 10, y: y + 3.5, "text-anchor": "end", "font-size": 11, fill: "#6b6b6b", text: v.toFixed(axis.decimals) }, root);
  }
  for (const t of [0, 0.1, 1, 10, 100].filter((t) => t <= gmax * 1.01)) {
    svg("line", { x1: X(t), x2: X(t), y1: m.t + ih, y2: m.t + ih + 4, stroke: "#c4c4c4", "shape-rendering": "crispEdges" }, root);
    svg("text", { x: X(t), y: m.t + ih + 18, "text-anchor": "middle", "font-size": 11, fill: "#6b6b6b", text: fmtGamma(t) }, root);
  }
  svg("text", { x: m.l + iw / 2, y: H - 6, "text-anchor": "middle", "font-size": 11.5, fill: "#3a3a3a", text: "dephasing rate γ  (log scale)" }, root);
  svg("text", { x: 12, y: m.t + ih / 2, "text-anchor": "middle", "font-size": 11.5, fill: "#3a3a3a",
    transform: `rotate(-90 12 ${m.t + ih / 2})`, text: c.yLabel }, root);
  if (c.chance !== null) svg("line", { x1: m.l, x2: W - m.r, y1: Y(c.chance), y2: Y(c.chance), stroke: "#cfcfcf", "stroke-width": 1 }, root);
  for (const b of c.baselines) {
    svg("line", { x1: m.l, x2: W - m.r, y1: Y(b.value), y2: Y(b.value), stroke: "#8a8a8a", "stroke-width": 1, "stroke-dasharray": "1.5 3" }, root);
    svg("text", { x: W - m.r + 6, y: Y(b.value) + 3.5, "font-size": 10.5, fill: "#6b6b6b", text: `${{ "proximity to active site": "proximity", "contact degree": "degree" }[b.label] || b.label} ${fmt(b.value)}` }, root);
  }
  const pg = c.peakGamma;
  svg("line", { x1: X(pg), x2: X(pg), y1: m.t - 6, y2: m.t + ih, stroke: "#a3a3a3", "stroke-width": 1, "stroke-dasharray": "3 3" }, root);
  svg("text", { x: X(pg) + (X(pg) > W - m.r - 110 ? -6 : 6), y: m.t - 10, "text-anchor": X(pg) > W - m.r - 110 ? "end" : "start",
    "font-size": 11, fill: "#6b6b6b", text: `best γ = ${fmtGamma(pg)}` }, root);

  const main = c.series[0];
  const pts = (vals) => gammas.map((g, i) => [X(g), Y(vals[i])]);
  if (!c.auc) {
    const mp = pts(main.values);
    svg("path", { d: `M${mp[0][0]},${Y(axis.min)} ` + mp.map((p) => `L${p[0]},${p[1]}`).join(" ") + ` L${mp[mp.length - 1][0]},${Y(axis.min)} Z`,
      fill: "#161616", "fill-opacity": 0.05 }, root);
  }
  const dots = [];
  [...c.series].reverse().forEach((s) => {
    const isMain = s === main;
    if (s.sd) {
      const up = gammas.map((g, i) => [X(g), Y(s.values[i] + s.sd[i])]), dn = gammas.map((g, i) => [X(g), Y(s.values[i] - s.sd[i])]).reverse();
      svg("path", { d: "M" + up.concat(dn).map((q) => q.join(",")).join(" L") + " Z", fill: "#8a8a8a", "fill-opacity": 0.16 }, root);
    }
    const p = pts(s.values);
    svg("path", { d: "M" + p.map((q) => q.join(",")).join(" L"), fill: "none", stroke: isMain ? "#161616" : "#8a8a8a",
      "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round", "stroke-dasharray": isMain ? "none" : "5 4" }, root);
    if (isMain) dots.push(p.map((q) => svg("circle", { cx: q[0], cy: q[1], r: 3.2, fill: "#161616", stroke: "#fff", "stroke-width": 1.5 }, root)));
  });

  const cross = svg("line", { y1: m.t, y2: m.t + ih, stroke: "#161616", "stroke-opacity": 0.25, "stroke-width": 1, visibility: "hidden" }, root);
  const overlay = svg("rect", { x: m.l, y: m.t - 8, width: iw, height: ih + 8, fill: "transparent", class: "chart-overlay", tabindex: 0,
    "aria-label": "Use left and right arrows to read values" }, root);
  let active = -1;
  function setActive(i, cx, cy) {
    active = i;
    dots.forEach((set) => set.forEach((d, j) => d.setAttribute("r", j === i ? 5 : 3.2)));
    if (i < 0) { cross.setAttribute("visibility", "hidden"); hideTip(); return; }
    cross.setAttribute("x1", X(gammas[i])); cross.setAttribute("x2", X(gammas[i]));
    cross.setAttribute("visibility", "visible");
    const box = root.getBoundingClientRect();
    const px = cx ?? box.left + (X(gammas[i]) / W) * box.width, py = cy ?? box.top + (m.t / H) * box.height + 40;
    showTip(px, py, `γ = ${fmtGamma(gammas[i])}`, c.series.map((s) => ({ key: s.style === "control" ? "dashed" : "solid",
      value: fmt(s.values[i], c.digits) + (s.sd ? ` ± ${fmt(s.sd[i], c.digits)}` : ""), label: s.style === "control" ? "control" : r0(c) })));
  }
  overlay.addEventListener("pointermove", (e) => {
    const box = root.getBoundingClientRect(), x = ((e.clientX - box.left) / box.width) * W;
    let best = 0, bd = Infinity;
    gammas.forEach((g, i) => { const d = Math.abs(X(g) - x); if (d < bd) { bd = d; best = i; } });
    setActive(best, e.clientX, e.clientY);
  });
  overlay.addEventListener("pointerleave", () => setActive(-1));
  overlay.addEventListener("blur", () => setActive(-1));
  overlay.addEventListener("keydown", (e) => {
    if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return;
    e.preventDefault();
    const start = gammas.indexOf(pg);
    setActive(active < 0 ? start : Math.min(gammas.length - 1, Math.max(0, active + (e.key === "ArrowRight" ? 1 : -1))));
  });
}
const r0 = (c) => (c.auc ? "AUC" : "signal");

/* ---------- residue signal map ---------- */
function mapCard(r) {
  const mp = r.map, others = mp.nodes.filter((n) => n.role !== "source").map((n) => n.score);
  const lo = Math.min(...others), hi = Math.max(...others);
  const ringKey = (kind) => {
    const s = svg("svg", { width: 18, height: 18, viewBox: "0 0 18 18" });
    svg("circle", { cx: 9, cy: 9, r: 4, fill: "#8a8a8a" }, s);
    svg("circle", { cx: 9, cy: 9, r: 7.5, fill: "none", stroke: kind === "source" ? "#161616" : "#6b6b6b", "stroke-width": 1.5,
      "stroke-dasharray": kind === "source" ? "none" : "2.5 2" }, s);
    return s;
  };
  const fromActive = mp.from === "active site";
  return el("section", { class: "card" },
    el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Residue signal map" }),
        el("p", { class: "card-sub" }, `Signal each residue receives from the ${fromActive ? "active site" : "source"} at `,
          el("span", { class: "mono", text: `γ = ${fmtGamma(mp.gamma)}` }),
          `${r.allosteric ? " (the best-AUC rate)" : " (the transport peak)"}. Darker means more signal. Residues sit at their 3D positions projected onto the protein's two main axes; lines are contacts.`))),
    el("div", { class: "map-legend" },
      el("span", { class: "ramp" }, el("span", { class: "mono", text: fmt(lo) }), el("span", { class: "ramp-bar" }), el("span", { class: "mono", text: fmt(hi) })),
      el("span", { class: "ring-key" }, ringKey("source"), el("span", { text: fromActive ? "active site" : `start ${r.summary.walk_from}` })),
      r.allosteric ? el("span", { class: "ring-key" }, ringKey("other"), el("span", { text: "known allosteric" })) : null),
    el("div", { class: "map-wrap", id: "map" }));
}

function shade(t) {                     // light grey -> near-black
  const v = Math.round(236 + (22 - 236) * Math.min(1, Math.max(0, t)));
  return `rgb(${v},${v},${v})`;
}

function drawMap() {
  const r = current, host = $("#map");
  if (!r || !host) return;
  host.replaceChildren();
  const mp = r.map, N = mp.nodes.length;
  const W = Math.max(300, host.clientWidth), pad = 22, inner = W - 2 * pad;
  const H = Math.round(Math.min(inner * Math.max(mp.aspect, 0.35), 560) + 2 * pad);
  const sx = Math.min(inner, (H - 2 * pad) / Math.max(mp.aspect, 1e-3));
  const offX = pad + (inner - sx) / 2;
  const P = mp.nodes.map((n) => [offX + n.x * sx, pad + n.y * sx]);
  const others = mp.nodes.filter((n) => n.role !== "source").map((n) => n.score);
  const lo = Math.min(...others), hi = Math.max(...others);
  const rad = Math.max(3, Math.min(7, 150 / Math.sqrt(N)));
  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
    "aria-label": `Contact network of ${N} residues shaded by signal received` }, host);
  const eg = svg("g", {}, root);
  const edgeEls = mp.edges.map(([a, b]) => svg("line", { x1: P[a][0], y1: P[a][1], x2: P[b][0], y2: P[b][1], stroke: "#e6e6e6", "stroke-width": 1 }, eg));
  const edgesOf = mp.nodes.map(() => []);
  mp.edges.forEach(([a, b], k) => { edgesOf[a].push(k); edgesOf[b].push(k); });

  const ng = svg("g", {}, root);
  const order = mp.nodes.map((_, i) => i).sort((a, b) => mp.nodes[a].score - mp.nodes[b].score);
  const nodeEls = [];
  for (const i of order) {
    const n = mp.nodes[i];
    const t = n.role === "source" ? 1 : (n.score - lo) / (hi - lo || 1);
    nodeEls[i] = svg("circle", { cx: P[i][0], cy: P[i][1], r: rad, fill: shade(t), stroke: "#fff", "stroke-width": 1.5 }, ng);
  }
  const rg = svg("g", {}, root);
  const ring = (i, dashed) => svg("circle", { cx: P[i][0], cy: P[i][1], r: rad + 4, fill: "none", stroke: dashed ? "#6b6b6b" : "#161616",
    "stroke-width": 1.4, "stroke-dasharray": dashed ? "3 2.5" : "none" }, rg);
  mp.nodes.forEach((n, i) => { if (n.role === "source") ring(i, false); else if (n.role === "allosteric") ring(i, true); });

  const hit = svg("rect", { x: 0, y: 0, width: W, height: H, fill: "transparent" }, root);
  let active = -1;
  function setActive(i, e) {
    if (active >= 0) {
      nodeEls[active].setAttribute("r", rad); nodeEls[active].setAttribute("stroke", "#fff");
      edgesOf[active].forEach((k) => edgeEls[k].setAttribute("stroke", "#e6e6e6"));
    }
    active = i;
    if (i < 0) { hideTip(); return; }
    nodeEls[i].setAttribute("r", rad + 1.5); nodeEls[i].setAttribute("stroke", "#161616");
    edgesOf[i].forEach((k) => edgeEls[k].setAttribute("stroke", "#9a9a9a"));
    const n = mp.nodes[i];
    const role = n.role === "source" ? (mp.from === "active site" ? "  ·  active site" : "  ·  source")
      : n.role === "allosteric" ? "  ·  known allosteric" : "";
    showTip(e.clientX, e.clientY, `${n.id}  ${n.resname}${role}`,
      [{ value: fmt(n.score, 4), label: "signal" }, { value: String(n.degree), label: "contacts" }]);
  }
  hit.addEventListener("pointermove", (e) => {
    const box = root.getBoundingClientRect();
    const x = ((e.clientX - box.left) / box.width) * W, y = ((e.clientY - box.top) / box.height) * H;
    let best = -1, bd = 14 * 14;
    P.forEach((p, i) => { const d = (p[0] - x) ** 2 + (p[1] - y) ** 2; if (d < bd) { bd = d; best = i; } });
    setActive(best, e);
  });
  hit.addEventListener("pointerleave", () => setActive(-1));
}

/* ---------- ranked table ---------- */
function rankingCard(r) {
  const top = r.ranking, max = Math.max(...top.map((x) => x.score));
  const labelled = !!r.allosteric;
  const hits = labelled ? top.filter((x) => x.known_allosteric).length : 0;
  const rows = top.map((x) => el("tr", {},
    el("td", { class: "rank", text: String(x.rank) }),
    el("td", { text: x.id }),
    el("td", { text: x.resname }),
    el("td", { class: "num" }, el("div", { class: "bar-cell" },
      el("span", { text: fmt(x.score, 4) }),
      el("span", { class: "bar" }, el("span", { style: `width:${(100 * x.score / max).toFixed(1)}%` })))),
    el("td", { class: "num", text: String(x.degree) }),
    labelled ? el("td", { class: "known" }, x.known_allosteric ? el("span", { class: "dot-known", title: "known allosteric" }) : null) : null));
  return el("details", { class: "card card-fold", open: true },
    el("summary", {}, el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Top residues by signal" }),
        el("p", { class: "card-sub" }, `Top ${top.length} of ${r.ranking_total} at `, el("span", { class: "mono", text: `γ = ${fmtGamma(r.map.gamma)}` }),
          `, ${r.summary.walk_from === "active site" ? "active-site residues" : "start residue"} excluded${labelled ? ` · ${hits} of these are known allosteric residues` : ""}`)),
      el("span", { class: "chev", "aria-hidden": "true" }))),
    el("div", { class: "table-scroll" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", { text: "#" }), el("th", { text: "Residue" }), el("th", { text: "Name" }),
        el("th", { class: "num", text: "Signal" }), el("th", { class: "num", text: "Contacts" }), labelled ? el("th", { text: "Known" }) : null)),
      el("tbody", {}, rows))),
    el("div", { class: "btn-row" }, download(r.files.ranking_csv, `Full ranking, ${r.ranking_total} residues`, "CSV")));
}

/* ---------- downloads + reproducibility ---------- */
function filesCard(r) {
  const f = r.files, name = r.name;
  const item = (href, fname, desc) => href ? el("li", {},
    el("div", {}, el("div", { class: "fname", text: fname }), el("div", { class: "fdesc", text: desc })),
    download(href, "Download")) : null;
  const v = r.parameters.code_version;
  const ver = v.git_commit ? `git ${v.git_commit.slice(0, 12)}${v.git_dirty ? " (pipeline files modified)" : ""}  ·  pipeline sha256 ${v.pipeline_sha256}`
    : `pipeline sha256 ${v.pipeline_sha256}`;
  return el("details", { class: "card card-fold", open: true },
    el("summary", {}, el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Files and reproducibility" }),
        el("p", { class: "card-sub", text: "Everything this run produced, plus the exact settings used." })),
      el("span", { class: "chev", "aria-hidden": "true" }))),
    el("ul", { class: "files" },
      item(f.graphml, `${name}.graphml`, "Residue contact network"),
      item(f.result, `${name}_result.json`, "Every number behind the figures"),
      item(f.parameters, `${name}_parameters.json`, "Every setting used, the labels and a code version marker")),
    el("div", { class: "version", text: ver }),
    el("details", { class: "fold" }, el("summary", { text: "Show parameters" }),
      el("pre", { class: "params", text: JSON.stringify(r.parameters, null, 2) })));
}

/* redraw the SVGs at the new width */
let rs = null;
window.addEventListener("resize", () => { clearTimeout(rs); rs = setTimeout(() => { hideTip(); redraw(); }, 120); });
