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
const fmt = (v, d = 3) => Number(v).toFixed(d);
const fmtInt = (v) => Number(v).toLocaleString("en-US");
const fmtGamma = (g) => String(+g);

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
function renderResults(r) {
  current = r;
  const root = $("#results");
  root.replaceChildren(
    summaryCard(r), chartCard(r), mapCard(r), rankingCard(r), filesCard(r),
  );
  [...root.children].forEach((c, i) => { c.classList.add("reveal"); c.style.animationDelay = `${i * 40}ms`; });
  drawChart(); drawMap();
  root.firstElementChild.scrollIntoView({ behavior: "smooth", block: "start" });
}

function verdictText(v, gamma) {
  if (v === "hump") return `Hump at γ = ${fmtGamma(gamma)}`;
  if (v === "quantum") return "No hump · best fully quantum";
  return "No hump · best at the classical end";
}

function qmodText(r) {
  const q = r.qmod, s = r.summary;
  switch (q.status) {
    case "written": return `written · ${fmtInt(q.pauli_terms)} Pauli terms · ${fmtInt(q.size_kb)} KB`;
    case "skipped_too_large": return `Qmod skipped, protein too large (${s.residues} residues > ${s.max_residues}). The graph is unaffected.`;
    case "skipped_by_user": return "skipped (graph only)";
    case "unavailable": return "Qmod unavailable: the classiq package is not installed. The graph is unaffected.";
    default: return `Qmod could not be written (${q.error || "unknown error"}). The graph is unaffected.`;
  }
}

function summaryCard(r) {
  const s = r.summary;
  const pill = el("span", { class: "pill " + (s.verdict === "hump" ? "solid" : "outline"),
    text: verdictText(s.verdict, s.peak_gamma) });
  const stat = (label, value) => el("div", { class: "stat" },
    el("div", { class: "stat-label", text: label }), el("div", { class: "stat-value", text: value }));
  const qmodDd = el("dd", {}, el("span", { class: q(r) ? "mono" : "", text: qmodText(r) }));
  return el("section", { class: "card" },
    el("div", { class: "summary-top" },
      el("div", {},
        el("h2", { class: "run-name", text: r.name }),
        el("div", { class: "run-meta", text: `chains ${s.chains.join(", ")}  ·  ${fmtInt(s.contacts)} contacts  ·  ran in ${r.elapsed_s} s` })),
      pill),
    el("div", { class: "stats" },
      stat("Residues", fmtInt(s.residues)), stat("Qubits", String(s.qubits)),
      stat("Peak γ", fmtGamma(s.peak_gamma)), stat("Peak signal", fmt(s.peak_signal))),
    el("dl", { class: "facts" },
      el("dt", { text: "Walk" }),
      el("dd", {}, el("span", { class: "mono", text: s.source }), el("span", { class: "arrow", text: "→" }),
        el("span", { class: "mono", text: s.target }),
        el("span", { class: "tag", text: (s.source_auto ? "source auto: most connected" : "source set by you") + ", target auto: farthest" })),
      el("dt", { text: "Site energies" }),
      el("dd", {}, el("span", { class: "mono", text: `${s.site_energy}, scale ${s.scale}` })),
      el("dt", { text: "Qmod" }), qmodDd));
}
const q = (r) => r.qmod.status === "written" || r.qmod.status === "skipped_by_user";

function download(href, label, ext, dark = false) {
  return el("a", { class: "btn" + (dark ? " dark" : ""), href, download: "" },
    el("span", { text: label }), ext ? el("span", { class: "ext", text: ext }) : null);
}

function chartCard(r) {
  const s = r.summary;
  const card = el("section", { class: "card" },
    el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Transport versus noise" }),
        el("p", { class: "card-sub" }, "Signal reaching ", el("span", { class: "mono", text: s.target }),
          " from ", el("span", { class: "mono", text: s.source }), ` across ${r.gammas.length} dephasing rates`))));
  if (r.curves.length > 1) {
    card.append(el("div", { class: "legend" },
      el("span", { class: "legend-item" }, el("span", { class: "key-line" }), el("span", { text: r.curves[0].label })),
      el("span", { class: "legend-item" }, el("span", { class: "key-line dashed" }), el("span", { text: r.curves[1].label }))));
  }
  card.append(el("div", { class: "chart-wrap", id: "chart" }));
  if (r.curves.length > 1) {
    const c = r.curves[1];
    card.append(el("p", { class: "note" }, "Random control: ",
      el("span", { class: "mono", text: verdictText(c.verdict, r.gammas[c.peak_index]).toLowerCase() }),
      `, peak ${fmt(c.values[c.peak_index])}. If the control also peaks at intermediate noise, the hump is not specific to the hydropathy energies.`));
  }
  card.append(el("div", { class: "btn-row" },
    download(r.files.figure_png, "Figure for paper", "PNG 300 dpi", true),
    download(r.files.figure_svg, "Figure", "SVG"),
    download(r.files.hump_csv, "Curve data", "CSV"),
    r.files.control_csv ? download(r.files.control_csv, "Control data", "CSV") : null));

  const head = el("tr", {}, el("th", { text: "γ" }), ...r.curves.map((c) => el("th", { class: "num", text: c.key === "random" && r.curves.length > 1 ? "control" : "signal" })));
  const rows = r.gammas.map((g, i) => el("tr", {}, el("td", { text: fmtGamma(g) }),
    ...r.curves.map((c) => el("td", { class: "num", text: fmt(c.values[i], 4) }))));
  card.append(el("details", { class: "fold" }, el("summary", { text: "Data table" }),
    el("div", { class: "table-scroll" }, el("table", { class: "data-table" }, el("thead", {}, head), el("tbody", {}, rows)))));
  return card;
}

/* hump chart: x on a log(1 + γ/0.05) axis so the low-noise points are readable */
const gx = (g) => Math.log10(1 + g / 0.05);
function niceAxis(v) {                 // round step (1, 2 or 5 x 10^k) giving about 4-5 ticks
  const raw = v / 4, p = 10 ** Math.floor(Math.log10(raw)), n = raw / p;
  const step = (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * p;
  const decimals = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
  return { step, max: Math.ceil(v / step - 1e-9) * step, decimals };
}

function drawChart() {
  const r = current, host = $("#chart");
  if (!r || !host) return;
  host.replaceChildren();
  const W = Math.max(300, host.clientWidth), H = Math.round(Math.min(340, Math.max(220, W * 0.46)));
  const m = { l: 64, r: 16, t: 24, b: 44 };
  const iw = W - m.l - m.r, ih = H - m.t - m.b;
  const gammas = r.gammas, gmax = gammas[gammas.length - 1];
  const all = r.curves.flatMap((c) => c.values);
  const axis = niceAxis(Math.max(...all) * 1.08), yMax = axis.max;
  const X = (g) => m.l + (gx(g) / gx(gmax)) * iw;
  const Y = (v) => m.t + ih - (v / yMax) * ih;

  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
    "aria-label": "Signal reaching the target residue versus dephasing rate" }, host);

  // faint horizontal gridlines + y ticks
  const yTicks = Math.round(yMax / axis.step);
  for (let i = 0; i <= yTicks; i++) {
    const v = axis.step * i, y = Y(v);
    svg("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: i === 0 ? "#d4d4d4" : "#f0f0f0", "stroke-width": 1, "shape-rendering": "crispEdges" }, root);
    svg("text", { x: m.l - 10, y: y + 3.5, "text-anchor": "end", "font-size": 11, fill: "#6b6b6b", text: v.toFixed(axis.decimals) }, root);
  }
  // x ticks
  for (const t of [0, 0.1, 0.5, 1, 5, 10, 50].filter((t) => t <= gmax)) {
    svg("line", { x1: X(t), x2: X(t), y1: m.t + ih, y2: m.t + ih + 4, stroke: "#c4c4c4", "shape-rendering": "crispEdges" }, root);
    svg("text", { x: X(t), y: m.t + ih + 18, "text-anchor": "middle", "font-size": 11, fill: "#6b6b6b", text: fmtGamma(t) }, root);
  }
  svg("text", { x: m.l + iw / 2, y: H - 6, "text-anchor": "middle", "font-size": 11.5, fill: "#3a3a3a", text: "dephasing rate γ  (log scale)" }, root);
  svg("text", { x: 12, y: m.t + ih / 2, "text-anchor": "middle", "font-size": 11.5, fill: "#3a3a3a",
    transform: `rotate(-90 12 ${m.t + ih / 2})`, text: "signal reaching target" }, root);

  // peak marker
  const main = r.curves[0], pg = gammas[main.peak_index];
  svg("line", { x1: X(pg), x2: X(pg), y1: m.t - 6, y2: m.t + ih, stroke: "#a3a3a3", "stroke-width": 1, "stroke-dasharray": "3 3" }, root);
  svg("text", { x: X(pg) + (X(pg) > W - 110 ? -6 : 6), y: m.t - 10, "text-anchor": X(pg) > W - 110 ? "end" : "start",
    "font-size": 11, fill: "#6b6b6b", text: `peak γ = ${fmtGamma(pg)}` }, root);

  // area under the main curve
  const pts = (c) => gammas.map((g, i) => [X(g), Y(c.values[i])]);
  const mp = pts(main);
  svg("path", { d: `M${mp[0][0]},${Y(0)} ` + mp.map((p) => `L${p[0]},${p[1]}`).join(" ") + ` L${mp[mp.length - 1][0]},${Y(0)} Z`,
    fill: "#161616", "fill-opacity": 0.06 }, root);

  // lines (control first so the main line sits on top)
  const dots = [];
  [...r.curves].reverse().forEach((c) => {
    const isMain = c === main, p = pts(c);
    svg("path", { d: "M" + p.map((q) => q.join(",")).join(" L"), fill: "none", stroke: isMain ? "#161616" : "#8a8a8a",
      "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round", "stroke-dasharray": isMain ? "none" : "5 4" }, root);
    const set = p.map((q) => svg("circle", { cx: q[0], cy: q[1], r: 4, fill: isMain ? "#161616" : "#8a8a8a", stroke: "#fff", "stroke-width": 2 }, root));
    dots.push({ c, set });
  });

  // crosshair + tooltip, snapping to the nearest noise level
  const cross = svg("line", { y1: m.t, y2: m.t + ih, stroke: "#161616", "stroke-opacity": 0.25, "stroke-width": 1, visibility: "hidden" }, root);
  const overlay = svg("rect", { x: m.l, y: m.t - 8, width: iw, height: ih + 8, fill: "transparent", class: "chart-overlay", tabindex: 0,
    "aria-label": "Use left and right arrows to read values" }, root);
  let active = -1;
  function setActive(i, cx, cy) {
    active = i;
    dots.forEach(({ set }) => set.forEach((d, j) => d.setAttribute("r", j === i ? 5.5 : 4)));
    if (i < 0) { cross.setAttribute("visibility", "hidden"); hideTip(); return; }
    cross.setAttribute("x1", X(gammas[i])); cross.setAttribute("x2", X(gammas[i]));
    cross.setAttribute("visibility", "visible");
    const box = root.getBoundingClientRect();
    const px = cx ?? box.left + (X(gammas[i]) / W) * box.width, py = cy ?? box.top + (m.t / H) * box.height + 40;
    showTip(px, py, `γ = ${fmtGamma(gammas[i])}`,
      r.curves.map((c, k) => ({ key: k === 0 ? "solid" : "dashed", value: fmt(c.values[i], 4), label: r.curves.length > 1 ? (k === 0 ? c.key : "control") : "signal" })));
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
    const next = active < 0 ? main.peak_index : Math.min(gammas.length - 1, Math.max(0, active + (e.key === "ArrowRight" ? 1 : -1)));
    setActive(next);
  });
}

/* ---------- residue signal map ---------- */
function mapCard(r) {
  const mp = r.map, others = mp.nodes.filter((_, i) => i !== mp.source).map((n) => n.score);
  const lo = Math.min(...others), hi = Math.max(...others);
  const ringKey = (dashed) => {
    const s = svg("svg", { width: 18, height: 18, viewBox: "0 0 18 18" });
    svg("circle", { cx: 9, cy: 9, r: 4, fill: "#8a8a8a" }, s);
    svg("circle", { cx: 9, cy: 9, r: 7.5, fill: "none", stroke: dashed ? "#6b6b6b" : "#161616", "stroke-width": 1.5, "stroke-dasharray": dashed ? "2.5 2" : "none" }, s);
    return s;
  };
  return el("section", { class: "card" },
    el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Residue signal map" }),
        el("p", { class: "card-sub" }, "Signal each residue receives at ", el("span", { class: "mono", text: `γ = ${fmtGamma(mp.gamma)}` }),
          ". Darker means more signal: candidate communication hotspots. Residues are placed at their 3D positions projected onto the protein's two main axes; lines are contacts."))),
    el("div", { class: "map-legend" },
      el("span", { class: "ramp" }, el("span", { class: "mono", text: fmt(lo) }), el("span", { class: "ramp-bar" }), el("span", { class: "mono", text: fmt(hi) })),
      el("span", { class: "ring-key" }, ringKey(false), el("span", {}, "source ", el("span", { class: "mono", text: r.summary.source }))),
      el("span", { class: "ring-key" }, ringKey(true), el("span", {}, "target ", el("span", { class: "mono", text: r.summary.target })))),
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
  const W = Math.max(300, host.clientWidth);
  const pad = 22;
  const inner = W - 2 * pad;
  const H = Math.round(Math.min(inner * Math.max(mp.aspect, 0.35), 560) + 2 * pad);
  const sx = Math.min(inner, (H - 2 * pad) / Math.max(mp.aspect, 1e-3));
  const offX = pad + (inner - sx) / 2;
  const P = mp.nodes.map((n) => [offX + n.x * sx, pad + n.y * sx]);
  const others = mp.nodes.filter((_, i) => i !== mp.source).map((n) => n.score);
  const lo = Math.min(...others), hi = Math.max(...others);
  const rad = Math.max(3, Math.min(7, 150 / Math.sqrt(N)));

  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
    "aria-label": `Contact network of ${N} residues shaded by signal received` }, host);
  const eg = svg("g", {}, root);
  const nbrs = mp.nodes.map(() => []);
  const edgeEls = mp.edges.map(([a, b]) => {
    nbrs[a].push(b); nbrs[b].push(a);
    return svg("line", { x1: P[a][0], y1: P[a][1], x2: P[b][0], y2: P[b][1], stroke: "#e6e6e6", "stroke-width": 1 }, eg);
  });
  const edgesOf = mp.nodes.map(() => []);
  mp.edges.forEach(([a, b], k) => { edgesOf[a].push(k); edgesOf[b].push(k); });

  const ng = svg("g", {}, root);
  // draw low-signal residues first so dark hotspots sit on top
  const order = mp.nodes.map((_, i) => i).sort((a, b) => mp.nodes[a].score - mp.nodes[b].score);
  const nodeEls = [];
  for (const i of order) {
    const t = i === mp.source ? 1 : (mp.nodes[i].score - lo) / (hi - lo || 1);
    nodeEls[i] = svg("circle", { cx: P[i][0], cy: P[i][1], r: rad, fill: shade(t), stroke: "#fff", "stroke-width": 1.5 }, ng);
  }
  const ring = (i, dashed, label) => {
    svg("circle", { cx: P[i][0], cy: P[i][1], r: rad + 4.5, fill: "none", stroke: dashed ? "#6b6b6b" : "#161616",
      "stroke-width": 1.5, "stroke-dasharray": dashed ? "3 2.5" : "none" }, root);
    const right = P[i][0] < W - 90;
    svg("text", { x: P[i][0] + (right ? rad + 9 : -rad - 9), y: P[i][1] + 3.5, "text-anchor": right ? "start" : "end",
      "font-size": 11, fill: "#3a3a3a", "paint-order": "stroke", stroke: "#fff", "stroke-width": 3, text: label }, root);
  };
  ring(mp.target, true, mp.nodes[mp.target].id);
  ring(mp.source, false, mp.nodes[mp.source].id);

  // hover: nearest residue within a generous radius; its contacts darken
  const hit = svg("rect", { x: 0, y: 0, width: W, height: H, fill: "transparent" }, root);
  let active = -1;
  function setActive(i, e) {
    if (active >= 0) {
      nodeEls[active].setAttribute("r", rad); nodeEls[active].setAttribute("stroke", "#fff");
      edgesOf[active].forEach((k) => { edgeEls[k].setAttribute("stroke", "#e6e6e6"); });
    }
    active = i;
    if (i < 0) { hideTip(); return; }
    nodeEls[i].setAttribute("r", rad + 1.5); nodeEls[i].setAttribute("stroke", "#161616");
    edgesOf[i].forEach((k) => { edgeEls[k].setAttribute("stroke", "#9a9a9a"); });
    const n = mp.nodes[i];
    const role = i === mp.source ? "  ·  source" : i === mp.target ? "  ·  target" : "";
    showTip(e.clientX, e.clientY, `${n.id}  ${n.resname}${role}`,
      [{ value: fmt(n.score, 4), label: "signal" }, { value: String(n.degree), label: "contacts" }]);
  }
  hit.addEventListener("pointermove", (e) => {
    const box = root.getBoundingClientRect();
    const x = ((e.clientX - box.left) / box.width) * W, y = ((e.clientY - box.top) / box.height) * H;
    let best = -1, bd = 14 * 14;
    P.forEach((p, i) => { const d = (p[0] - x) ** 2 + (p[1] - y) ** 2; if (d < bd) { bd = d; best = i; } });
    if (best !== active) setActive(best, e);
    else if (best >= 0) setActive(best, e);
  });
  hit.addEventListener("pointerleave", () => setActive(-1));
}

/* ---------- ranked table ---------- */
function rankingCard(r) {
  const top = r.ranking, max = Math.max(...top.map((x) => x.score));
  const rows = top.map((x) => el("tr", {},
    el("td", { class: "rank", text: String(x.rank) }),
    el("td", { text: x.id }),
    el("td", { text: x.resname }),
    el("td", { class: "num" }, el("div", { class: "bar-cell" },
      el("span", { text: fmt(x.score, 4) }),
      el("span", { class: "bar" }, el("span", { style: `width:${(100 * x.score / max).toFixed(1)}%` })))),
    el("td", { class: "num", text: String(x.degree) })));
  return el("details", { class: "card card-fold", open: true },
    el("summary", {}, el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: "Top residues by signal" }),
        el("p", { class: "card-sub" }, `Top ${top.length} of ${r.ranking_total} at `, el("span", { class: "mono", text: `γ = ${fmtGamma(r.map.gamma)}` }), ", source excluded")),
      el("span", { class: "chev", "aria-hidden": "true" }))),
    el("div", { class: "table-scroll" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", { text: "#" }), el("th", { text: "Residue" }), el("th", { text: "Name" }),
        el("th", { class: "num", text: "Signal" }), el("th", { class: "num", text: "Contacts" }))),
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
      item(f.qmod, `${name}.qmod`, "Classiq quantum-walk circuit (noise-free; add dephasing at run time)"),
      item(f.graphml, `${name}.graphml`, "Residue contact network"),
      item(f.hump_csv, `${name}_hump.csv`, "Transport versus noise"),
      item(f.control_csv, `${name}_control_hump.csv`, "Random site-energy control"),
      item(f.parameters, `${name}_parameters.json`, "Every setting used, with a code version marker")),
    el("div", { class: "version", text: ver }),
    el("details", { class: "fold" }, el("summary", { text: "Show parameters" }),
      el("pre", { class: "params", text: JSON.stringify(r.parameters, null, 2) })));
}

/* redraw the SVGs at the new width */
let rs = null;
window.addEventListener("resize", () => { clearTimeout(rs); rs = setTimeout(() => { hideTip(); drawChart(); drawMap(); }, 120); });
