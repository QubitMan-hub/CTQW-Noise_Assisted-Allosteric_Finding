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
/* theme colours, read from the CSS tokens so charts follow light and dark mode */
const C = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
function rgbOf(c) {
  if (c.startsWith("#")) {
    const h = c.length === 4 ? c.slice(1).split("").map((x) => x + x).join("") : c.slice(1, 7);
    return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16));
  }
  return (c.match(/[\d.]+/g) || [0, 0, 0]).slice(0, 3).map(Number);
}
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
  btn.disabled = true; btn.textContent = "Walking";
  $("#results").replaceChildren();
  $("#loading").hidden = false;
  startLoadingFun();
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
    else { renderResults(body); loadRecent(body.run_id); loadTable(); }
  } catch {
    showError("Could not reach the local server. Is web/app.py still running?");
  } finally {
    clearInterval(timer);
    stopLoadingFun();
    $("#loading").hidden = true;
    btn.disabled = false; btn.textContent = "Run the walk";
  }
});

/* ---------- results ---------- */
let current = null;
const charts = [];
function renderResults(r) {
  current = r;
  charts.length = 0;
  const root = $("#results");
  const pair = (...cards) => el("div", { class: "pair" + (cards.length === 2 ? " two" : "") }, cards);
  root.replaceChildren(...[summaryCard(r),
    r.allosteric ? (r.allosteric.adjusted ? pair(allostericCard(r), adjustedCard(r)) : allostericCard(r)) : null,
    pair(transportCard(r), mapCard(r)), rankingCard(r), filesCard(r)].filter(Boolean));
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

/* ---------- plain-language help ---------- */
function fmtP(p) { return p == null ? "–" : p < 0.001 ? "< 0.001" : p.toFixed(p < 0.01 ? 3 : 2); }
function guide(items) {
  return el("details", { class: "fold guide" }, el("summary", { text: "How to read this" }),
    el("ul", { class: "guide-list" }, items.map(([lead, rest]) => el("li", {}, el("strong", { text: lead }), " " + rest))));
}

// the run's answers to the questions people actually ask, from its own numbers
function takeaways(r) {
  const out = [], t = r.transport, a = r.allosteric;
  const mark = (kind) => el("span", { class: "tk-mark " + kind, text: kind === "yes" ? "yes" : kind === "no" ? "no" : "partly" });
  const line = (kind, q, a2) => out.push(el("li", {}, mark(kind), el("div", {}, el("div", { class: "tk-q", text: q }), el("div", { class: "tk-a", text: a2 }))));
  if (a && a.adjusted) {
    const d = a.adjusted, st = d.stats, best = st.peak_value, sig = d.significance, pv = sig ? sig.p_value : null;
    const signif = pv !== null && pv < 0.05;
    line(best >= 0.6 && (pv === null || signif) ? "yes" : best > 0.55 ? "partly" : "no", "Does the walk point at the known allosteric site?",
      best >= 0.6 && pv !== null && !signif ? `Possibly. It scores ${fmt(best)} beyond plain closeness (0.5), but with ${a.n_allosteric} known allosteric residues that could still be luck (p = ${fmtP(pv)}).`
      : best >= 0.6 ? `Yes. Among residues equally far from the active site, it ranks the allosteric ones higher (best ${fmt(best)}, where 0.5 is what closeness alone gives${pv !== null ? `; p = ${fmtP(pv)}` : ""}).`
      : best > 0.55 ? `Only a little. Beyond plain closeness it scores ${fmt(best)} (0.5 is closeness alone).`
      : `No. Once distance from the active site is taken out it scores ${fmt(best)}, about chance. Whatever the raw score shows here is just closeness.`);
    if (d.control) {
      const seeds = d.control.best_seeds || [], beat = seeds.filter((x) => best > x).length;
      line(best >= 0.55 && beat === seeds.length ? "yes" : best >= 0.55 && beat > 0 ? "partly" : "no", "Is it about this protein's chemistry?",
        best < 0.55 ? "Nothing to explain: there is no signal beyond distance."
        : beat === seeds.length ? `Yes. Its best beats the best of all ${seeds.length} random-energy runs, so the hydropathy pattern matters.`
        : beat > 0 ? `Partly. It beats ${beat} of ${seeds.length} random-energy runs, so much of it comes from the contact network itself.`
        : `No. Random energies do as well, so the contact network alone explains it.`);
    }
    const cl = a.classical;
    if (cl) {
      const cb = Math.max(...cl.auc_adjusted.filter((x) => x !== null)), diff = best - cb;
      line(diff >= 0.02 && best > 0.55 ? "yes" : diff > -0.02 && best > 0.55 ? "partly" : "no", "Does it need to be quantum?",
        best <= 0.55 ? `Nothing to explain here. A classical random walk scores ${fmt(cb)}.`
        : diff >= 0.02 ? `It helps. The quantum walk's ${fmt(best)} beats the best classical random walk on the same contacts (${fmt(cb)}, over ${cl.rates.length} hopping speeds).`
        : diff > -0.02 ? `Not clearly. A classical random walk on the same contacts does about as well (${fmt(cb)}).`
        : `No. A classical random walk on the same contacts does better (${fmt(cb)}).`);
    }
    const helps = st.verdict === "hump" && st.gain_abs >= 0.01 && best > 0.55;
    line(helps ? "yes" : "no", "Does a little noise help find it?",
      helps ? `Yes. The best is at γ = ${fmtGamma(st.peak_gamma)}, ${fmt(st.gain_abs)} above the better of no noise (${fmt(st.quantum_value)}) and heavy noise (${fmt(st.dephased_value)}).`
      : st.verdict === "quantum" ? "No. It does best with no noise at all." : best <= 0.55 ? "Not meaningfully: there is little signal to improve."
      : `Not meaningfully: the peak is less than 0.01 above the ends.`);
  }
  const ts = t.stats;
  if (ts.verdict === "hump") {
    const c = t.control;
    line(c && c.exceeds_seeds < c.n_seeds ? "partly" : "yes", "Does noise help the signal travel?",
      `Yes, it peaks at γ = ${fmtGamma(ts.peak_gamma)} (${pct(ts.gain_rel)} over the better end).` +
      (c ? (c.exceeds_seeds < c.n_seeds ? ` But random energies show it too (${c.seed_verdicts.filter((v) => v === "hump").length} of ${c.n_seeds}), so this alone says little about this protein.` : " And it is stronger than in every random-energy run.")
         : " Without the control you cannot tell whether that is special: most networks do this."));
  } else {
    line("no", "Does noise help the signal travel?", ts.verdict === "quantum" ? "No, it travels best with no noise." : "No, more noise always did better here.");
  }
  if (!a) out.push(el("li", { class: "tk-hint" }, "To test site finding, pick a protein with a known allosteric site (the key button above) or give your own active-site and allosteric residues."));
  return el("div", { class: "takeaways" }, el("div", { class: "tk-title", text: "What this run says" }), el("ul", {}, out));
}

function summaryCard(r) {
  const s = r.summary, t = r.transport.stats, a = r.allosteric;
  const pill = (st, what) => el("span", { class: "pill " + (st.verdict === "hump" ? "solid" : "outline"), text: verdictText(st, what) });
  const stat = (label, value, hint) => el("div", { class: "stat", title: hint || "" },
    el("div", { class: "stat-label", text: label }), el("div", { class: "stat-value", text: value }));
  const notes = r.notes.length ? el("ul", { class: "notes" }, r.notes.map((n) => el("li", { text: n }))) : null;
  return el("section", { class: "card" },
    el("div", { class: "summary-top" },
      el("div", {},
        el("h2", { class: "run-name", text: r.name }),
        el("div", { class: "run-meta", text: `chains ${s.chains.join(", ")}  ·  ran in ${r.elapsed_s} s` })),
      el("div", { class: "pills" }, a ? (a.adjusted ? pill(a.adjusted.stats, "beyond distance") : pill(a.stats, "allosteric AUC")) : null, pill(t, "transport"))),
    el("div", { class: "stats" },
      stat("Residues", fmtInt(s.residues), "Amino acids in the network (one node each)"),
      stat("Contacts", fmtInt(s.contacts), "Residue pairs closer than 8 Å; the walker hops along these"),
      a && a.adjusted && a.adjusted.significance ? stat("p-value", fmtP(a.adjusted.significance.p_value), "Chance of scoring this well by luck: the allosteric labels were shuffled within each distance group 10,000 times. Below 0.05 is the usual bar.")
        : a ? stat("Best AUC", fmt(a.stats.peak_value), "How well the raw signal ranks the known allosteric residues: 0.5 is a coin flip, 1 is perfect")
        : stat("Transport peak γ", fmtGamma(t.peak_gamma), "The noise level at which most signal reaches the far side"),
      a && a.adjusted ? stat("Best AUC, same distance", fmt(a.adjusted.stats.peak_value), "The same, but only comparing residues equally far from the active site: 0.5 means closeness explains everything")
        : stat("Gain over ends", t.verdict === "hump" ? pct(t.gain_rel) : "none", "How much better the best noise level is than both no noise and heavy noise")),
    takeaways(r),
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
  if (!ctrl) return el("p", { class: "note", text: "Ran without the random-energy control. Humps like this show up in almost any network, so on its own it is not evidence about this protein." });
  const humps = ctrl.seed_verdicts ? `${ctrl.seed_verdicts.filter((v) => v === "hump").length} of ${ctrl.seed_verdicts.length} seeds show a hump; ` : "";
  return el("p", { class: "note" }, "Random-energy control: ",
    el("span", { class: "mono", text: verdictText(ctrl.stats, "mean curve").toLowerCase() }), `. ${humps}`,
    `${humps ? "the" : "The"} real curve's gain exceeds ${ctrl.exceeds_seeds} of ${ctrl.n_seeds} seed gains`,
    ctrl.z !== null ? ` (z = ${fmt(ctrl.z, 2)}).` : ".",
    st.verdict === "hump" && ctrl.exceeds_seeds < ctrl.n_seeds
      ? " A hump that random energies reproduce is not specific to the hydropathy model." : "");
}

function aucCard(r, blk, o) {
  const a = r.allosteric, st = blk.stats;
  const series = [{ values: blk.auc, label: `${r.summary.site_energy} site energies`, style: "main" }];
  if (blk.control) series.push({ values: blk.control.mean, sd: blk.control.sd, label: "random energies (mean ± sd)", style: "control" });
  const chart = { id: o.id, gammas: r.gammas, series, yLabel: o.yLabel, peakGamma: st.peak_gamma,
    baselines: Object.entries(blk.baselines).map(([k, v]) => ({ value: v, label: k })), chance: 0.5, auc: true, digits: 3 };
  charts.push(chart);
  return el("section", { class: "card" },
    el("div", { class: "card-head" },
      el("div", {},
        el("h3", { class: "card-title", text: o.title }),
        el("p", { class: "card-sub", text: o.sub }))),
    legend(series, blk.baselines),
    el("div", { class: "chart-wrap", id: chart.id }),
    statsRow(st, "auc"),
    el("div", { class: "metric-row" }, Object.entries(blk.baselines).map(([k, v]) => el("span", {}, `baseline: ${k}`, el("strong", { text: fmt(v) })))),
    o.extra || null,
    controlNote(blk.control, st),
    o.guide ? guide(o.guide) : null,
    el("div", { class: "btn-row" },
      download(r.files[o.tag + "_png"], "Download figure", "PNG 300 dpi", true),
      download(r.files[o.tag + "_svg"], "Figure", "SVG"),
      download(r.files[o.tag + "_csv"], "Data", "CSV")),
    dataTable(r.gammas, series, 3));
}

function allostericCard(r) {
  const a = r.allosteric, miss = [...a.missing.active, ...a.missing.allosteric];
  return aucCard(r, a, { id: "chart-allo", tag: "allosteric", yLabel: "ROC AUC, allosteric residues",
    title: "Allosteric recovery versus noise",
    guide: [
      ["The solid line:", "for each noise level, how well the signal a residue receives picks out the known allosteric residues. 0.5 is a coin flip, 1.0 would put every allosteric residue first."],
      ["The dotted lines:", "simple guesses to beat. \u2018Proximity\u2019 ranks residues by how few contacts they are from the active site; \u2018degree\u2019 by how many contacts they have; the classical walk is an ordinary random walk on the same contacts."],
      ["The grey band:", "the same test with random site energies (mean \u00b1 sd over 5 runs). A solid line inside the band is not about this protein's chemistry."],
      ["Careful:", "signal fades with distance, so this raw score mostly measures closeness. The \u2018Beyond distance\u2019 chart is the fairer test."],
    ],
    sub: `Walk from the ${a.n_active} active-site residues; how well the signal each residue receives ranks the ${a.n_allosteric} known allosteric residues among ${a.n_candidates} candidates (ROC AUC; 0.5 is chance).`,
    extra: miss.length || a.name_mismatch.length ? el("p", { class: "note", text:
      `${miss.length ? `Not in the network: ${miss.join(", ")}. ` : ""}${a.name_mismatch.length ? `Name mismatch: ${a.name_mismatch.join("; ")}.` : ""}` }) : null });
}

function adjustedCard(r) {
  const d = r.allosteric.adjusted;
  return aucCard(r, d, { id: "chart-adj", tag: "adjusted", yLabel: "ROC AUC, same distance",
    title: "Beyond distance",
    extra: d.significance ? el("p", { class: "note" }, "Significance: ", el("strong", { text: `p = ${fmtP(d.significance.p_value)}` }),
      `. The allosteric labels were shuffled within each distance group ${d.significance.n_perm.toLocaleString()} times, and each time the best score over all noise levels was taken, so picking the best γ is accounted for (shuffles reach ${fmt(d.significance.null_95)} one time in twenty). Classical random walk: best ${fmt(r.allosteric.classical.significance_adjusted.observed_best)}, p = ${fmtP(r.allosteric.classical.significance_adjusted.p_value)}.`) : null,
    guide: [
      ["The idea:", "residues are grouped by how many contacts they are from the active site, and each is only ranked against its own group. Being close no longer helps."],
      ["0.5 line:", "what closeness alone scores here. Above it, the walk sees something extra about where the allosteric site is; below it, the walk points away from it."],
      ["Grey band:", "random site energies. Solid line above the band means the protein's hydropathy pattern is doing the work, not just the shape of the network."],
      ["A hump:", "a peak between the two ends means a moderate amount of noise helps the walk find the site. Flat or highest at an end means it does not."],
      ["Classical walk line:", "an ordinary random walk on the same contacts, at its best hopping speed. Beating it is the case that the quantum part matters."],
      ["p-value:", "how often shuffled labels score as well. With only a handful of allosteric residues, even 0.7 can happen by luck, so check it."],
    ],
    sub: `The same test, but each residue is only compared with residues equally far from the active site (${d.shells} distance shells). Plain closeness scores 0.5 here, so anything above it is what the walk adds.` });
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
    guide([
      ["The solid line:", "how much of the walker's signal reaches the far part of the protein (residues at least half the network's width away), for each noise level."],
      ["Left edge, γ = 0:", "a perfectly quantum walk. Waves interfere and can get stuck near the start."],
      ["Right edge, large γ:", "so much noise that the walker keeps getting \u2018measured\u2019 and barely moves (the quantum Zeno effect)."],
      ["A hump in between:", "a little noise breaks up the interference and helps the signal spread. This happens in almost any disordered network, which is why the grey random-energy band matters."],
    ]),
    el("div", { class: "btn-row" },
      download(r.files.hump_png, "Download figure", "PNG 300 dpi", !r.allosteric),
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
  const m = { l: 64, r: 16, t: 24, b: 44 };
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
    svg("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: Math.abs(v - axis.min) < 1e-12 ? C("--line-2") : C("--grid"), "stroke-width": 1, "shape-rendering": "crispEdges" }, root);
    svg("text", { x: m.l - 10, y: y + 3.5, "text-anchor": "end", "font-size": 11, fill: C("--muted"), text: v.toFixed(axis.decimals) }, root);
  }
  for (const t of [0, 0.1, 1, 10, 100].filter((t) => t <= gmax * 1.01)) {
    svg("line", { x1: X(t), x2: X(t), y1: m.t + ih, y2: m.t + ih + 4, stroke: C("--line-2"), "shape-rendering": "crispEdges" }, root);
    svg("text", { x: X(t), y: m.t + ih + 18, "text-anchor": "middle", "font-size": 11, fill: C("--muted"), text: fmtGamma(t) }, root);
  }
  svg("text", { x: m.l + iw / 2, y: H - 6, "text-anchor": "middle", "font-size": 11.5, fill: C("--ink-2"), text: "dephasing rate γ  (log scale)" }, root);
  svg("text", { x: 12, y: m.t + ih / 2, "text-anchor": "middle", "font-size": 11.5, fill: C("--ink-2"),
    transform: `rotate(-90 12 ${m.t + ih / 2})`, text: c.yLabel }, root);
  if (c.chance !== null) svg("line", { x1: m.l, x2: W - m.r, y1: Y(c.chance), y2: Y(c.chance), stroke: C("--line-2"), "stroke-width": 1 }, root);
  const short = { "proximity to active site": "proximity", "contact degree": "degree", "classical walk, best rate": "classical walk",
    "contact degree, same distance": "degree, same distance" };
  let lastY = -Infinity;                  // top to bottom, nudging labels apart so close baselines stay readable
  for (const b of [...c.baselines].sort((p, q) => q.value - p.value)) {
    svg("line", { x1: m.l, x2: W - m.r, y1: Y(b.value), y2: Y(b.value), stroke: C("--ctrl"), "stroke-width": 1, "stroke-dasharray": "1.5 3" }, root);
    const ly = Math.max(Y(b.value) - 5, lastY + 13);
    lastY = ly;
    svg("text", { x: W - m.r - 4, y: ly, "text-anchor": "end", "font-size": 10.5, fill: C("--muted"),
      "paint-order": "stroke", stroke: C("--card"), "stroke-width": 3,
      text: `${short[b.label] || b.label} ${fmt(b.value)}` }, root);
  }
  const pg = c.peakGamma;
  svg("line", { x1: X(pg), x2: X(pg), y1: m.t - 6, y2: m.t + ih, stroke: C("--faint"), "stroke-width": 1, "stroke-dasharray": "3 3" }, root);
  svg("text", { x: X(pg) + (X(pg) > W - m.r - 110 ? -6 : 6), y: m.t - 10, "text-anchor": X(pg) > W - m.r - 110 ? "end" : "start",
    "font-size": 11, fill: C("--muted"), text: `best γ = ${fmtGamma(pg)}` }, root);

  const main = c.series[0];
  const pts = (vals) => gammas.map((g, i) => [X(g), Y(vals[i])]);
  if (!c.auc) {
    const mp = pts(main.values);
    svg("path", { d: `M${mp[0][0]},${Y(axis.min)} ` + mp.map((p) => `L${p[0]},${p[1]}`).join(" ") + ` L${mp[mp.length - 1][0]},${Y(axis.min)} Z`,
      fill: C("--ink"), "fill-opacity": 0.05 }, root);
  }
  const dots = [];
  [...c.series].reverse().forEach((s) => {
    const isMain = s === main;
    if (s.sd) {
      const up = gammas.map((g, i) => [X(g), Y(s.values[i] + s.sd[i])]), dn = gammas.map((g, i) => [X(g), Y(s.values[i] - s.sd[i])]).reverse();
      svg("path", { d: "M" + up.concat(dn).map((q) => q.join(",")).join(" L") + " Z", fill: C("--ctrl"), "fill-opacity": 0.16 }, root);
    }
    const p = pts(s.values);
    svg("path", { d: "M" + p.map((q) => q.join(",")).join(" L"), fill: "none", stroke: isMain ? C("--ink") : C("--ctrl"),
      "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round", "stroke-dasharray": isMain ? "none" : "5 4" }, root);
    if (isMain) dots.push(p.map((q) => svg("circle", { cx: q[0], cy: q[1], r: 3.2, fill: C("--ink"), stroke: C("--card"), "stroke-width": 1.5 }, root)));
  });

  const cross = svg("line", { y1: m.t, y2: m.t + ih, stroke: C("--ink"), "stroke-opacity": 0.25, "stroke-width": 1, visibility: "hidden" }, root);
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
    svg("circle", { cx: 9, cy: 9, r: 4, fill: C("--ctrl") }, s);
    svg("circle", { cx: 9, cy: 9, r: 7.5, fill: "none", stroke: kind === "source" ? C("--ink") : C("--muted"), "stroke-width": 1.5,
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
          `${r.allosteric ? " (the best-AUC rate)" : " (the transport peak)"}. Stronger shades mean more signal. Residues sit at their 3D positions projected onto the protein's two main axes; lines are contacts.`))),
    el("div", { class: "map-legend" },
      el("span", { class: "ramp" }, el("span", { class: "mono", text: fmt(lo) }), el("span", { class: "ramp-bar" }), el("span", { class: "mono", text: fmt(hi) })),
      el("span", { class: "ring-key" }, ringKey("source"), el("span", { text: fromActive ? "active site" : `start ${r.summary.walk_from}` })),
      r.allosteric ? el("span", { class: "ring-key" }, ringKey("other"), el("span", { text: "known allosteric" })) : null),
    el("div", { class: "map-wrap", id: "map" }),
    guide([
      ["Each dot:", "one residue, placed at its real 3D position flattened onto the protein's two longest axes. Lines are contacts."],
      ["Shade:", "how much signal that residue collected over the whole walk. Hover a dot for its name and value."],
      ["Rings:", `solid rings mark where the walk starts${r.allosteric ? "; dashed rings mark the known allosteric residues. Dashed rings in strong shades, far from the start, are what a good prediction looks like" : ""}.`],
    ]));
}

function shade(t) {                     // theme ramp: little signal -> lots of signal
  const lo = rgbOf(C("--ramp-lo")), hi = rgbOf(C("--ramp-hi")), k = Math.min(1, Math.max(0, t));
  return `rgb(${lo.map((v, i) => Math.round(v + (hi[i] - v) * k)).join(",")})`;
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
  const edgeEls = mp.edges.map(([a, b]) => svg("line", { x1: P[a][0], y1: P[a][1], x2: P[b][0], y2: P[b][1], stroke: C("--edge"), "stroke-width": 1 }, eg));
  const edgesOf = mp.nodes.map(() => []);
  mp.edges.forEach(([a, b], k) => { edgesOf[a].push(k); edgesOf[b].push(k); });

  const ng = svg("g", {}, root);
  const order = mp.nodes.map((_, i) => i).sort((a, b) => mp.nodes[a].score - mp.nodes[b].score);
  const nodeEls = [];
  for (const i of order) {
    const n = mp.nodes[i];
    const t = n.role === "source" ? 1 : (n.score - lo) / (hi - lo || 1);
    nodeEls[i] = svg("circle", { cx: P[i][0], cy: P[i][1], r: rad, fill: shade(t), stroke: C("--card"), "stroke-width": 1.5 }, ng);
  }
  const rg = svg("g", {}, root);
  const ring = (i, dashed) => svg("circle", { cx: P[i][0], cy: P[i][1], r: rad + 4, fill: "none", stroke: dashed ? C("--muted") : C("--ink"),
    "stroke-width": 1.4, "stroke-dasharray": dashed ? "3 2.5" : "none" }, rg);
  mp.nodes.forEach((n, i) => { if (n.role === "source") ring(i, false); else if (n.role === "allosteric") ring(i, true); });

  const hit = svg("rect", { x: 0, y: 0, width: W, height: H, fill: "transparent" }, root);
  let active = -1;
  function setActive(i, e) {
    if (active >= 0) {
      nodeEls[active].setAttribute("r", rad); nodeEls[active].setAttribute("stroke", C("--card"));
      edgesOf[active].forEach((k) => edgeEls[k].setAttribute("stroke", C("--edge")));
    }
    active = i;
    if (i < 0) { hideTip(); return; }
    nodeEls[i].setAttribute("r", rad + 1.5); nodeEls[i].setAttribute("stroke", C("--ink"));
    edgesOf[i].forEach((k) => edgeEls[k].setAttribute("stroke", C("--muted")));
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
    el("div", { class: "btn-row" }, download(r.files.ranking_csv, `Full ranking, ${r.ranking_total} residues`, "CSV")),
    guide([
      ["Signal:", "the total time the walker spent on that residue at the noise level shown. The bar compares it with the top residue."],
      ["Contacts:", "how many neighbours the residue has. Well-connected residues tend to collect more signal."],
      ["Read with care:", "the top of this list is mostly residues right next to the start, because signal fades with distance. It is a list of where the signal goes, not a finished prediction."],
    ]));
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
      el("pre", { class: "params", text: JSON.stringify(r.parameters, (k, v) => (k === "citation" ? undefined : v), 2) })));
}

/* redraw the SVGs at the new width */
let rs = null;
window.addEventListener("resize", () => { clearTimeout(rs); rs = setTimeout(() => { hideTip(); redraw(); }, 120); });

/* ======================================================================
   QubitMan extras: pixel art, the dice, the id preview, the loading fun
   ====================================================================== */
PX.wordmark($("#wordmark"));
$("#logo").innerHTML = LOGO.svg(104);           // trusted, built from constants in logo.js
$("#loading-logo").innerHTML = LOGO.svg(40);
PX.icons();

/* ---------- light / dark ---------- */
function isDark() {
  const t = document.documentElement.dataset.theme;
  return t ? t === "dark" : window.matchMedia("(prefers-color-scheme: dark)").matches;
}
function paintToggle() {
  const dark = isDark();
  PX.drawBitmap($("#theme-icon"), dark ? PX.ICONS.sun : PX.ICONS.moon, 3);
  $("#theme-label").textContent = dark ? "light" : "dark";
}
$("#theme-toggle").addEventListener("click", () => {
  const next = isDark() ? "light" : "dark";
  document.documentElement.dataset.theme = next;
  try { localStorage.setItem("qm-theme", next); } catch (e) { /* private mode: fine, just not remembered */ }
  paintToggle();
  if (current) { hideTip(); redraw(); }
});
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
  if (!document.documentElement.dataset.theme) { paintToggle(); if (current) redraw(); }
});
paintToggle();
PX.background($("#bg"));
const walkerBar = PX.walkerBar($("#walker"));

/* ---------- loading: hopping walker + facts ---------- */
const FACTS = [
  "A quantum walker is spread over many residues at once, not hopping between them one at a time.",
  "Dephasing noise scrambles the walker's phases. Too little and interference can trap it; too much and it freezes in place (the quantum Zeno effect).",
  "Noise-assisted transport was first proposed to explain how energy moves efficiently through photosynthetic light-harvesting complexes.",
  "Allosteric sites can sit far from the active site, yet a molecule binding there still changes what the protein does.",
  "Each dot in the network is one amino acid; two are linked when their β-carbons are within 8 Å of each other.",
  "The Protein Data Bank holds over 200,000 experimentally determined structures of proteins and nucleic acids.",
  "With no noise, the walk is solved exactly by diagonalising the network's Hamiltonian; with noise, it is stepped through time.",
  "Every noise level runs as its own job, so more CPU cores means a faster walk.",
];
let factTimer = null, factIdx = Math.floor(Math.random() * FACTS.length);
function showFact() {
  const f = $("#fact");
  f.style.opacity = "0";
  setTimeout(() => { f.textContent = FACTS[factIdx % FACTS.length]; factIdx += 1; f.style.opacity = "1"; }, 250);
}
function startLoadingFun() {
  walkerBar.start();
  showFact(); clearInterval(factTimer); factTimer = setInterval(showFact, 6000);
}
function stopLoadingFun() { walkerBar.stop(); clearInterval(factTimer); }

/* ---------- protein info panel ---------- */
function infoPanel(info, { compact = false, onRoll = null } = {}) {
  const multi = info.chains.length > 1;
  const top = el("div", { class: "info-top" },
    el("span", { class: "info-id", text: info.id }),
    info.classification ? el("span", { class: "pill-tag", text: info.classification }) : null,
    info.year ? el("span", { class: "pill-tag", text: String(info.year) }) : null,
    info.known_site ? el("span", { class: "pill-tag dark", text: "known allosteric site" }) : null);
  const method = [info.method, info.resolution ? `${info.resolution.toFixed(2)} Å` : null].filter(Boolean).join(" · ");
  const mols = info.molecules.map((m) => el("div", {}, el("span", { text: m.name }),
    el("span", { class: "muted", text: `  ·  chain${m.chains.length > 1 ? "s" : ""} ${m.chains.join(", ")}  ·  ${m.residues} residues` })));
  const card = el("div", { class: "info" + (compact ? " compact" : "") }, top,
    info.title ? el("p", { class: "info-title", text: info.title }) : null);

  if (compact) {
    const bits = [info.organism, method, `${info.chains.length} protein chain${info.chains.length === 1 ? "" : "s"}`].filter(Boolean);
    card.append(el("p", { class: "muted", style: "margin:8px 0 0;font-size:13px", text: bits.join("  ·  ") }));
    if (info.known_site) {
      card.append(el("p", { style: "margin:8px 0 0;font-size:13px", text:
        `Known allosteric site (ligand ${info.known_site.ligand}): the walk starts at the active site and gets scored automatically.` }));
    } else if (multi && info.suggested_chains) {
      const btn = el("button", { type: "button", class: "linkish", text: `use chain ${info.suggested_chains} only` });
      const msg = el("span", { class: "muted" });
      btn.addEventListener("click", () => { $("#chains").value = info.suggested_chains; msg.textContent = `  ·  chains set to ${info.suggested_chains}`; });
      card.append(el("p", { style: "margin:8px 0 0;font-size:13px" },
        `This file holds ${info.chains.length} chains. Tip: `, btn, msg));
    }
    return card;
  }

  const grid = el("dl", { class: "info-grid" },
    info.organism ? [el("dt", { text: "Organism" }), el("dd", { text: info.organism })] : null,
    method ? [el("dt", { text: "Method" }), el("dd", { text: method })] : null,
    [el("dt", { text: "Protein" }), el("dd", {}, mols)],
    info.ligands.length ? [el("dt", { text: "Ligands" }), el("dd", {},
      info.ligands.map((l) => el("span", { class: "lig", title: l.name || l.id, text: l.id })))] : null,
    info.known_site ? [el("dt", { text: "Known site" }), el("dd", {
      text: `${info.known_site.n_allosteric} allosteric residues around ligand ${info.known_site.ligand}. The walk will start at the active site and be scored against them.` })] : null);
  card.append(grid);

  const run = el("button", { type: "button", class: "btn-primary small", text: "Run this one" });
  run.addEventListener("click", () => {
    setFile(null);
    pdbInput.value = info.id;
    $("#chains").value = info.known_site ? "" : (info.suggested_chains || "");
    form.querySelector('input[name="labels"][value="auto"]').checked = true;
    $("#custom-labels").hidden = true;
    showPreview(info);
    form.requestSubmit();
  });
  // estimate_s is for the plain walk; the random-energy control adds 5 more sweeps
  const secs = info.estimate_s && info.estimate_s * ($("#control").checked ? CONTROL_FACTOR : 1);
  const est = secs ? `roughly ${secs < 90 ? Math.round(secs) + " s" : Math.round(secs / 60) + " min"} on 4 cores` : "";
  const chainNote = !info.known_site && multi && info.suggested_chains ? `chain ${info.suggested_chains} (${info.run_residues} residues)` : `${info.run_residues} residues`;
  card.append(el("div", { class: "info-actions" }, run,
    onRoll ? (() => { const b = el("button", { type: "button", class: "linkish", text: "roll again" }); b.addEventListener("click", onRoll); return b; })() : null,
    el("span", { class: "muted", style: "font-size:13px", text: [chainNote, est].filter(Boolean).join("  ·  ") })));
  return card;
}

const CONTROL_FACTOR = 4;   // measured: 1IWH 19 s plain, 71 s with the control

/* ---------- the dice ---------- */
async function roll(kind, button) {
  const host = $("#dice-result");
  document.querySelectorAll(".btn-roll").forEach((b) => { b.disabled = true; });
  button.classList.add("rolling", "pressed");
  host.replaceChildren(el("div", { class: "info" },
    el("span", { class: "skeleton-line", style: "width:30%" }), el("span", { class: "skeleton-line", style: "width:85%" }),
    el("span", { class: "skeleton-line", style: "width:60%" })));
  try {
    const res = await fetch(`/api/random${kind === "allosteric" ? "?kind=allosteric" : ""}`);
    const info = await res.json();
    if (!res.ok || info.error) throw new Error(info.error || "The dice fell off the table. Roll again.");
    host.replaceChildren(infoPanel(info, { onRoll: () => roll(kind, button) }));
  } catch (err) {
    host.replaceChildren(el("p", { class: "error", text: err.message || "Could not reach the server." }));
  } finally {
    document.querySelectorAll(".btn-roll").forEach((b) => { b.disabled = false; });
    button.classList.remove("rolling", "pressed");
  }
}
document.querySelectorAll(".btn-roll").forEach((b) => b.addEventListener("click", () => roll(b.dataset.kind, b)));

/* ---------- live preview while typing an id ---------- */
let previewToken = 0, previewTimer = null;
function showPreview(info) {
  const hit = tableMatch(info.id || pdbInput.value.trim().toUpperCase());
  let note = null;
  if (hit) {
    const open = el("button", { type: "button", class: "linkish", text: "open the stored result" });
    open.addEventListener("click", () => reopen(hit.run_id));
    note = el("div", { class: "in-table" }, el("span", { class: "chip-mini solid", text: "in the table" }),
      el("span", { text: `walked ${ago(hit.finished_utc)}${hit.beyond != null ? `, beyond distance ${fmt(hit.beyond)}` : ""}.` }), open,
      el("span", { class: "muted", text: "or run it again with new settings." }));
  }
  $("#id-preview").replaceChildren(infoPanel(info, { compact: true }), note || "");
}
pdbInput.addEventListener("input", () => {
  clearTimeout(previewTimer);
  const id = pdbInput.value;
  const token = ++previewToken;
  if (!/^[0-9][A-Z0-9]{3}$/.test(id)) { $("#id-preview").replaceChildren(); return; }
  previewTimer = setTimeout(async () => {
    $("#id-preview").replaceChildren(el("div", { class: "info compact" },
      el("span", { class: "skeleton-line", style: "width:40%;margin-top:0" }), el("span", { class: "skeleton-line", style: "width:75%" })));
    try {
      const res = await fetch(`/api/info/${id}`);
      const info = await res.json();
      if (token !== previewToken) return;
      if (!res.ok || info.error) { $("#id-preview").replaceChildren(el("p", { class: "muted", style: "margin:10px 0 0;font-size:13px", text: info.error })); return; }
      showPreview(info);
    } catch {
      if (token === previewToken) $("#id-preview").replaceChildren();
    }
  }, 350);
});

/* ---------- sidebar: recent walks + facts ---------- */
function ago(iso) {
  if (!iso) return "";
  const s = Math.max(0, (Date.now() - Date.parse(iso)) / 1000);
  if (s < 60) return "just now";
  if (s < 3600) return `${Math.floor(s / 60)} min ago`;
  if (s < 86400) return `${Math.floor(s / 3600)} h ago`;
  return `${Math.floor(s / 86400)} d ago`;
}
function verdictShort(v, g) { return v === "hump" ? `hump at γ ${fmtGamma(g)}` : v === "quantum" ? "best at γ 0" : "best most dephased"; }
let activeRun = null;
async function loadRecent(markId) {
  if (markId) activeRun = markId;
  let items = [];
  try { items = await (await fetch("/api/runs")).json(); } catch { return; }
  const list = $("#recent");
  if (!Array.isArray(items) || !items.length) {
    list.replaceChildren(el("li", { class: "recent-empty", text: "Your walks show up here. Click one to reopen its results instantly." }));
    return;
  }
  list.replaceChildren(...items.map((it) => {
    const main = it.auc ? `AUC ${fmt(it.auc.best)} · ${verdictShort(it.transport.verdict, it.transport.peak_gamma)}`
      : verdictShort(it.transport.verdict, it.transport.peak_gamma);
    const btn = el("button", { type: "button", class: "recent-item" + (it.run_id === activeRun ? " active" : ""),
      title: `Reopen ${it.name}` },
      el("span", { class: "recent-id", text: it.name.slice(0, 8) }),
      el("span", { class: "recent-main", text: main }),
      el("span", { class: "recent-sub", text: `${it.residues} res · ${it.control ? "with control · " : ""}${ago(it.finished_utc)}` }));
    btn.addEventListener("click", () => reopen(it.run_id));
    return el("li", {}, btn);
  }));
}
async function reopen(id) {
  clearError();
  try {
    const res = await fetch(`/api/runs/${id}`);
    const body = await res.json();
    if (!res.ok || body.error) { showError(body.error || "Could not reopen that walk."); loadRecent(); return; }
    activeRun = id;
    renderResults(body);
    document.querySelectorAll(".recent-item").forEach((b) => b.classList.remove("active"));
    loadRecent();
  } catch { showError("Could not reach the local server. Is web/app.py still running?"); }
}
let sideIdx = Math.floor(Math.random() * FACTS.length);
function sideFact() {
  const f = $("#side-fact");
  f.style.opacity = "0";
  setTimeout(() => { f.textContent = FACTS[sideIdx % FACTS.length]; sideIdx += 1; f.style.opacity = "1"; }, 200);
}
$("#next-fact").addEventListener("click", sideFact);
sideFact();
setInterval(() => { if (!document.hidden) sideFact(); }, 12000);
loadRecent();

/* ---------- Walk | Table ---------- */
let tableRows = [], tableSort = { key: "finished_utc", dir: -1 }, openRow = null;
function setView(view, { scroll = true } = {}) {
  const table = view === "table";
  $("#view-walk").hidden = table;
  $("#view-table").hidden = !table;
  document.querySelectorAll(".view-tab").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.view === view)));
  try { localStorage.setItem("qm-view", view); } catch {}
  if (table) loadTable();
  else redraw();                          // charts measure their box, so draw them once visible
  if (scroll) $(".view-switch").scrollIntoView({ behavior: "smooth", block: "start" });
}
document.querySelectorAll(".view-tab").forEach((b) => b.addEventListener("click", () => setView(b.dataset.view)));

const TABLE_COLS = [
  { key: "name", label: "Protein" },
  { key: "residues", label: "Residues", num: true },
  { key: "beyond", label: "Beyond distance", num: true, hint: "0.5 = closeness alone" },
  { key: "p_value", label: "p", num: true },
  { key: "classical", label: "Classical", num: true },
  { key: "beyond_gamma", label: "Noise helps" },
  { key: "transport_verdict", label: "Transport" },
  { key: "control", label: "Control" },
  { key: "finished_utc", label: "When" },
];

async function loadTable() {
  try { tableRows = await (await fetch("/api/table")).json(); } catch { return; }
  if (!Array.isArray(tableRows)) tableRows = [];
  $("#table-count").textContent = tableRows.length ? String(tableRows.length) : "";
  if (!$("#view-table").hidden) renderTable();
}

function sortValue(r, key) {
  const v = r[key];
  if (key === "beyond_gamma") return r.beyond_verdict === "hump" ? v : null;
  if (typeof v === "boolean") return v ? 1 : 0;
  return v === undefined ? null : v;
}

function renderTable() {
  const q = $("#table-filter").value.trim().toLowerCase(), onlyLab = $("#table-labelled").checked;
  let rows = tableRows.filter((r) => (!onlyLab || r.labelled) &&
    (!q || r.name.toLowerCase().includes(q) || (r.title || "").toLowerCase().includes(q)));
  const { key, dir } = tableSort;
  rows = rows.slice().sort((a, b) => {
    const x = sortValue(a, key), y = sortValue(b, key);
    if (x === null && y === null) return 0;
    if (x === null) return 1;                     // blanks always last
    if (y === null) return -1;
    return (x < y ? -1 : x > y ? 1 : 0) * dir;
  });
  const lab = tableRows.filter((r) => r.beyond != null);
  const withP = lab.filter((r) => r.p_value != null), sig = withP.filter((r) => r.p_value < 0.05);
  const withC = lab.filter((r) => r.classical != null), beat = withC.filter((r) => r.beyond > r.classical);
  $("#table-summary").replaceChildren(...[
    el("span", {}, "walks", el("strong", { text: String(tableRows.length) })),
    el("span", {}, "proteins", el("strong", { text: String(new Set(tableRows.map((r) => r.name)).size) })),
    lab.length ? el("span", {}, "with a known site", el("strong", { text: String(lab.length) })) : null,
    withP.length ? el("span", {}, "p < 0.05", el("strong", { text: `${sig.length} of ${withP.length}` })) : null,
    withC.length ? el("span", {}, "beat classical", el("strong", { text: `${beat.length} of ${withC.length}` })) : null].filter(Boolean));
  $("#table-empty").hidden = tableRows.length > 0;
  const head = el("tr", {}, TABLE_COLS.map((c) => {
    const th = el("th", { class: c.num ? "num" : "", "data-key": c.key, title: c.hint || "" },
      c.label, el("span", { class: "sort", text: key === c.key ? (dir > 0 ? "↑" : "↓") : "" }));
    th.addEventListener("click", () => {
      tableSort = { key: c.key, dir: tableSort.key === c.key ? -tableSort.dir : (c.num ? -1 : 1) };
      renderTable();
    });
    return th;
  }));
  const dash = "–";
  const body = rows.map((r) => {
    const isOpen = openRow === r.run_id;
    const tr = el("tr", { class: (r.stored ? "" : "gone") + (isOpen ? " open" : ""), "aria-expanded": String(isOpen),
      title: isOpen ? "Click to close" : "Click for a description" },
      el("td", { class: "name" }, r.name + (r.source === "file" ? " (file)" : ""),
        el("span", { class: "sub", text: r.title || (r.labelled ? "known site" : `from ${r.walk_from}`) })),
      el("td", { class: "num", text: String(r.residues) }),
      el("td", { class: "num" + (r.beyond != null && r.beyond >= 0.6 ? " strong-cell" : ""), text: r.beyond != null ? fmt(r.beyond) : dash }),
      el("td", { class: "num" + (r.p_value != null && r.p_value < 0.05 ? " strong-cell" : ""), text: r.p_value != null ? fmtP(r.p_value) : dash }),
      el("td", { class: "num", text: r.classical != null ? fmt(r.classical) : dash }),
      el("td", {}, r.beyond_verdict === "hump" && r.beyond_gain >= 0.01 ? el("span", { class: "chip-mini solid", text: `γ ${fmtGamma(r.beyond_gamma)}` }) : dash),
      el("td", {}, el("span", { class: "chip-mini" + (r.transport_verdict === "hump" ? "" : ""), text: verdictShort(r.transport_verdict, r.transport_gamma) })),
      el("td", { text: r.control ? "yes" : "no" }),
      el("td", { text: ago(r.finished_utc) }));
    tr.addEventListener("click", () => { openRow = isOpen ? null : r.run_id; renderTable(); });
    return isOpen ? [tr, el("tr", { class: "row-detail" }, el("td", { colspan: String(TABLE_COLS.length) }, rowDescription(r)))] : tr;
  });
  $("#run-table").replaceChildren(el("thead", {}, head), el("tbody", {}, body.flat()));
}
$("#table-filter").addEventListener("input", renderTable);
$("#table-labelled").addEventListener("change", renderTable);

// the same plain-language answers as the summary card, from the row's numbers
function rowDescription(r) {
  const lines = [];
  const say = (kind, text) => lines.push(el("li", {}, el("span", { class: "tk-mark " + kind, text: kind === "yes" ? "yes" : kind === "no" ? "no" : "partly" }), el("span", { text })));
  if (r.beyond != null) {
    const sig = r.p_value != null && r.p_value < 0.05;
    say(r.beyond >= 0.6 && (r.p_value == null || sig) ? "yes" : r.beyond > 0.55 ? "partly" : "no",
      r.beyond <= 0.55 ? `Finds the allosteric site beyond closeness: no (${fmt(r.beyond)}, where 0.5 is closeness alone).`
      : `Finds the allosteric site beyond closeness: scores ${fmt(r.beyond)}${r.p_value != null ? ` (p = ${fmtP(r.p_value)}${sig ? ", unlikely to be luck" : ", could still be luck"})` : ""}.`);
    if (r.classical != null) {
      const d = r.beyond - r.classical;
      say(d >= 0.02 && r.beyond > 0.55 ? "yes" : d > -0.02 && r.beyond > 0.55 ? "partly" : "no",
        `Needs the quantum walk: the best classical random walk scores ${fmt(r.classical)}${d >= 0.02 ? ", lower" : d > -0.02 ? ", about the same" : ", higher"}.`);
    }
    if (r.control_best != null)
      say(r.beyond > r.control_best ? "yes" : "no", `About this protein's chemistry: the best random-energy run scores ${fmt(r.control_best)}${r.beyond > r.control_best ? ", lower" : ", as high or higher"}.`);
    const helps = r.beyond_verdict === "hump" && r.beyond_gain >= 0.01 && r.beyond > 0.55;
    say(helps ? "yes" : "no", helps ? `A little noise helps find it: best at γ = ${fmtGamma(r.beyond_gamma)}, ${fmt(r.beyond_gain)} above both ends.` : "A little noise helps find it: no.");
  } else {
    lines.push(el("li", { class: "tk-hint", text: "No known allosteric site, so only transport was measured. Pick a protein with a known site to test site finding." }));
  }
  say(r.transport_verdict === "hump" ? "partly" : "no", r.transport_verdict === "hump"
    ? `Noise helps the signal travel (peak at γ = ${fmtGamma(r.transport_gamma)}), which happens in most networks.` : "Noise does not help the signal travel.");
  const open = el("button", { type: "button", class: "btn-primary small", text: r.stored ? "Open full results" : "Set up a new run" });
  open.addEventListener("click", (e) => { e.stopPropagation(); openFromTable(r); });
  const settings = `${r.residues} residues, ${r.contacts} contacts, chains ${(r.chains || []).join(", ")} · walk from ${r.walk_from} · ${r.site_energy} energies, scale ${r.scale}${r.cutoff ? `, cutoff ${r.cutoff} Å` : ""} · control ${r.control ? "on" : "off"} · took ${r.elapsed_s} s`;
  return el("div", { class: "row-desc" },
    r.title ? el("div", { class: "row-title", text: r.title }) : null,
    el("ul", { class: "row-lines" }, lines),
    el("div", { class: "row-foot" }, open, el("span", { class: "muted", text: r.stored ? settings : `Full results cleared. ${settings}` })));
}

async function openFromTable(r) {
  if (r.stored) {
    setView("walk", { scroll: false });
    await reopen(r.run_id);
    return;
  }
  // full result cleared: set the form up to run it again with the same main settings
  setView("walk", { scroll: false });
  if (r.source === "pdb") {
    pdbInput.value = r.name;
    pdbInput.dispatchEvent(new Event("input"));
  }
  const ctl = $("#control"); if (ctl) ctl.checked = !!r.control;
  showError(`${r.name}'s full results were cleared (only the last 30 are kept). ${r.source === "pdb" ? "The form is filled in; press Run the walk to redo it." : "Upload the file again to redo it."}`);
  form.scrollIntoView({ behavior: "smooth", block: "start" });
}

// typing a PDB id that is already in the table offers its stored result
function tableMatch(id) {
  return tableRows.find((r) => r.name === id && r.source === "pdb" && r.stored);
}

try { const v = localStorage.getItem("qm-view"); if (v === "table") setView("table", { scroll: false }); } catch {}
loadTable();
