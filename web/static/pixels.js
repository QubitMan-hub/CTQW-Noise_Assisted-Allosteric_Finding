"use strict";
/* Pixel art for QubitMan: the wordmark, the icons, the background network and the
   loading walker. Everything is greyscale and drawn from small bitmaps below
   ('#' ink pixel, '+' mid grey, '.' empty); colours come from the CSS theme. */

const PX = (() => {
  const NS = "http://www.w3.org/2000/svg";
  const reduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* ---------- bitmap font (7 rows) ---------- */
  const GLYPHS = {
    Q: [".###.", "#...#", "#...#", "#...#", "#.#.#", "#..#.", ".##.#"],
    u: [".....", ".....", "#...#", "#...#", "#...#", "#..##", ".##.#"],
    b: ["#....", "#....", "####.", "#...#", "#...#", "#...#", "####."],
    i: ["#", ".", "#", "#", "#", "#", "#"],
    t: [".#..", ".#..", "####", ".#..", ".#..", ".#.#", "..#."],
    M: ["#...#", "##.##", "#.#.#", "#.#.#", "#...#", "#...#", "#...#"],
    a: [".....", ".....", ".###.", "....#", ".####", "#...#", ".####"],
    n: [".....", ".....", "####.", "#...#", "#...#", "#...#", "#...#"],
  };

  /* ---------- 8x8 icons ---------- */
  const ICONS = {
    dice:  ["########", "#......#", "#.#..#.#", "#......#", "#..##..#", "#......#", "#.#..#.#", "########"],
    key:   ["..###...", ".#...#..", ".#...#..", "..###...", "...#....", "...##...", "...#....", "...##..."],
    file:  ["#####...", "#...##..", "#....##.", "#.####.#", "#......#", "#.####.#", "#......#", "########"],
    helix: [".######.", "#......#", ".######.", "#......#", ".######.", "#......#", ".######.", "........"],
    net:   ["##....##", "##+..+##", "..+..+..", "...##...", "...##...", "..+..+..", "##+..+##", "##....##"],
    wave:  ["........", ".##.....", "#..#....", "....#..#", ".....##.", "........", "+.+.+.+.", "........"],
    qubit: ["..####..", ".#....#.", "#...#..#", "#..###.#", "#.#.#.##", "#...#..#", ".#....#.", "..####.."],
    moon:  ["..####..", ".###....", "###.....", "###.....", "###.....", "###.....", ".###....", "..####.."],
    sun:   ["...#....", ".#.#.#..", "..###...", "#######.", "..###...", ".#.#.#..", "...#....", "........"],
  };

  function svgEl(tag, attrs, parent) {
    const el = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, v);
    if (parent) parent.appendChild(el);
    return el;
  }

  function drawBitmap(svg, rows, px, { gap = 0, animate = false, delayStep = 6 } = {}) {
    const w = Math.max(...rows.map((r) => r.length)), h = rows.length;
    svg.setAttribute("viewBox", `0 0 ${w * px} ${h * px}`);
    svg.setAttribute("shape-rendering", "crispEdges");
    svg.replaceChildren();
    const cells = [];
    rows.forEach((row, y) => [...row].forEach((ch, x) => {
      if (ch === "." || ch === " ") return;
      const r = svgEl("rect", { x: x * px + gap / 2, y: y * px + gap / 2, width: px - gap, height: px - gap,
        class: ch === "+" ? "px mid" : "px" }, svg);
      if (animate && !reduced) {
        r.style.animationDelay = `${(x + y) * delayStep}ms`;
        r.classList.add("px-in");
      }
      cells.push({ el: r, x, y });
    }));
    return { cells, w, h };
  }

  /* ---------- wordmark ---------- */
  function wordmark(svg, text = "QubitMan", px = 7) {
    const rows = Array.from({ length: 7 }, () => "");
    for (const ch of text) {
      const g = GLYPHS[ch];
      g.forEach((line, i) => { rows[i] += line + "."; });
    }
    const { cells } = drawBitmap(svg, rows.map((r) => r.slice(0, -1)), px, { gap: 1, animate: true, delayStep: 9 });
    if (reduced) return;
    // pixels near the pointer lift a little: a ripple you can play with
    svg.addEventListener("pointermove", (e) => {
      const box = svg.getBoundingClientRect();
      const scale = box.width / svg.viewBox.baseVal.width;
      const mx = (e.clientX - box.left) / scale / px, my = (e.clientY - box.top) / scale / px;
      for (const c of cells) {
        const d = Math.hypot(c.x + 0.5 - mx, c.y + 0.5 - my);
        c.el.style.transform = d < 3.2 ? `translateY(${-(3.2 - d) * 1.2}px)` : "";
        c.el.style.fill = d < 1.6 ? "var(--muted)" : "";
      }
    });
    svg.addEventListener("pointerleave", () => cells.forEach((c) => { c.el.style.transform = ""; c.el.style.fill = ""; }));
  }

  function icons(root = document) {
    root.querySelectorAll("svg.px-icon[data-icon]").forEach((s) => drawBitmap(s, ICONS[s.dataset.icon] || ICONS.qubit, 3));
  }

  /* ---------- background: a faint contact network with a walker hopping on it ---------- */
  function background(canvas) {
    const ctx = canvas.getContext("2d");
    let W, H, nodes = [], walker = 0, trail = [], last = 0, hopAt = 0, running = true;
    function build() {
      const dpr = Math.min(2, window.devicePixelRatio || 1);
      W = canvas.width = Math.floor(innerWidth * dpr);
      H = canvas.height = Math.floor(innerHeight * dpr);
      canvas.style.width = innerWidth + "px"; canvas.style.height = innerHeight + "px";
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const n = Math.round(Math.min(90, (innerWidth * innerHeight) / 16000));
      nodes = Array.from({ length: n }, () => ({ x: Math.random() * innerWidth, y: Math.random() * innerHeight,
        vx: (Math.random() - 0.5) * 0.08, vy: (Math.random() - 0.5) * 0.08 }));
      walker = 0; trail = [];
    }
    function neighbours(i) {
      const a = nodes[i], out = [];
      nodes.forEach((b, j) => { if (j !== i && Math.hypot(a.x - b.x, a.y - b.y) < 150) out.push(j); });
      return out;
    }
    function frame(t) {
      if (!running) return;
      const ink = getComputedStyle(document.documentElement).getPropertyValue("--ink-rgb").trim() || "22, 22, 22";
      const dt = Math.min(50, t - (last || t)); last = t;
      ctx.clearRect(0, 0, innerWidth, innerHeight);
      for (const p of nodes) {
        if (!reduced) { p.x += p.vx * dt; p.y += p.vy * dt; }
        if (p.x < -20) p.x = innerWidth + 20; if (p.x > innerWidth + 20) p.x = -20;
        if (p.y < -20) p.y = innerHeight + 20; if (p.y > innerHeight + 20) p.y = -20;
      }
      ctx.lineWidth = 1;
      for (let i = 0; i < nodes.length; i++) for (let j = i + 1; j < nodes.length; j++) {
        const a = nodes[i], b = nodes[j], d = Math.hypot(a.x - b.x, a.y - b.y);
        if (d < 150) {
          ctx.strokeStyle = `rgba(${ink},${0.07 * (1 - d / 150)})`;
          ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
        }
      }
      if (!reduced && t > hopAt) {
        const nb = neighbours(walker);
        trail.push({ i: walker, t });
        walker = nb.length ? nb[Math.floor(Math.random() * nb.length)] : Math.floor(Math.random() * nodes.length);
        hopAt = t + 700;
      }
      trail = trail.filter((s) => t - s.t < 5000);
      for (const p of nodes) { ctx.fillStyle = `rgba(${ink},0.10)`; ctx.fillRect(p.x - 1.5, p.y - 1.5, 3, 3); }
      for (const s of trail) {
        const p = nodes[s.i], a = 0.22 * (1 - (t - s.t) / 5000);
        ctx.fillStyle = `rgba(${ink},${a})`; ctx.fillRect(p.x - 2.5, p.y - 2.5, 5, 5);
      }
      const w = nodes[walker];
      if (w) { ctx.fillStyle = `rgba(${ink},0.35)`; ctx.fillRect(w.x - 3.5, w.y - 3.5, 7, 7); }
      if (!reduced) requestAnimationFrame(frame);
    }
    build();
    addEventListener("resize", () => { build(); if (reduced) frame(performance.now()); });
    document.addEventListener("visibilitychange", () => {
      running = !document.hidden;
      if (running && !reduced) { last = 0; requestAnimationFrame(frame); }
    });
    requestAnimationFrame(frame);
  }

  /* ---------- loading: a pixel walker hopping along a chain of residues ---------- */
  function walkerBar(host, n = 28) {
    host.replaceChildren(...Array.from({ length: n }, () => document.createElement("span")));
    const cells = [...host.children];
    let pos = Math.floor(n / 2), timer = null;
    const heat = new Array(n).fill(0);
    function step() {
      pos = Math.max(0, Math.min(n - 1, pos + (Math.random() < 0.5 ? -1 : 1) * (Math.random() < 0.25 ? 2 : 1)));
      for (let i = 0; i < n; i++) heat[i] *= 0.82;
      heat[pos] = 1;
      cells.forEach((c, i) => { c.style.opacity = String(0.12 + 0.88 * heat[i]); });
    }
    return {
      start() { clearInterval(timer); if (!reduced) timer = setInterval(step, 110); step(); },
      stop() { clearInterval(timer); },
    };
  }

  return { wordmark, icons, background, walkerBar, drawBitmap, reduced, ICONS };
})();
