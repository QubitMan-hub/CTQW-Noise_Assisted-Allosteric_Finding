"use strict";
/* The QubitMen badge, redrawn as a crisp vector from the hand-made original:
   a grey disc, a qubit wave with scattered noise dots, the name, two ring dots and
   circuit traces converging from both sides. Colours come from CSS variables, so the
   badge follows the light/dark theme. */

const LOGO = (() => {
  // qubit wave: 6.5 periods across the top, drawn as a smooth polyline
  const wave = [];
  for (let x = 50; x <= 150; x += 1) {
    const y = 45 - 5 * Math.sin(((x - 50) / 100) * Math.PI * 2 * 6.5);
    wave.push(`${x},${y.toFixed(2)}`);
  }
  // dephasing noise around the wave (positions follow the original sketch)
  const noise = [[53, 50, 2.3], [59, 38, 2.9], [66, 51, 2.0], [72, 41, 1.8], [79, 48, 2.6], [87, 36, 3.1],
    [94, 50, 2.2], [101, 40, 1.9], [108, 48, 2.9], [115, 37, 2.1], [121, 50, 1.8], [129, 40, 2.7],
    [136, 48, 2.0], [143, 37, 2.4], [148, 49, 1.8]];
  // left half of the circuit traces; the right half is its mirror image
  const traces = [
    "M60 125 H74 L80 130 H92 L99.5 135.5",
    "M40 130 H64 L69 134 H84",
    "M13 137 H30 L35 141 H68 L73 137.5 H100",
    "M44 146 H66 L71 142 H100",
    "M52 151 H78 L83 146.5 H94",
    "M63 156 H85",
  ];
  const tr = traces.map((d) => `<path d="${d}"/>`).join("");

  function svg(size = 96, title = "QubitMen") {
    return `<svg class="logo-svg" viewBox="0 0 200 200" width="${size}" height="${size}" role="img" aria-label="${title}">
  <circle class="lg-disc" cx="100" cy="100" r="97"/>
  <circle class="lg-rim" cx="100" cy="100" r="92.5" fill="none" stroke-width="1"/>
  <polyline class="lg-line" points="${wave.join(" ")}" fill="none" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/>
  <g class="lg-noise">${noise.map(([x, y, r]) => `<circle cx="${x}" cy="${y}" r="${r}"/>`).join("")}</g>
  <text class="lg-text" x="100" y="83" text-anchor="middle" font-size="29">Qubit</text>
  <text class="lg-text" x="100" y="112" text-anchor="middle" font-size="29">Men</text>
  <circle class="lg-ring" cx="47" cy="92" r="5.5" stroke-width="2.2"/>
  <circle class="lg-ring" cx="153" cy="92" r="5.5" stroke-width="2.2"/>
  <g class="lg-line" fill="none" stroke-width="2.3" stroke-linecap="round" stroke-linejoin="round">
    ${tr}<g transform="translate(200 0) scale(-1 1)">${tr}</g>
  </g>
</svg>`;
  }
  return { svg };
})();
