/**
 * InSAR Multilooking Visualizer — main.js
 *
 * Interactive demonstration of InSAR multilooking concepts.
 * Based on: S. Samiei-Esfahany, "Exploitation of Distributed Scatterers
 * in Synthetic Aperture Radar Interferometry", PhD Thesis, TU Delft, 2017.
 *
 * Author: Generated with Claude (Anthropic)
 * License: MIT
 *
 * Physics implemented:
 *  - SAR SLC phase model: psi = psi_range + psi_atmo + psi_scat + psi_noise
 *  - Complex circular Gaussian noise model for distributed scatterers
 *  - Phase variance: sigma^2 = (1 - |gamma|^2) / (2*L*|gamma|^2)  [Cramér-Rao bound]
 *  - Boxcar multilooking: phi_hat = arg[ (1/L) sum(P_M * conj(P_S)) ]
 *  - Gaussian-weighted multilooking
 *  - Adaptive multilooking with brotherhood selection (parametric mean test)
 */

'use strict';

// ─────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────

/** Pixel grid dimension (N×N pixels) */
const N = 64;

// ─────────────────────────────────────────────────────────────
// State
// ─────────────────────────────────────────────────────────────

let rawPhase = null;
let mlPhase = null;
let signalPhase = null;
let noisePhase = null;
let brotherhoodMask = null; // Uint8Array for adaptive mode

/** Currently inspected pixel coordinates */
let selectedPx = { x: 24, y: 24 };

// ─────────────────────────────────────────────────────────────
// Phase → Color mapping (HSL wheel: -pi → purple, 0 → green, +pi → orange)
// ─────────────────────────────────────────────────────────────

/**
 * Convert a phase value in [-pi, pi] to an RGB color triple.
 * Uses a perceptually meaningful colormap (cyclic HSL-based).
 * @param {number} phi - Phase value in radians
 * @returns {number[]} [r, g, b] each in [0, 255]
 */
function phaseToColor(phi) {
  // Normalize to [0, 1]
  const t = (phi + Math.PI) / (2 * Math.PI);
  // Cyclic hue mapping
  const hue = t * 360;
  const s = 0.75;
  const l = 0.48;

  const h = hue / 60;
  const i = Math.floor(h);
  const f = h - i;
  const p = l * (1 - s);
  const q = l * (1 - s * f);
  const u = l * (1 - s * (1 - f));

  let r, g, b;
  switch (i % 6) {
    case 0: r = l; g = u; b = p; break;
    case 1: r = q; g = l; b = p; break;
    case 2: r = p; g = l; b = u; break;
    case 3: r = p; g = q; b = l; break;
    case 4: r = u; g = p; b = l; break;
    default: r = l; g = p; b = q; break;
  }

  return [
    Math.round(Math.max(0, Math.min(255, r * 255))),
    Math.round(Math.max(0, Math.min(255, g * 255))),
    Math.round(Math.max(0, Math.min(255, b * 255))),
  ];
}

// ─────────────────────────────────────────────────────────────
// Math utilities
// ─────────────────────────────────────────────────────────────

/**
 * Wrap a phase value to the interval (-pi, pi].
 * @param {number} phi
 * @returns {number}
 */
function wrap(phi) {
  while (phi > Math.PI) phi -= 2 * Math.PI;
  while (phi < -Math.PI) phi += 2 * Math.PI;
  return phi;
}

/**
 * Compute the circular standard deviation of a phase array.
 * Uses the circular mean as reference to avoid wrap-around artifacts.
 * @param {number[]} arr - Array of phase values in radians
 * @returns {number} Standard deviation in radians
 */
function circularStdDev(arr) {
  const n = arr.length;
  if (n === 0) return 0;
  const meanSin = arr.reduce((s, v) => s + Math.sin(v), 0) / n;
  const meanCos = arr.reduce((s, v) => s + Math.cos(v), 0) / n;
  const circMean = Math.atan2(meanSin, meanCos);
  const diffs = arr.map(v => wrap(v - circMean));
  return Math.sqrt(diffs.reduce((s, v) => s + v * v, 0) / n);
}

/**
 * Generate Box-Muller normally distributed random number.
 * @returns {number} Standard normal sample
 */
function randn() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// ─────────────────────────────────────────────────────────────
// Signal generation
// ─────────────────────────────────────────────────────────────

/**
 * Generate the true deformation signal phase field.
 * Models a smooth sinusoidal deformation pattern (e.g., subsidence bowl).
 *
 * @param {number} freq - Spatial frequency of the signal
 * @returns {Float32Array} N×N signal phase in radians
 */
function generateSignalPhase(freq) {
  const signal = new Float32Array(N * N);
  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const fx = x / N;
      const fy = y / N;
      // Superposition of two sinusoids in range and azimuth directions
      signal[y * N + x] = Math.PI * freq * (
        0.7 * Math.sin(2 * Math.PI * fx) +
        0.5 * Math.cos(2 * Math.PI * fy * 0.8)
      );
    }
  }
  return signal;
}

/**
 * Generate complex circular Gaussian noise based on coherence model.
 *
 * For a distributed scatterer with coherence |gamma|, the noise
 * standard deviation follows from the Cramér-Rao bound:
 *   sigma ≈ sqrt((1 - |gamma|^2) / (2 * |gamma|^2))
 *
 * @param {number} coherence - Absolute coherence |gamma| in [0, 1]
 * @returns {Float32Array} N×N phase noise in radians
 */
function generateNoise(coherence) {
  const noise = new Float32Array(N * N);
  // Noise sigma from coherence (single-look, L=1)
  const sigma = Math.sqrt((1 - coherence * coherence) / (2 * coherence * coherence + 1e-9));
  for (let i = 0; i < N * N; i++) {
    // Complex circular Gaussian: noise in real + imaginary parts
    const noiseRe = randn() * sigma;
    const noiseIm = randn() * sigma;
    // Phase noise = angle of (1 + noise_complex) ≈ noise_Im for small noise
    noise[i] = Math.atan2(noiseIm, 1 + noiseRe) * 2.5;
  }
  return noise;
}

// ─────────────────────────────────────────────────────────────
// Multilooking algorithms
// ─────────────────────────────────────────────────────────────

/**
 * Standard boxcar (uniform) complex multilooking.
 *
 * Computes the complex average of all pixels within a rectangular window:
 *   <I_k> = (1/L) * sum_{i in window} P_M,i * conj(P_S,i)
 *
 * The multilooked phase is then:
 *   phi_hat = arg(<I_k>)
 *
 * This is the Maximum Likelihood estimator of the interferometric phase
 * for distributed scatterers (Rodriguez & Martin, 1992).
 *
 * @param {Float32Array} rawRe - Real parts of complex interferogram (N×N)
 * @param {Float32Array} rawIm - Imaginary parts of complex interferogram (N×N)
 * @param {number} winSize - Window side length (odd integer)
 * @param {Float32Array|null} kernel - Optional weight kernel (winSize×winSize)
 * @returns {{ re: Float32Array, im: Float32Array }}
 */
function boxcarMultilook(rawRe, rawIm, winSize, kernel = null) {
  const out = { re: new Float32Array(N * N), im: new Float32Array(N * N) };
  const half = Math.floor(winSize / 2);

  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      let sumRe = 0, sumIm = 0, sumW = 0;

      for (let dy = -half; dy <= half; dy++) {
        for (let dx = -half; dx <= half; dx++) {
          const nx = x + dx, ny = y + dy;
          if (nx < 0 || nx >= N || ny < 0 || ny >= N) continue;

          const w = kernel
            ? kernel[(dy + half) * winSize + (dx + half)]
            : 1.0;
          const idx = ny * N + nx;
          sumRe += rawRe[idx] * w;
          sumIm += rawIm[idx] * w;
          sumW += w;
        }
      }

      const idx = y * N + x;
      out.re[idx] = sumW > 0 ? sumRe / sumW : 0;
      out.im[idx] = sumW > 0 ? sumIm / sumW : 0;
    }
  }

  return out;
}

/**
 * Build a 2D Gaussian kernel for weighted multilooking.
 * Center pixel receives maximum weight; weights decay as exp(-r²/2σ²).
 *
 * @param {number} size - Kernel side length (odd integer)
 * @returns {Float32Array} Normalized Gaussian kernel (size×size)
 */
function makeGaussianKernel(size) {
  const kernel = new Float32Array(size * size);
  const half = Math.floor(size / 2);
  const sigma = size / 3.0;
  let sum = 0;

  for (let dy = -half; dy <= half; dy++) {
    for (let dx = -half; dx <= half; dx++) {
      const w = Math.exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
      kernel[(dy + half) * size + (dx + half)] = w;
      sum += w;
    }
  }

  // Normalize so weights sum to 1
  for (let i = 0; i < kernel.length; i++) kernel[i] /= sum;
  return kernel;
}

/**
 * Adaptive multilooking with brotherhood selection.
 *
 * For each pixel p, tests neighboring pixels q to determine if they
 * share the same amplitude statistics (null hypothesis: same distribution).
 * Uses the parametric mean test (Jiang et al., 2014a):
 *
 *   H0: E{A_bar_q} = mu_Ap  (q is a "brother" of p)
 *   T_mean = sqrt(N) * (A_bar_q - mu_hat_Ap) / sigma_hat_Ap
 *   Reject if |T_mean| > k_alpha (critical value, here alpha=0.05 → k=1.96)
 *
 * Only brotherhood pixels are included in the complex average.
 * This prevents mixing heterogeneous terrain types.
 *
 * Reference: Ferretti et al. (2011b) SqueeSAR; Jiang et al. (2014a)
 *
 * @param {Float32Array} rawRe - Real parts of complex interferogram
 * @param {Float32Array} rawIm - Imaginary parts
 * @param {number} winSize - Window side length
 * @param {Float32Array} amplitudes - Per-pixel amplitude estimates
 * @param {number} selX - Selected pixel x (for brotherhood mask output)
 * @param {number} selY - Selected pixel y
 * @returns {{ re: Float32Array, im: Float32Array, mask: Uint8Array }}
 */
function adaptiveMultilook(rawRe, rawIm, winSize, amplitudes, selX, selY) {
  const out = { re: new Float32Array(N * N), im: new Float32Array(N * N) };
  const half = Math.floor(winSize / 2);
  const maskSize = winSize * winSize;
  const mask = new Uint8Array(maskSize); // brotherhood mask for selected pixel
  const kAlpha = 1.96; // 5% significance level (two-tailed z-test)

  for (let y = 0; y < N; y++) {
    for (let x = 0; x < N; x++) {
      const pidx = y * N + x;
      const muP = amplitudes[pidx];
      // Coefficient of variation for Rayleigh distribution ≈ 0.52
      const sigmaP = muP * 0.52 + 1e-9;

      let sumRe = 0, sumIm = 0, count = 0;
      const isSelected = (x === selX && y === selY);

      for (let dy = -half; dy <= half; dy++) {
        for (let dx = -half; dx <= half; dx++) {
          const nx = x + dx, ny = y + dy;
          if (nx < 0 || nx >= N || ny < 0 || ny >= N) continue;

          const qidx = ny * N + nx;
          const muQ = amplitudes[qidx];

          // Parametric mean test statistic
          const T = Math.abs((muQ - muP) * Math.sqrt(N / 4)) / sigmaP;
          const isBrother = (T <= kAlpha);

          const kidx = (dy + half) * winSize + (dx + half);
          if (isSelected) mask[kidx] = isBrother ? 1 : 0;

          if (isBrother) {
            sumRe += rawRe[qidx];
            sumIm += rawIm[qidx];
            count++;
          }
        }
      }

      if (count === 0) {
        // Fallback: use just the pixel itself
        out.re[pidx] = rawRe[pidx];
        out.im[pidx] = rawIm[pidx];
      } else {
        out.re[pidx] = sumRe / count;
        out.im[pidx] = sumIm / count;
      }
    }
  }

  return { re: out.re, im: out.im, mask };
}

// ─────────────────────────────────────────────────────────────
// Canvas rendering
// ─────────────────────────────────────────────────────────────

/**
 * Draw a phase image onto a canvas.
 * @param {HTMLCanvasElement} canvas
 * @param {Float32Array} phases - N×N phase array in radians
 * @param {{ x: number, y: number }|null} highlight - Pixel to highlight
 */
function drawPhaseImage(canvas, phases, highlight = null) {
  const ctx = canvas.getContext('2d');
  const imgData = ctx.createImageData(N, N);

  for (let i = 0; i < N * N; i++) {
    const [r, g, b] = phaseToColor(phases[i]);
    imgData.data[i * 4 + 0] = r;
    imgData.data[i * 4 + 1] = g;
    imgData.data[i * 4 + 2] = b;
    imgData.data[i * 4 + 3] = 255;
  }

  ctx.putImageData(imgData, 0, 0);

  if (highlight) {
    const scale = canvas.width / N;
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.9)';
    ctx.lineWidth = 1.5 / scale;
    ctx.strokeRect(highlight.x, highlight.y, 1, 1);
    // Draw crosshair
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.4)';
    ctx.lineWidth = 0.5 / scale;
    ctx.beginPath();
    ctx.moveTo(highlight.x + 0.5, 0);
    ctx.lineTo(highlight.x + 0.5, N);
    ctx.moveTo(0, highlight.y + 0.5);
    ctx.lineTo(N, highlight.y + 0.5);
    ctx.stroke();
  }
}

/**
 * Draw the zoomed window showing the pixel neighborhood.
 * @param {HTMLCanvasElement} canvas
 * @param {number[][]} window2d - 2D array of phase values (winSize × winSize)
 * @param {Uint8Array|null} brotherhoodMask - 1 = brother, 0 = excluded
 */
function drawZoomWindow(canvas, window2d, brotherhoodMask) {
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const rows = window2d.length;
  const cols = window2d[0].length;
  const cellW = Math.floor(W / cols);
  const cellH = Math.floor(H / rows);

  ctx.clearRect(0, 0, W, H);

  for (let dy = 0; dy < rows; dy++) {
    for (let dx = 0; dx < cols; dx++) {
      const phi = window2d[dy][dx];
      const [r, g, b] = phaseToColor(phi);
      const kidx = dy * cols + dx;
      const isBrother = brotherhoodMask ? brotherhoodMask[kidx] === 1 : true;

      if (isBrother) {
        ctx.fillStyle = `rgb(${r},${g},${b})`;
      } else {
        // Excluded pixels: dimmed with reddish overlay
        ctx.fillStyle = `rgba(${r},${g},${b},0.2)`;
      }

      ctx.fillRect(dx * cellW, dy * cellH, cellW, cellH);

      if (!isBrother) {
        ctx.strokeStyle = 'rgba(224,82,82,0.6)';
        ctx.lineWidth = 1;
        ctx.strokeRect(dx * cellW + 0.5, dy * cellH + 0.5, cellW - 1, cellH - 1);
      }
    }
  }

  // Grid lines
  ctx.strokeStyle = 'rgba(0, 0, 0, 0.25)';
  ctx.lineWidth = 0.5;
  for (let i = 0; i <= rows; i++) {
    ctx.beginPath(); ctx.moveTo(0, i * cellH); ctx.lineTo(W, i * cellH); ctx.stroke();
  }
  for (let j = 0; j <= cols; j++) {
    ctx.beginPath(); ctx.moveTo(j * cellW, 0); ctx.lineTo(j * cellW, H); ctx.stroke();
  }

  // Highlight center pixel (the selected pixel)
  const cx = Math.floor(cols / 2);
  const cy = Math.floor(rows / 2);
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.95)';
  ctx.lineWidth = 2;
  ctx.strokeRect(cx * cellW + 1, cy * cellH + 1, cellW - 2, cellH - 2);
}

/**
 * Draw the complex phasor diagram.
 * Shows each pixel in the window as a unit complex vector,
 * plus the resultant (multilooked estimate) as the orange vector.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {number[][]} window2d - Phase values in averaging window
 * @param {Uint8Array|null} brotherhoodMask
 */
function drawPhasorDiagram(canvas, window2d, brotherhoodMask) {
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const cx = W / 2, cy = H / 2;
  const R = Math.min(W, H) * 0.38;

  // Unit circle
  ctx.strokeStyle = 'rgba(255,255,255,0.08)';
  ctx.lineWidth = 0.5;
  ctx.beginPath(); ctx.arc(cx, cy, R, 0, Math.PI * 2); ctx.stroke();

  // Axes
  ctx.strokeStyle = 'rgba(255,255,255,0.1)';
  ctx.lineWidth = 0.5;
  ctx.beginPath();
  ctx.moveTo(cx - R * 1.15, cy); ctx.lineTo(cx + R * 1.15, cy);
  ctx.moveTo(cx, cy - R * 1.15); ctx.lineTo(cx, cy + R * 1.15);
  ctx.stroke();

  // Axis labels
  ctx.fillStyle = 'rgba(255,255,255,0.3)';
  ctx.font = '10px JetBrains Mono, monospace';
  ctx.textAlign = 'left';
  ctx.fillText('Re', cx + R * 1.05, cy - 4);
  ctx.textAlign = 'center';
  ctx.fillText('Im', cx, cy - R * 1.12);

  // Individual phasors
  let sumRe = 0, sumIm = 0, brotherCount = 0;
  const rows = window2d.length;
  const cols = window2d[0].length;

  window2d.forEach((row, dy) => {
    row.forEach((phi, dx) => {
      const kidx = dy * cols + dx;
      const isBrother = brotherhoodMask ? brotherhoodMask[kidx] === 1 : true;
      const alpha = isBrother ? 0.5 : 0.1;

      ctx.strokeStyle = `rgba(79, 142, 247, ${alpha})`;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(cx, cy);
      ctx.lineTo(cx + Math.cos(phi) * R * 0.78, cy - Math.sin(phi) * R * 0.78);
      ctx.stroke();

      if (isBrother) {
        sumRe += Math.cos(phi);
        sumIm += Math.sin(phi);
        brotherCount++;
      }
    });
  });

  // Resultant (multilooked) vector
  if (brotherCount > 0) {
    const magnitude = Math.sqrt(sumRe * sumRe + sumIm * sumIm) / brotherCount;
    const resultAngle = Math.atan2(sumIm, sumRe);
    const endX = cx + Math.cos(resultAngle) * R * 0.85;
    const endY = cy - Math.sin(resultAngle) * R * 0.85;

    // Arrow shaft
    ctx.strokeStyle = 'rgba(240, 164, 41, 0.95)';
    ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(endX, endY); ctx.stroke();

    // Arrowhead
    const angle = Math.atan2(endY - cy, endX - cx);
    const arrowLen = 8;
    ctx.fillStyle = 'rgba(240, 164, 41, 0.95)';
    ctx.beginPath();
    ctx.moveTo(endX, endY);
    ctx.lineTo(endX - arrowLen * Math.cos(angle - 0.4), endY - arrowLen * Math.sin(angle - 0.4));
    ctx.lineTo(endX - arrowLen * Math.cos(angle + 0.4), endY - arrowLen * Math.sin(angle + 0.4));
    ctx.closePath(); ctx.fill();

    // Info text
    ctx.fillStyle = 'rgba(255,255,255,0.35)';
    ctx.font = '10px JetBrains Mono, monospace';
    ctx.textAlign = 'left';
    ctx.fillText(`n = ${brotherCount} brothers`, 6, H - 18);
    ctx.fillText(`|result| = ${magnitude.toFixed(3)}`, 6, H - 6);

    // Update DOM
    document.getElementById('phasorN').textContent = `n = ${brotherCount} brothers`;
    document.getElementById('phasorMag').textContent = `|result| = ${magnitude.toFixed(3)}`;
  }
}

/**
 * Draw a phase histogram for the pixels in the averaging window.
 * Shows the distribution of raw phases (blue bars) and marks
 * the multilooked estimate (orange vertical line).
 *
 * @param {HTMLCanvasElement} canvas
 * @param {number[][]} window2d - Phase values
 * @param {number} mlPhaseValue - Final multilooked phase estimate
 */
function drawPhaseHistogram(canvas, window2d, mlPhaseValue) {
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const BINS = 24;
  const binWidth = (2 * Math.PI) / BINS;
  const counts = new Array(BINS).fill(0);

  window2d.flat().forEach(phi => {
    const bin = Math.min(BINS - 1, Math.floor((phi + Math.PI) / binWidth));
    counts[bin]++;
  });

  const maxCount = Math.max(...counts, 1);
  const padL = 10, padR = 10, padT = 10, padB = 28;
  const plotW = W - padL - padR;
  const plotH = H - padT - padB;
  const bW = plotW / BINS;

  // Bars
  counts.forEach((count, i) => {
    const barH = (count / maxCount) * plotH;
    ctx.fillStyle = 'rgba(79,142,247,0.4)';
    ctx.fillRect(padL + i * bW, padT + plotH - barH, bW - 1, barH);
    ctx.strokeStyle = 'rgba(79,142,247,0.7)';
    ctx.lineWidth = 0.5;
    ctx.strokeRect(padL + i * bW, padT + plotH - barH, bW - 1, barH);
  });

  // X-axis baseline
  ctx.strokeStyle = 'rgba(255,255,255,0.1)';
  ctx.lineWidth = 0.5;
  ctx.beginPath();
  ctx.moveTo(padL, padT + plotH); ctx.lineTo(W - padR, padT + plotH);
  ctx.stroke();

  // ML estimate line
  const mlX = padL + ((mlPhaseValue + Math.PI) / (2 * Math.PI)) * plotW;
  ctx.strokeStyle = 'rgba(240,164,41,0.95)';
  ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(mlX, padT); ctx.lineTo(mlX, padT + plotH); ctx.stroke();

  // Labels
  ctx.fillStyle = 'rgba(255,255,255,0.3)';
  ctx.font = '10px JetBrains Mono, monospace';
  ctx.textAlign = 'left';  ctx.fillText('−π', padL, H - 4);
  ctx.textAlign = 'right'; ctx.fillText('+π', W - padR, H - 4);
  ctx.textAlign = 'center'; ctx.fillText('0', W / 2, H - 4);

  ctx.fillStyle = 'rgba(240,164,41,0.85)';
  ctx.fillText('ML', mlX, padT + plotH + 13);
}

/**
 * Draw the phase colorbar legend.
 * @param {HTMLCanvasElement} canvas
 */
function drawColorbar(canvas) {
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const imgData = ctx.createImageData(W, H);

  for (let x = 0; x < W; x++) {
    const phi = -Math.PI + (x / W) * 2 * Math.PI;
    const [r, g, b] = phaseToColor(phi);
    for (let y = 0; y < H; y++) {
      const i = (y * W + x) * 4;
      imgData.data[i] = r; imgData.data[i+1] = g; imgData.data[i+2] = b; imgData.data[i+3] = 255;
    }
  }

  ctx.putImageData(imgData, 0, 0);
}

// ─────────────────────────────────────────────────────────────
// Zoom / detail panel update
// ─────────────────────────────────────────────────────────────

/**
 * Update all detail panels (zoom, phasor, histogram) for the selected pixel.
 * @param {number} winSize - Window side length
 * @param {string} method - 'boxcar' | 'gaussian' | 'adaptive'
 */
function updateDetailPanels(winSize, method) {
  const half = Math.floor(winSize / 2);
  const px = selectedPx.x, py = selectedPx.y;

  // Build 2D window array of raw phases
  const window2d = [];
  for (let dy = -half; dy <= half; dy++) {
    const row = [];
    for (let dx = -half; dx <= half; dx++) {
      const nx = Math.max(0, Math.min(N - 1, px + dx));
      const ny = Math.max(0, Math.min(N - 1, py + dy));
      row.push(rawPhase[ny * N + nx]);
    }
    window2d.push(row);
  }

  // Brotherhood mask (only meaningful for adaptive method)
  const maskForDisplay = method === 'adaptive' ? brotherhoodMask : null;

  drawZoomWindow(document.getElementById('canvasZoom'), window2d, maskForDisplay);
  drawPhasorDiagram(document.getElementById('canvasPhasor'), window2d, maskForDisplay);
  drawPhaseHistogram(document.getElementById('canvasHist'), window2d, mlPhase[py * N + px]);

  // Update pixel info text
  document.getElementById('pixX').textContent = px;
  document.getElementById('pixY').textContent = py;
  document.getElementById('rawPhiVal').textContent = rawPhase[py * N + px].toFixed(3);
  document.getElementById('mlPhiVal').textContent = mlPhase[py * N + px].toFixed(3);

  // Show/hide brotherhood note
  document.getElementById('brotherhoodNote').style.display =
    method === 'adaptive' ? 'inline' : 'none';
}

// ─────────────────────────────────────────────────────────────
// Main render function
// ─────────────────────────────────────────────────────────────

/**
 * Full re-render: regenerate signal+noise, run selected multilooking
 * algorithm, update all canvases and statistics.
 */
function render() {
  const winSize = parseInt(document.getElementById('winSize').value);
  const coherence = parseFloat(document.getElementById('coherence').value);
  const sigFreq = parseFloat(document.getElementById('sigFreq').value);
  const method = document.querySelector('input[name="mlMethod"]:checked').value;

  // Update control labels
  document.getElementById('winSizeVal').textContent = `${winSize}×${winSize}`;
  document.getElementById('coherenceVal').textContent = coherence.toFixed(2);
  document.getElementById('sigFreqVal').textContent = sigFreq.toFixed(1);
  document.getElementById('mlPanelTitle').textContent =
    `Multilooked phase (${winSize}×${winSize} ${method})`;
  document.getElementById('arrowLabel').textContent =
    `${winSize}×${winSize}\n${method}`;

  // ── 1. Generate signal + noise ──────────────────────────────
  signalPhase = generateSignalPhase(sigFreq);
  noisePhase = generateNoise(coherence);
  rawPhase = new Float32Array(N * N);

  const rawRe = new Float32Array(N * N);
  const rawIm = new Float32Array(N * N);
  const amplitudes = new Float32Array(N * N); // for adaptive brotherhood test

  for (let i = 0; i < N * N; i++) {
    rawPhase[i] = wrap(signalPhase[i] + noisePhase[i]);
    rawRe[i] = Math.cos(rawPhase[i]);
    rawIm[i] = Math.sin(rawPhase[i]);
    // Simulate amplitude: mean ~ 1, Rayleigh distributed
    amplitudes[i] = Math.abs(rawIm[i]) * 0.6 + 0.4 + Math.abs(noisePhase[i]) * 0.15;
  }

  // ── 2. Run multilooking ────────────────────────────────────
  brotherhoodMask = null;
  let mlResult;

  if (method === 'boxcar') {
    mlResult = boxcarMultilook(rawRe, rawIm, winSize);
  } else if (method === 'gaussian') {
    const kernel = makeGaussianKernel(winSize);
    mlResult = boxcarMultilook(rawRe, rawIm, winSize, kernel);
  } else {
    // Adaptive
    const res = adaptiveMultilook(rawRe, rawIm, winSize, amplitudes, selectedPx.x, selectedPx.y);
    mlResult = res;
    brotherhoodMask = res.mask;
  }

  // Extract multilooked phases
  mlPhase = new Float32Array(N * N);
  for (let i = 0; i < N * N; i++) {
    mlPhase[i] = Math.atan2(mlResult.im[i], mlResult.re[i]);
  }

  // ── 3. Compute residual ────────────────────────────────────
  const residPhase = new Float32Array(N * N);
  for (let i = 0; i < N * N; i++) {
    residPhase[i] = wrap(rawPhase[i] - mlPhase[i]);
  }

  // ── 4. Draw phase images ───────────────────────────────────
  drawPhaseImage(document.getElementById('canvasRaw'), rawPhase, selectedPx);
  drawPhaseImage(document.getElementById('canvasML'), mlPhase, selectedPx);
  drawPhaseImage(document.getElementById('canvasResid'), residPhase);

  // ── 5. Update statistics ───────────────────────────────────
  const rawArr = Array.from(rawPhase);
  const mlArr = Array.from(mlPhase);
  const rawStd = circularStdDev(rawArr);
  const mlStd = circularStdDev(mlArr);
  const noiseReduction = rawStd > 0 ? ((rawStd - mlStd) / rawStd * 100) : 0;
  const effLooks = rawStd > 0 && mlStd > 0 ? (rawStd * rawStd) / (mlStd * mlStd) : 1;

  document.getElementById('stat-raw').textContent = rawStd.toFixed(3);
  document.getElementById('stat-ml').textContent = mlStd.toFixed(3);
  document.getElementById('stat-nr').textContent = Math.max(0, noiseReduction).toFixed(1);
  document.getElementById('stat-L').textContent = '~' + Math.max(1, effLooks).toFixed(1);

  // ── 6. Update detail panels ────────────────────────────────
  updateDetailPanels(winSize, method);
}

// ─────────────────────────────────────────────────────────────
// Event listeners
// ─────────────────────────────────────────────────────────────

/** Attach a click handler to a canvas to select a pixel for inspection. */
function attachCanvasClickHandler(canvasId) {
  const canvas = document.getElementById(canvasId);
  canvas.addEventListener('click', (e) => {
    const rect = canvas.getBoundingClientRect();
    const x = Math.floor((e.clientX - rect.left) / rect.width * N);
    const y = Math.floor((e.clientY - rect.top) / rect.height * N);
    selectedPx = {
      x: Math.max(0, Math.min(N - 1, x)),
      y: Math.max(0, Math.min(N - 1, y)),
    };
    render();
  });
}

// ─────────────────────────────────────────────────────────────
// Initialization
// ─────────────────────────────────────────────────────────────

document.addEventListener('DOMContentLoaded', () => {
  // Sliders
  ['winSize', 'coherence', 'sigFreq'].forEach(id => {
    document.getElementById(id).addEventListener('input', render);
  });

  // Radio buttons
  document.querySelectorAll('input[name="mlMethod"]').forEach(radio => {
    radio.addEventListener('change', render);
  });

  // Canvas click handlers
  attachCanvasClickHandler('canvasRaw');
  attachCanvasClickHandler('canvasML');

  // Draw colorbar
  drawColorbar(document.getElementById('colorbarCanvas'));

  // Initial render
  render();

  // Remove "click to inspect" overlay after first real click
  let hasClicked = false;
  ['canvasRaw', 'canvasML'].forEach(id => {
    document.getElementById(id).addEventListener('click', () => {
      if (!hasClicked) {
        hasClicked = true;
        document.getElementById('overlayRaw').style.opacity = '0';
        document.getElementById('overlayML').style.opacity = '0';
      }
    });
  });
});
