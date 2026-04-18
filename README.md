# InSAR Multilooking Visualizer

An interactive browser-based tool for understanding **complex multilooking** in Interferometric Synthetic Aperture Radar (InSAR) processing.

Live demo: open `index.html` directly in any modern browser — **no server, no build step, no dependencies**.

---

![screenshot placeholder](https://via.placeholder.com/900x480/0d0f12/4f8ef7?text=InSAR+Multilooking+Visualizer)

---

## Background

This visualizer is based on concepts from:

> S. Samiei-Esfahany, *"Exploitation of Distributed Scatterers in Synthetic Aperture Radar Interferometry"*, PhD Thesis, Delft University of Technology, 2017.
> DOI: [10.4233/uuid:22d46f1e-9061-46b0-9726-760c41404b6f](https://doi.org/10.4233/uuid:22d46f1e-9061-46b0-9726-760c41404b6f)

### What is multilooking?

In InSAR, each pixel's interferometric phase is the sum of a deformation signal and noise. For **distributed scatterers (DS)** — natural land surfaces such as farmland, deserts, or bare soil — the noise is high (low coherence). Multilooking spatially averages the complex interferogram over *L* neighboring pixels to reduce phase noise by a factor of ~√L:

```
φ̂_k = arg[ (1/L) · Σ P_M,i · P_S,i* ]
```

This is the **Maximum Likelihood estimator** of the interferometric phase for distributed scatterers (Rodriguez & Martin, 1992).

The **Cramér-Rao lower bound** for phase variance after L looks is:

```
σ_φ² ≥ (1 − |γ|²) / (2L|γ|²)
```

where |γ| ∈ [0, 1] is the coherence.

---

## Features

| Feature | Description |
|---|---|
| **Synthetic SAR scene** | Generates a realistic N×N interferogram with configurable coherence and signal frequency |
| **Three multilooking methods** | Boxcar (uniform), Gaussian (weighted), Adaptive (brotherhood) |
| **Live statistics** | Raw vs. multilooked phase std dev, noise reduction %, effective looks L |
| **Pixel inspection** | Click any pixel to see its window, phasor diagram, and phase histogram |
| **Complex phasor diagram** | Visualizes individual phasors vs. resultant vector — the core of complex averaging |
| **Brotherhood mask** | Adaptive mode shows which neighbors pass the statistical homogeneity test |
| **Theory panel** | Equations and context from the thesis rendered inline |

---

## Multilooking Methods Implemented

### 1. Boxcar (uniform)
Standard rectangular window — all pixels in the N×N window get equal weight. Used in original SBAS (Berardino et al., 2002).

### 2. Gaussian (weighted)
Center pixel gets highest weight; weights decay as `exp(-r²/2σ²)`. Reduces edge artifacts at the cost of slightly reduced effective looks.

### 3. Adaptive (SqueeSAR brotherhood)
Based on **SqueeSAR** (Ferretti et al., 2011b) and the parametric mean test (Jiang et al., 2014a).

For each pixel *p*, neighboring pixels *q* are tested for statistical homogeneity:

```
H0 : E{Ā_q} = μ_Ap    (q is a "brother" of p)
T  = √N · (Ā_q − μ̂_Ap) / σ̂_Ap
Reject if |T| > 1.96   (α = 5%)
```

Only **statistically homogeneous pixels (SHPs)** — the "brotherhood" — are included in the complex average. This prevents mixing heterogeneous terrain types (e.g., a road pixel mixed with a field pixel).

---

## File Structure

```
insar-multilooking/
├── index.html       # Main HTML — layout, panels, controls
├── style.css        # Dark-theme styling
├── main.js          # All physics, algorithms, and rendering
└── README.md        # This file
```

**`main.js`** is fully documented with JSDoc. Key functions:

| Function | Purpose |
|---|---|
| `generateSignalPhase(freq)` | Sinusoidal deformation signal model |
| `generateNoise(coherence)` | Complex circular Gaussian noise (Goodman, 1976) |
| `boxcarMultilook(...)` | Uniform + Gaussian weighted multilooking |
| `makeGaussianKernel(size)` | Normalized 2D Gaussian kernel |
| `adaptiveMultilook(...)` | Brotherhood selection + adaptive averaging |
| `drawPhaseImage(...)` | Phase → HSL colormap rendering |
| `drawPhasorDiagram(...)` | Complex phasor vector display |
| `drawPhaseHistogram(...)` | Raw vs. ML phase distribution |
| `circularStdDev(arr)` | Circular (wrapped) standard deviation |

---

## Usage

1. Clone or download the repository
2. Open `index.html` in any modern browser (Chrome, Firefox, Safari, Edge)
3. Use the sidebar controls to adjust parameters
4. Click anywhere on the phase images to inspect that pixel in detail

No internet connection required after initial font load.

---

## Physics Notes

### SAR SLC Phase Model
```
ψ = W{ ψ_range + ψ_atmo + ψ_scat + ψ_noise }
```

### Interferometric Phase Decomposition
```
φ_MS = W{ ψ_M − ψ_S }
     = φ_flat + φ_topo + φ_defo + φ_atmo + φ_scat + φ_noise
```

### Deformation-Phase Relationship
```
φ_defo = −(4π/λ) · D_LOS
```
where λ ≈ 5.6 cm (C-band Sentinel-1). A 1 mm displacement gives ~0.11 rad of phase change.

### Coherence Sources
| Factor | Symbol | Cause |
|---|---|---|
| Temporal decorrelation | γ_T | Changes in surface scattering over time |
| Baseline decorrelation | γ_geom | Different incidence angles → different phasor superposition |
| Doppler decorrelation | γ_dc | Different squint angles |
| Thermal noise | γ_thermal | System SNR |
| Processing noise | γ_proc | Coregistration and resampling errors |

Total: `γ_total = γ_T · γ_geom · γ_dc · γ_thermal · γ_proc`

---

## References

- Samiei-Esfahany, S. (2017). *Exploitation of Distributed Scatterers in Synthetic Aperture Radar Interferometry*. PhD Thesis, TU Delft.
- Ferretti, A., Fumagalli, A., Novali, F., Prati, C., Rocca, F., & Rucci, A. (2011b). A new algorithm for processing interferometric data-stacks: SqueeSAR. *IEEE TGRS*, 49(9), 3460–3470.
- Berardino, P., Fornaro, G., Lanari, R., & Sansosti, E. (2002). A new algorithm for surface deformation monitoring based on small baseline differential SAR interferograms. *IEEE TGRS*, 40(11), 2375–2383.
- Rodriguez, E., & Martin, J. M. (1992). Theory and design of interferometric synthetic aperture radars. *IEE Proceedings F*, 139(2), 147–159.
- Jiang, M., Ding, X., Hanssen, R. F., Malhotra, R., & Chang, L. (2014a). Fast statistically homogeneous pixel selection for covariance matrix estimation for multitemporal InSAR. *IEEE TGRS*, 53(3), 1213–1224.
- Goodman, J. W. (1976). Some fundamental properties of speckle. *JOSA*, 66(11), 1145–1150.

---

## License

MIT License — free to use, modify, and distribute.

---
