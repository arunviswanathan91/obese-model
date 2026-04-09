# GigaTIME Virtual Spatial Proteomics — CPTAC-PDA

This pipeline applies [GigaTIME](https://github.com/prov-gigatime/GigaTIME) to whole-slide H&E images from the CPTAC pancreatic ductal adenocarcinoma (PDAC) cohort, producing virtual multiplexed immunofluorescence (mIF) maps and comparing immune architecture across BMI groups (Normal / Overweight / Obese).

**Scripts:**

- `nci_idc_data_access.py` — data download from NCI Imaging Data Commons
- `installation_and_qc.py` — GigaTIME setup, model loading, and slide QC
- `gigatime_and_figures.py` — continuous canvas inference, spatial metrics, and figures

---

## Overview

The pipeline runs in three sequential phases:

| Phase | Script | Purpose |
| ----- | ------ | ------- |
| 1 | `nci_idc_data_access.py` | Query and download CPTAC WSI slides from NCI IDC |
| 2 | `installation_and_qc.py` | Install GigaTIME, QC slides, classify usable patients |
| 3 | `gigatime_and_figures.py` | Run inference, compute spatial metrics, generate figures |

---

## Phase 1 — Data Access

Connects to the [NCI Imaging Data Commons](https://datacommons.cancer.gov/repository/imaging-data-commons) via `idc-index` and downloads all CPTAC PDAC whole-slide images (Modality = `SM`) for the `cptac_pda` collection.

**Key steps:**

- SQL query against the IDC index for `SM` modality slides
- Crash-safe batch download with `download_progress.json` resume file
- Folder structure: `wsi_raw/<collection>/<PatientID>/<StudyUID>/<SeriesUID>/`
- Thumbnail visualisation of the first slide to verify the download

**Outputs:**

| File / Folder | Description |
| ------------- | ----------- |
| `wsi_raw/` | Downloaded DICOM series directories |
| `download_progress.json` | Patient-level download checkpoint |

---

## Phase 2 — Installation & QC

Installs GigaTIME from source, downloads model weights from Hugging Face, and performs a multi-stage quality-control census across all ~168 patients before building the final run list.

### Model Configuration

| Setting | Value |
| ------- | ----- |
| Repository | `prov-gigatime/GigaTIME` (HuggingFace) |
| Output channels | 23 (21 protein markers + TRITC + Cy5) |
| Input window | 256 × 256 px |
| Sanity check | `[1, 23, 256, 256]` output tensor |

### Protein Channels

| Excluded (background) | Included in analysis |
| --------------------- | -------------------- |
| TRITC, Cy5 | DAPI, PD-1, CD14, CD4, T-bet, CD34, CD68, CD16, CD11c, CD138, CD20, CD3, CD8, PD-L1, CK, Ki67, Tryptase, Actin-D, Caspase3-D, PHH3-B, Transgelin |

### Tissue Threshold QC

Patches with mean pixel value ≥ 220 are classified as background. The threshold is validated by:

1. Sampling 200 patches from a small region and plotting the pixel-mean distribution
2. Sampling 500 random patches from across the full slide to confirm the bimodal tissue / background separation

### Slide Census

Two census passes are run across all patients:

**Quick census (thumbnail-based):**

- Thumbnail tissue fraction (pixel mean < 220)
- Visual grid of all patient thumbnails saved to `dataset_census_all_patients.png`

**Full census (patch-sampled):**

- 200 random patches per slide → tissue coverage, mean/std pixel intensity, RGB channel means, high-texture fraction (σ > 30)
- Results saved to `dataset_census_full.csv`

### Slide Classification

| Status | Criteria |
| ------ | --------- |
| `GOOD` | tissue fraction ≥ 7% and high-texture fraction ≥ 70% |
| `WARN` | tissue fraction 3–7% or high-texture fraction 50–70% (usable, flagged) |
| `BAD` | tissue fraction < 3% or high-texture fraction < 50% |
| `ERROR` | `wsidicom` failed to open — fallback series attempted |

ERROR patients are recovered by trying all available DICOM series for that patient (largest first, minimum 5,000 × 5,000 px). Successfully fixed patients are reclassified as `WARN`.

**Outputs:**

| File | Description |
| ---- | ----------- |
| `slide_file_stats.csv` | Per-patient file size and series count |
| `dataset_census.csv` | Thumbnail-level tissue fraction |
| `dataset_census_full.csv` | Full patch-sampled QC metrics |
| `dataset_census_classified.csv` | Per-patient GOOD / WARN / BAD / ERROR status |
| `dataset_census_final.csv` | Final census after ERROR recovery |
| `patients_to_run.json` | GOOD + WARN patients — input to Phase 3 |
| `patients_to_skip.json` | BAD + unrecoverable ERROR patients |

---

## Phase 3 — Inference & Figures

### Continuous Canvas Inference

GigaTIME is applied to each slide using a sliding-window strategy with a global logit accumulation canvas:

| Parameter | Value |
| --------- | ----- |
| Extraction window | 256 × 256 px (WSI coordinates) |
| Model input | 512 × 512 px (upsampled via albumentations Resize + ImageNet normalise) |
| Output downsampled back | 256 × 256 px (bilinear interpolation) |
| Stride | 128 px (50% overlap) |
| Background skip | patch mean ≥ 220 → skipped |
| Canvas dtype | float16 × 23 channels |

Logits from overlapping windows are averaged, then thresholded at sigmoid > 0.5 to produce a binary canvas (`uint8`, shape `[23, H, W]`). One WSI pixel = one canvas pixel.

### Fast CK Probe for Section Selection

Large slides often contain multiple tissue fragments. Before running the full canvas (which can require several GB of RAM), a fast CK probe samples 30 random patches from each candidate section (~2 seconds vs ~10 minutes for a full canvas) and returns the mean CK activation density.

**Section selection logic:**

1. Coarse 512 px patch extraction → connected-component clustering → sections
2. Geometric filter: sections with < 50 patches or tissue density < 30% are skipped
3. Fast CK probe: sections with probe CK < 0.01 (no tumour signal) are skipped
4. Best section (highest probe CK) → full canvas inference

### Pipeline Loop

The main inference loop is crash-safe:

- `done_patients.json` records completed patients; the loop resumes after any interruption
- Canvases are saved as compressed `.npz` (`canvas` array + `bbox` in WSI coordinates)
- Per-patient scores (protein activation densities) are appended to `virtual_mif_scores.csv`

**Outputs:**

| File / Folder | Description |
| ------------- | ----------- |
| `canvases/<pid>_canvas.npz` | Binary activation canvas (23, H, W) + bbox |
| `virtual_mif_scores.csv` | Activation density per protein per patient |
| `done_patients.json` | Inference checkpoint |

### Spatial Metrics

Computed from the binary canvases and merged with clinical data:

| Metric | Description |
| ------ | ----------- |
| `stromal_cd8_to_ck_dist` | Mean distance from stromal CD8+ pixels to nearest CK+ (tumour) pixel |
| `transgelin_cd8_dice` | Dice coefficient between Transgelin and CD8 masks |
| `cd11c_to_cd8_dist` | Mean distance from CD11c+ pixels to nearest CD8+ pixel (aDC–T cell axis) |

Results saved to `gigatime_spatial_metrics.csv`.

### BMI-Stratified Figure

The final publication figure compares protein scores and spatial metrics across three BMI groups using violin + jitter plots:

| BMI Group | Threshold |
| --------- | --------- |
| Normal | BMI < 25.0 |
| Overweight | 25.0 ≤ BMI < 30.0 |
| Obese | BMI ≥ 30.0 |

**Panels:**

| Panel | Type | Metric | Y-axis |
| ----- | ---- | ------ | ------ |
| A | Bulk | CD8 activation density | Activation Density |
| B | Bulk | Transgelin activation density | Activation Density |
| C | Spatial | Stromal CD8-to-tumour distance | Distance (pixels) |
| D | Spatial | Transgelin–CD8 Dice coefficient | Dice Coefficient |
| E | Spatial | CD11c-to-CD8 distance | Distance (pixels) |

**Statistics:** Kruskal-Wallis omnibus test + pairwise Mann-Whitney U with Benjamini-Hochberg correction. Significance brackets are drawn only for BH-corrected p < 0.05.

**Output files:** `gigatime_final_figure_v4.pdf` (300 dpi) and `.png` (200 dpi).

---

## Dependencies

```text
idc-index          # NCI IDC data access
wsidicom           # DICOM WSI reading
albumentations     # image augmentation / preprocessing
torch              # GigaTIME inference
huggingface_hub    # model weight download
scikit-learn       # spatial metrics
scipy              # statistics
tqdm               # progress bars
pandas, numpy, matplotlib
```

Install via:

```bash
pip install idc-index wsidicom albumentations huggingface_hub scikit-learn scipy tqdm
apt-get install -y openslide-tools libopenjp2-7
```
