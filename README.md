# Condensate-localized m7G signal immunofluorescence intensity analysis script

Python script for quantifying m7G-cap immunofluorescence within P-bodies and stress granules, implementing the automated local background subtraction and histogram-based peak detection described in the manuscript *mRNA decapping is partitioned by condensed and dilute phase mechanisms*.

## Overview

This script quantifies m7G signals within P-bodies and stress granules using a local background subtraction to account for spatial heterogeneity of m7G fluorescence across the cell. The program takes immunofluoresce micrographs of m7G intensity, ROI lists from imageJ of P-bodys and/or stress granules, and corresponding lists of ROIs that have been dilated. The program then uses the space surrounding each condensate created by the condensate ROIs and dilated ROIs to calculate a locat background for each condensate as described in the methods section. The program will output a .csv file containing background subtracted intesity values for each condensate in the provided ROI lists.

## System Requirements

### Operating Systems

- **macOS**: 10.14 (Mojave) or later (fully tested)
- **Linux**: Ubuntu 18.04+, CentOS 7+, or equivalent distributions
- **Windows**: Windows 10 or later (with Python installed)

### Hardware Requirements

- **Minimum**: 4 GB RAM, 2 CPU cores
- **Recommended**: 8+ GB RAM, 4+ CPU cores (for processing large datasets)
- **Storage**: Sufficient space for input images and output files (typically 2-3x input data size)

### Software Prerequisites

- **Python**: Version 3.8 or higher
- **pip**: Python package installer (included with Python)
- **ImageJ/Fiji**: Required for generating ROI files (not required for running the analysis script)

## Software Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | >=1.20.0 | Array operations and numerical computing |
| pillow | >=8.0.0 | Image file I/O |
| tifffile | >=2021.1.1 | TIFF file handling |
| matplotlib | >=3.4.0 | Visualization and histogram generation |
| scipy | >=1.7.0 | Gaussian filtering and peak detection |
| roifile | >=2021.6.6 | ImageJ ROI file reading |
| scikit-image | >=0.18.0 | Image processing (polygon drawing, contour detection) |

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/marcusjoshm/pb_cap_analysis.git
cd pb_cap_analysis
```

### 2. Create a Virtual Environment (Recommended)

**macOS/Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
```

**Windows:**
```bash
python -m venv venv
venv\Scripts\activate
```

### 3. Install Dependencies

```bash
pip install numpy pillow tifffile matplotlib scipy roifile scikit-image
```

Or install from a requirements file:
```bash
pip install -r requirements.txt
```

**Typical installation time:** 1–3 minutes on a standard desktop computer with a broadband internet connection. Most of this time is spent downloading packages; actual installation is fast.

### 4. Verify Installation

```bash
python -c "import numpy, PIL, tifffile, matplotlib, scipy, roifile, skimage; print('All dependencies installed successfully!')"
```

## Demo Instructions

### Sample Data

Demo data is included in the `m7G_demo_data/` directory. The dataset contains two conditions:

| Condition | Description | P-body ROIs | Stress Granule ROIs |
|-----------|-------------|-------------|---------------------|
| `Untreated/` | Control cells without stress | 78 | None (no SGs present) |
| `1Hr_NaAsO2/` | Cells treated with sodium arsenite for 1 hour | 224 | 306 |

Each condition folder contains:
- Intensity images: `*_Cap_Intensity.tif`, `*_DDX6_Intensity.tif`, and `*_G3BP1_Intensity.tif` (G3BP1 only in stressed condition)
- ROI files: `*_PB_Mask.zip`, `*_PB_Dilated_Mask.zip`, and SG equivalents (stressed condition only)

### Running the Demo

1. Ensure you have completed the installation steps above.

2. From the repository directory, run the analysis on the demo data:
```bash
python analyze_intensity.py m7G_demo_data/
```

   Optionally generate histogram plots for quality control:
```bash
python analyze_intensity.py m7G_demo_data/ --plot-histograms
```

3. Check the output files in each condition subdirectory:
   - `m7G_demo_data/Untreated/Untreated_PB_Cap_intensity_analysis.csv`
   - `m7G_demo_data/Untreated/Untreated_DDX6_intensity_analysis.csv`
   - `m7G_demo_data/1Hr_NaAsO2/1Hr_NaAsO2_PB_Cap_intensity_analysis.csv`
   - `m7G_demo_data/1Hr_NaAsO2/1Hr_NaAsO2_DDX6_intensity_analysis.csv`
   - `m7G_demo_data/1Hr_NaAsO2/1Hr_NaAsO2_SG_Cap_intensity_analysis.csv`
   - `m7G_demo_data/1Hr_NaAsO2/1Hr_NaAsO2_G3BP1_intensity_analysis.csv`

   If you used `--plot-histograms`, corresponding `*_background_histograms.png` files will also be generated.

**Expected run time:** Less than 1 minute for the demo dataset on a standard desktop computer. The demo includes 78 P-body ROIs (Untreated) and 224 P-body + 306 stress granule ROIs (1Hr_NaAsO2).

## Usage

### Basic Usage

```bash
python analyze_intensity.py /path/to/your/data
```

Or use the default path:
```bash
python analyze_intensity.py
```

### Command-Line Flags

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--roi-enlargement` | int | 0 | Pixels to further enlarge dilated ROIs by (binary dilation). Use when the default dilated mask doesn't provide enough background pixels. |
| `--max-background` | float | None | Maximum background threshold for peak detection; only peaks below this value are considered. Use when true background should be below a known value. |
| `--plot-histograms` | flag | off | Generate and save background intensity histogram PNGs for each configuration (see [Histogram plots](#histogram-plots) below). |

**Examples:**

```bash
# Default (no enlargement, no max background constraint, no histogram plots)
python analyze_intensity.py /path/to/data

# Enlarge background ROIs by 5 pixels
python analyze_intensity.py /path/to/data --roi-enlargement 5

# Constrain background peak to values below 100
python analyze_intensity.py /path/to/data --max-background 100

# Generate histogram plots (saved as *_background_histograms.png per config)
python analyze_intensity.py /path/to/data --plot-histograms

# Combine options
python analyze_intensity.py /path/to/data --roi-enlargement 5 --max-background 100 --plot-histograms
```

#### ROI Enlargement (`--roi-enlargement`)
- **Default**: 0 (no enlargement)
- **Purpose**: Expands background ROIs by the specified number of pixels using binary dilation
- **When to use**: When the default dilated mask doesn't provide enough background pixels
- **Effect**: Larger background regions for more robust background estimation

#### Maximum Background Threshold (`--max-background`)
- **Default**: None (no constraint)
- **Purpose**: Constrains background peak selection to only consider peaks below this threshold
- **When to use**: When you have prior knowledge that true background should be below a certain value
- **Effect**: Prevents selection of spurious high-intensity peaks as background

#### Histogram plots
- **Flag**: `--plot-histograms`
- **Default**: Off (no histogram PNGs are generated)
- **Purpose**: When set, saves a background intensity histogram figure for each successfully analyzed configuration
- **Output**: One PNG per config per condition, e.g. `{Condition}_PB_Cap_background_histograms.png`, written into the condition folder. Each figure shows up to 9 ROIs: histogram of background-region intensities, smoothed curve, detected peaks (green dashed lines), and selected background value (blue line).
- **When to use**: For quality control and inspecting peak detection; omit for batch runs when only CSV results are needed.

### Analysis Configurations

The script processes four analysis configurations for each dataset:

| Configuration | Intensity File | ROI Mask | Description |
|---------------|----------------|----------|-------------|
| PB_Cap | `*Cap*Intensity*` | `*PB*Mask*` | P-body cap enrichment |
| DDX6 | `*DDX6*Intensity*` | `*PB*Mask*` | DDX6 protein in P-bodies |
| SG_Cap | `*Cap*Intensity*` | `*SG*Mask*` | Stress granule cap enrichment |
| G3BP1 | `*G3BP1*Intensity*` | `*SG*Mask*` | G3BP1 protein in stress granules |

**Note**: SG configurations are automatically skipped if no SG files are present (e.g., untreated samples without stress granules).

## Input folder: contents, structure, and naming

### What needs to be in the input folder

The input path must be a **base directory** that contains **one subdirectory per condition** (e.g. per treatment or replicate). Do not put files directly in the base directory.

Inside each condition subdirectory you need:

**For P-body (PB) analyses (required for every condition):**
- One **intensity image** for m7G cap: TIFF whose name contains `Cap` and `Intensity`.
- One **intensity image** for DDX6: TIFF whose name contains `DDX6` and `Intensity`.
- One **P-body ROI set**: ImageJ ROI zip whose name contains `PB` and `Mask` but not `Dilated`.
- One **dilated P-body ROI set**: ImageJ ROI zip whose name contains `PB`, `Dilated`, and `Mask`.

**For stress granule (SG) analyses (only when stress granules are present):**
- One **intensity image** for m7G cap: same Cap/Intensity TIFF as above.
- One **intensity image** for G3BP1: TIFF whose name contains `G3BP1` and `Intensity`.
- One **stress granule ROI set**: ImageJ ROI zip whose name contains `SG` and `Mask` but not `Dilated`.
- One **dilated stress granule ROI set**: ImageJ ROI zip whose name contains `SG`, `Dilated`, and `Mask`.

All intensity images must be `.tif`; all ROI sets must be `.zip` (ImageJ ROI Manager export). Hidden/dot files (e.g. `.DS_Store`) are ignored.

### Directory structure

Use a two-level layout: **base directory → one folder per condition**. All condition-specific files live inside that condition’s folder. Example using the test dataset path:

```
/Users/leelab/Desktop/Test_m7G_cap_data/
├── Untreated/
│   ├── Untreated_Cap_Intensity.tif
│   ├── Untreated_DDX6_Intensity.tif
│   ├── Untreated_PB_Mask.zip
│   └── Untreated_PB_Dilated_Mask.zip
└── 1Hr_NaAsO2/
    ├── 1Hr_NaAsO2_Cap_Intensity.tif
    ├── 1Hr_NaAsO2_DDX6_Intensity.tif
    ├── 1Hr_NaAsO2_G3BP1_Intensity.tif
    ├── 1Hr_NaAsO2_PB_Mask.zip
    ├── 1Hr_NaAsO2_PB_Dilated_Mask.zip
    ├── 1Hr_NaAsO2_SG_Mask.zip
    └── 1Hr_NaAsO2_SG_Dilated_Mask.zip
```

- **Untreated**: only P-body files are required; SG files can be omitted and SG analyses will be skipped.
- **1Hr_NaAsO2** (stress condition): include both PB and SG mask/dilated-mask pairs to run all four analyses (PB_Cap, DDX6, SG_Cap, G3BP1).

### File naming

Filenames are matched **case-insensitively** by **keywords**. The condition folder name does not need to match the filename prefix; the script only requires that the right keywords appear in the filename.

**Pattern:** use a consistent prefix (e.g. condition name) plus the keywords below. Examples from the test dataset:

| What you need | Keywords in filename | Example (test data) |
|---------------|----------------------|----------------------|
| m7G cap intensity | `Cap`, `Intensity` | `Untreated_Cap_Intensity.tif` |
| DDX6 intensity | `DDX6`, `Intensity` | `Untreated_DDX6_Intensity.tif` |
| G3BP1 intensity | `G3BP1`, `Intensity` | `1Hr_NaAsO2_G3BP1_Intensity.tif` |
| P-body ROIs (condensate) | `PB`, `Mask` — must **not** contain `Dilated` | `Untreated_PB_Mask.zip` |
| P-body ROIs (dilated, for background) | `PB`, `Dilated`, `Mask` | `Untreated_PB_Dilated_Mask.zip` |
| Stress granule ROIs (condensate) | `SG`, `Mask` — must **not** contain `Dilated` | `1Hr_NaAsO2_SG_Mask.zip` |
| Stress granule ROIs (dilated) | `SG`, `Dilated`, `Mask` | `1Hr_NaAsO2_SG_Dilated_Mask.zip` |

**Rules:**
- Intensity files: extension `.tif`.
- ROI files: ImageJ ROI export as `.zip`.
- Mask and dilated-mask files must be **different** files (script will skip if the same file is used for both).
- For SG analyses you must have both the SG mask zip and the SG dilated mask zip in that condition folder.

### Keyword reference (for custom names)

If you use different naming schemes, the script matches as follows:

| File Type | Required Keywords | Exclude Keywords |
|-----------|-------------------|------------------|
| Cap Intensity | `Cap`, `Intensity` | - |
| DDX6 Intensity | `DDX6`, `Intensity` | - |
| G3BP1 Intensity | `G3BP1`, `Intensity` | - |
| PB Mask | `PB`, `Mask` | `Dilated` |
| PB Dilated Mask | `PB`, `Dilated`, `Mask` | - |
| SG Mask | `SG`, `Mask` | `Dilated` |
| SG Dilated Mask | `SG`, `Dilated`, `Mask` | - |

## Output Files

### Per-Configuration Output

For each successfully analyzed configuration, a CSV file is always written. A histogram PNG is written only when `--plot-histograms` is used.

#### CSV Statistics (`{dataset}_{config}_intensity_analysis.csv`)

| Column | Description |
|--------|-------------|
| `roi_id` | 1-indexed ROI identifier |
| `roi_name` | Original ROI name from ImageJ |
| `n_pixels` | Number of pixels in ROI |
| `n_background_pixels` | Number of pixels in background region |
| `mean_raw` | Mean raw intensity in ROI |
| `median_raw` | Median raw intensity in ROI |
| `background` | Detected background value |
| `mean_bg_subtracted` | Mean background-subtracted intensity |
| `median_bg_subtracted` | Median background-subtracted intensity |
| `background_mean` | Mean of background region intensities |
| `background_std` | Standard deviation of background region |

#### Histogram Visualization (`{dataset}_{config}_background_histograms.png`) — optional

Only generated when you pass `--plot-histograms`.

- Shows background intensity distributions for up to 9 ROIs
- Displays detected peaks (green dashed lines)
- Highlights selected background value (blue solid line)
- Includes smoothed histogram overlay (red line)

## Methods

### Gaussian Peak Detection Algorithm

1. **Histogram Creation**: Background region intensities are binned into 50 bins (range: 0 to max value)

2. **Smoothing**: Gaussian filter (sigma=2) applied to reduce noise

3. **Peak Detection**: scipy's `find_peaks` with 15% prominence threshold

4. **Background Selection**:
   - **Default mode**: Uses the leftmost (lowest intensity) detected peak
   - **Constrained mode** (with `max_background`):
     - Selects most prominent peak below threshold
     - Falls back to highest-frequency bin below threshold if no peaks detected

5. **Fallback**: If no peaks detected, uses the bin with maximum count

### ROI Enlargement

- Uses binary dilation with 4-connected structuring element
- Equivalent to ImageJ's RoiEnlarger functionality
- Contours are extracted from enlarged masks to create new ROI coordinates

### Background Subtraction

For each ROI:
1. Extract intensity values within ROI mask
2. Extract intensity values within dilated (background) mask
3. Detect background value using Gaussian peak detection
4. Subtract background from ROI intensities
5. Calculate statistics (mean, median) for raw and background-subtracted values

## Troubleshooting

### Common Issues

**"No subdirectories found"**
- Ensure the base directory contains subdirectories with data files
- Check that subdirectories are not hidden (don't start with `.`)

**"Skipping [config]: No intensity file found"**
- Verify intensity files contain the required keywords (e.g., `Cap`, `Intensity`)
- Check file extensions are `.tif`

**"Skipping [config]: No mask file found"**
- Ensure ROI zip files are present with correct naming
- Verify files contain required keywords (e.g., `PB`, `Mask`)

**Import errors**
- Reinstall dependencies: `pip install --upgrade numpy pillow tifffile matplotlib scipy roifile scikit-image`
- Verify Python version: `python --version` (must be 3.8+)

## License

MIT License