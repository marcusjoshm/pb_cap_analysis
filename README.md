# Microscopy Intensity Enrichment Analysis

Universal Python analysis pipeline for analyzing microscopy data to assess intensity enrichment in particles using per-particle background subtraction.

## Overview

This pipeline analyzes microscopy images to determine if particles are enriched with signal compared to background. It uses per-particle background subtraction based on perimeter intensity analysis with Gaussian peak detection. Originally developed for P-body Cap enrichment analysis, it has been refactored to work with any type of microscopy intensity data.

## Features

- **Universal Intensity Analysis**: Works with any intensity image - no longer limited to specific markers (Cap, G3BP1, etc.)
- **Background Subtraction**: Two methods available:
  - `minimum`: Uses minimum perimeter intensity
  - `gaussian_peaks`: Detects Gaussian peaks in perimeter histograms to identify true background (handles particles near enriched structures)
- **Per-Particle Analysis**: Each particle is analyzed individually with its own perimeter-based background
- **ROI Support**: Uses ImageJ ROI zip files for precise particle definitions
- **Visualization**: Generates intensity maps and perimeter histograms
- **Flexible File Naming**: Automatically detects any `.tif` file with "intensity" in the filename

## Setup

### Prerequisites

- Python 3.8 or higher
- pip

### Installation

1. Clone the repository:
```bash
git clone https://github.com/marcusjoshm/pb_cap_analysis.git
cd pb_cap_analysis
```

2. Create a virtual environment:
```bash
python3 -m venv venv
```

3. Activate the virtual environment:
```bash
source venv/bin/activate
```

4. Install dependencies:
```bash
pip install numpy pillow tifffile matplotlib scipy roifile scikit-image
```

## Usage

### 1. Extract Microscopy Data

Extract intensity values and create perimeter masks:

```bash
python extract_microscopy_data.py
```

This script:
- Loads mask and intensity images
- Creates perimeter masks (donut-shaped regions)
- Extracts intensity arrays for mask and perimeter regions
- Generates visualization plots with viridis colormap
- Saves perimeter masks as 8-bit TIFF files (values: 0 and 255)

### 2. Analyze Enrichment

Perform per-particle background subtraction and enrichment analysis:

```bash
python analyze_enrichment.py /path/to/your/data
```

Or use the default path (if not modified):

```bash
python analyze_enrichment.py
```

**Command-line options:**

- No arguments: Uses default path `/Volumes/NX-01-A/2025-10-08_test_data`
- One argument: Custom base directory path containing subdirectories with data
- `--bg-factor`: Multiplication factor for background subtraction (default: 1.0)
- `--max-background`: Maximum background value threshold. Only peaks below this value will be considered for background (default: None, no constraint)
- `--enlarge-rois`: Number of pixels to enlarge dilated ROIs using binary dilation (default: 0, no enlargement)

**Examples with command-line options:**

```bash
# Apply 0.8x multiplication to background values (more conservative background subtraction)
python analyze_enrichment.py /path/to/your/data --bg-factor 0.8

# Apply 1.2x multiplication to background values (more aggressive background subtraction)
python analyze_enrichment.py /path/to/your/data --bg-factor 1.2

# Force background peak selection below 100 (useful when you know true background is low)
python analyze_enrichment.py /path/to/your/data --max-background 100

# Enlarge perimeter ROIs by 5 pixels to capture more background
python analyze_enrichment.py /path/to/your/data --enlarge-rois 5

# Combine multiple options
python analyze_enrichment.py /path/to/your/data --bg-factor 1.0 --max-background 100 --enlarge-rois 5

# Default behavior (1.0x factor, no constraint, no enlargement)
python analyze_enrichment.py /path/to/your/data
```

**Parameter descriptions:**

- **`--bg-factor`**: Adjusts background subtraction stringency. Values < 1.0 result in more conservative background estimates (less background subtracted), while values > 1.0 result in more aggressive background subtraction.

- **`--max-background`**: Constrains background peak selection to only consider peaks below this threshold. When set:
  - Selects the most prominent peak below the threshold
  - If no peaks detected below threshold, uses the histogram bin with highest frequency below threshold
  - Useful for datasets where you know the true background should be low (e.g., < 100)
  - When not set (default), uses the leftmost (lowest intensity) detected peak

- **`--enlarge-rois`**: Expands the perimeter ROIs by the specified number of pixels using binary dilation. This increases the region used for background estimation, which can be helpful when the default perimeter is too narrow.

This script:

- Automatically scans all subdirectories in the base directory
- Loads ImageJ ROI files for particles and perimeters
- Analyzes each particle individually
- Detects Gaussian peaks in perimeter intensity histograms
- Performs background subtraction using appropriate peak values
- Exports CSV files with detailed statistics
- Creates histogram visualizations for each particle
- Skips subdirectories missing required files with informative messages

## Input Data Structure

Expected directory structure:
```
/path/to/data/
├── Condition1/
│   ├── *Mask*.tif                  # Particle mask
│   ├── *Dilated*Mask*.tif          # Dilated mask (for perimeter, optional)
│   ├── *Intensity*.tif             # ANY intensity image (e.g., Cap, GFP, mCherry, etc.)
│   ├── *Mask*.zip                  # ImageJ ROIs for particles
│   └── *Dilated*.zip               # ImageJ ROIs for perimeters (dilated regions)
└── Condition2/
    └── (same structure)
```

**Key Requirements:**
- Intensity images must contain "intensity" in the filename (case-insensitive)
- The script now accepts **any** `.tif` file with "intensity" in the name
- No specific marker name (like "Cap" or "G3BP1") is required in the filename
- Examples of valid intensity filenames:
  - `sample_intensity.tif`
  - `GFP_Intensity.tif`
  - `Cap_intensity_channel1.tif`
  - `intensity_map.tif`

## Output Files

### Data Extraction
- `{dataset}_Cap_Full.png` - Full Cap intensity visualization (if using Cap data)
- `{dataset}_Cap_Perimeter.png` - Cap intensity at perimeter only (if using Cap data)
- `{dataset}_G3BP1_Full.png` - Full G3BP1 intensity visualization (if using G3BP1 data)
- `{dataset}_G3BP1_Perimeter.png` - G3BP1 intensity at perimeter only (if using G3BP1 data)
- `*Perimeter Mask.tif` - Binary perimeter masks (8-bit)

### Enrichment Analysis
- `{dataset}_enrichment_analysis.csv` - Per-particle statistics including:
  - Particle ID and ROI name
  - Raw and background-subtracted intensities (columns now prefixed with `image_` instead of `cap_`)
  - Background values
  - Perimeter statistics
  - **CSV Column Names** (updated to be marker-agnostic):
    - `image_mean_raw`, `image_median_raw`
    - `image_background`
    - `image_mean_bg_subtracted`, `image_median_bg_subtracted`
    - `image_perimeter_mean`, `image_perimeter_std`
  - **Note**: Negative values (from background subtraction) are replaced with empty cells in the CSV
- `{dataset}_perimeter_histograms.png` - Histogram visualizations showing detected peaks

## Methods

### Background Subtraction Approaches

1. **Gaussian Peak Detection** (default):
   - Analyzes perimeter intensity histogram with range starting at 0 to capture near-zero peaks
   - Uses 15% of maximum histogram count as minimum prominence threshold to detect significant peaks
   - Detects peaks in smoothed histogram (Gaussian filter with sigma=2)
   - **Default mode** (no `--max-background`): Uses the leftmost (lowest intensity) detected peak as background
   - **Constrained mode** (with `--max-background`): Only considers peaks below the specified threshold
     - Selects the most prominent peak below the threshold
     - If no peaks detected, uses the histogram bin with highest frequency below the threshold
     - Ensures background values stay within expected range
   - Multiple peaks → assumes higher peak is from nearby enriched structure
   - Fallback: If no peaks detected at all, uses the bin with maximum count as background

2. **Minimum Value**:
   - Uses minimum perimeter intensity as background
   - Simpler but may overestimate enrichment if minimum is anomalously low

### ROI Enlargement

The `--enlarge-rois` parameter allows expanding perimeter ROIs for background estimation:
- Uses binary dilation with 4-connected structuring element (equivalent to ImageJ's RoiEnlarger)
- Expands the perimeter region by the specified number of pixels
- Useful when the default dilated mask doesn't provide enough perimeter pixels for reliable background estimation
- Example: `--enlarge-rois 5` expands each perimeter ROI by 5 pixels in all directions

### Background Multiplication Factor

The `--bg-factor` parameter allows fine-tuning of background subtraction:
- **Values < 1.0** (e.g., 0.8): More conservative - reduces the background value before subtraction, resulting in lower enrichment scores
- **Values = 1.0**: Default - uses the detected background value as-is
- **Values > 1.0** (e.g., 1.2): More aggressive - increases the background value before subtraction, resulting in higher enrichment scores

This is useful when you want to adjust for systematic over- or under-estimation of background.

### Maximum Background Constraint

The `--max-background` parameter constrains background peak selection:
- Ensures background values don't exceed a known threshold for your dataset
- Useful when you have prior knowledge about expected background intensity ranges
- Example: If you know true background is around 20 and never exceeds 100, use `--max-background 100`
- Helps avoid selecting spurious high-intensity peaks as background
- When used, the algorithm prioritizes peaks below the threshold based on prominence
- Fallback uses the highest frequency histogram bin below the threshold if no peaks are detected

### Negative Value Handling

After background subtraction, some particles may have negative enrichment values (particle intensity < background). In the output CSV:
- **Negative values are replaced with empty cells**
- This prevents misleading negative enrichment scores
- Empty cells make it easy to identify particles where signal is below background
- All other numeric and non-numeric values are preserved as-is

## Dependencies

- numpy - Array operations
- pillow - Image I/O
- tifffile - TIFF file handling
- matplotlib - Visualization
- scipy - Peak detection and signal processing
- roifile - ImageJ ROI file reading
- scikit-image - Image processing (polygon drawing)

## License

MIT License

## Author

Joshua Marcus
