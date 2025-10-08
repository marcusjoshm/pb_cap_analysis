# P-Body Cap Enrichment Analysis

Python analysis pipeline for analyzing microscopy data to assess Cap enrichment in P-bodies (processing bodies).

## Overview

This pipeline analyzes microscopy images to determine if P-bodies are enriched with Cap signal compared to background. It uses per-particle background subtraction based on perimeter intensity analysis with Gaussian peak detection.

## Features

- **Data Extraction**: Extracts Cap and G3BP1 intensity values from mask regions and perimeters
- **Background Subtraction**: Two methods available:
  - `minimum`: Uses minimum perimeter intensity
  - `gaussian_peaks`: Detects Gaussian peaks in perimeter histograms to identify true background (handles particles near Cap-enriched structures)
- **Per-Particle Analysis**: Each P-body is analyzed individually with its own perimeter-based background
- **ROI Support**: Uses ImageJ ROI zip files for precise particle definitions
- **Visualization**: Generates intensity maps and perimeter histograms

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
├── Untreated/
│   ├── *Mask*.tif                  # P-body mask
│   ├── *Dilated*Mask*.tif          # Dilated mask (for perimeter)
│   ├── *Cap*Intensity*.tif         # Cap intensity image
│   ├── *G3BP1*Intensity*.tif       # G3BP1 intensity image
│   ├── *Mask*.zip                  # ImageJ ROIs for particles
│   └── *Perimeter*.zip             # ImageJ ROIs for perimeters
└── As Treated/
    └── (same structure)
```

## Output Files

### Data Extraction
- `{dataset}_Cap_Full.png` - Full Cap intensity visualization
- `{dataset}_Cap_Perimeter.png` - Cap intensity at perimeter only
- `{dataset}_G3BP1_Full.png` - Full G3BP1 intensity visualization
- `{dataset}_G3BP1_Perimeter.png` - G3BP1 intensity at perimeter only
- `*Perimeter Mask.tif` - Binary perimeter masks (8-bit)

### Enrichment Analysis
- `{dataset}_enrichment_analysis.csv` - Per-particle statistics including:
  - Particle ID and ROI name
  - Raw and background-subtracted intensities
  - Background values
  - Perimeter statistics
- `{dataset}_perimeter_histograms.png` - Histogram visualizations showing detected peaks

## Methods

### Background Subtraction Approaches

1. **Gaussian Peak Detection** (default):
   - Analyzes perimeter intensity histogram
   - Detects Gaussian peaks
   - Single peak → uses that peak as background
   - Multiple peaks → uses lower peak as background (assumes higher peak is from nearby enriched structure)

2. **Minimum Value**:
   - Uses minimum perimeter intensity as background
   - Simpler but may overestimate enrichment if minimum is anomalously low

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
