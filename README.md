# Microscopy Intensity Analysis

Universal Python analysis pipeline for analyzing microscopy intensity data using per-ROI local background subtraction with Gaussian peak detection.

## Overview

This pipeline analyzes microscopy images to measure intensity in regions of interest (ROIs) with local background subtraction. It uses per-ROI background subtraction based on dilated region intensity analysis with Gaussian peak detection. Originally developed for P-body (PB) and Stress Granule (SG) m7G Cap enrichment analysis, it has been refactored to work with any type of microscopy intensity measurement.

The pipeline gracefully handles cases where some analysis configurations may not have all required files (e.g., SG analysis when no stress granules are present in untreated samples).

## Features

- **Multi-Configuration Analysis**: Simultaneously analyzes multiple intensity/ROI combinations:
  - PB_Cap: P-body cap enrichment
  - DDX6: P-body marker protein intensity
  - SG_Cap: Stress granule cap enrichment
  - G3BP1: Stress granule marker protein intensity
- **Gaussian Peak Detection**: Detects peaks in background region histograms to identify true background intensity
- **Per-ROI Local Background Subtraction**: Each ROI is analyzed with its own local background region
- **Flexible File Matching**: Keyword-based file matching works with various naming conventions
- **Graceful Handling of Missing Files**: Automatically skips configurations when required files aren't present
- **Dot File Filtering**: Ignores macOS hidden files (`.DS_Store`, `._*` files)
- **Interactive Parameter Configuration**: Prompts for per-configuration analysis parameters
- **Visualization**: Generates background intensity histograms showing detected peaks
- **CSV Export**: Detailed per-ROI statistics with all measurements

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

Or install from a requirements file (if provided):
```bash
pip install -r requirements.txt
```

### 4. Verify Installation

```bash
python -c "import numpy, PIL, tifffile, matplotlib, scipy, roifile, skimage; print('All dependencies installed successfully!')"
```

## Demo Instructions

### Sample Data

Sample data will be provided in a future release. The sample dataset will include:
- Pre-processed intensity images (`.tif` files)
- ImageJ ROI files (`.zip` files) for P-bodies and stress granules
- Expected output files for validation

### Running the Demo

Once sample data is available:

1. Download and extract the sample data to a directory (e.g., `~/sample_data/`)

2. Run the analysis:
```bash
python analyze_intensity.py ~/sample_data/
```

3. Follow the interactive prompts to configure analysis parameters (or press Enter to accept defaults)

4. Check the output files in each condition subdirectory:
   - `*_intensity_analysis.csv` - Per-ROI statistics
   - `*_background_histograms.png` - Visualization of detected peaks

## Usage

### Basic Usage

```bash
python analyze_intensity.py /path/to/your/data
```

Or use the default path:
```bash
python analyze_intensity.py
```

### Interactive Configuration

When you run the script, you'll be prompted to configure parameters for each analysis type:

```
============================================================
CONFIGURATION PARAMETERS
============================================================
Please enter parameters for each analysis configuration.
Press Enter to use default values shown in brackets.

--- PB_Cap Configuration ---
  ROI enlargement pixels [0]:
  Maximum background threshold [None]:

--- DDX6 Configuration ---
  ROI enlargement pixels [0]:
  Maximum background threshold [None]: 100

--- SG_Cap Configuration ---
  ROI enlargement pixels [0]: 5
  Maximum background threshold [None]:

--- G3BP1 Configuration ---
  ROI enlargement pixels [0]: 5
  Maximum background threshold [None]: 150
```

### Parameter Descriptions

#### ROI Enlargement (`enlarge_rois`)
- **Default**: 0 (no enlargement)
- **Purpose**: Expands background ROIs by the specified number of pixels using binary dilation
- **When to use**: When the default dilated mask doesn't provide enough background pixels
- **Effect**: Larger background regions for more robust background estimation

#### Maximum Background Threshold (`max_background`)
- **Default**: None (no constraint)
- **Purpose**: Constrains background peak selection to only consider peaks below this threshold
- **When to use**: When you have prior knowledge that true background should be below a certain value
- **Effect**: Prevents selection of spurious high-intensity peaks as background

### Analysis Configurations

The script processes four analysis configurations for each dataset:

| Configuration | Intensity File | ROI Mask | Description |
|---------------|----------------|----------|-------------|
| PB_Cap | `*Cap*Intensity*` | `*PB*Mask*` | P-body cap enrichment |
| DDX6 | `*DDX6*Intensity*` | `*PB*Mask*` | DDX6 protein in P-bodies |
| SG_Cap | `*Cap*Intensity*` | `*SG*Mask*` | Stress granule cap enrichment |
| G3BP1 | `*G3BP1*Intensity*` | `*SG*Mask*` | G3BP1 protein in stress granules |

**Note**: SG configurations are automatically skipped if no SG files are present (e.g., untreated samples without stress granules).

## Input Data Structure

### Expected Directory Structure

```
/path/to/data/
├── Condition1/
│   ├── Condition1_Cap_Intensity.tif       # Cap channel intensity image
│   ├── Condition1_DDX6_Intensity.tif      # DDX6 channel intensity image
│   ├── Condition1_G3BP1_Intensity.tif     # G3BP1 channel intensity (if SG present)
│   ├── Condition1_PB_Mask.zip             # P-body ROIs
│   ├── Condition1_PB_Dilated_Mask.zip     # Dilated P-body ROIs (background)
│   ├── Condition1_SG_Mask.zip             # Stress granule ROIs (if SG present)
│   └── Condition1_SG_Dilated_Mask.zip     # Dilated SG ROIs (if SG present)
├── Condition2/
│   └── (same structure)
└── Condition3/
    └── (same structure, SG files optional)
```

### File Naming Requirements

Files are matched using case-insensitive keyword matching:

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

For each successfully analyzed configuration, two files are generated:

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

**Note**: Negative values (when ROI intensity < background) are replaced with empty cells in the CSV.

#### Histogram Visualization (`{dataset}_{config}_background_histograms.png`)

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

### Performance Tips

- Process smaller batches if memory is limited
- Close other applications when processing large datasets
- Use SSD storage for faster file I/O

## ImageJ Macros

The `Macros/` directory contains ImageJ macros for preprocessing:

- **PB_Macro.ijm**: P-body analysis with stress granule subtraction
- **SG_Macro.ijm**: Stress granule ROI extraction
- **copy_masks_and_raw.sh**: Organizes raw data into analysis directory structure
- **copy_ch0_to_processed.sh**: Copies Cap intensity files to processed directories

These macros prepare data for the Python analysis pipeline.

## Workflow Integration

This script is designed to work with outputs from:
1. **PerCell** analysis software (combined_masks and raw_data directories)
2. **ImageJ/Fiji** ROI Manager exports
3. Custom preprocessing pipelines that generate compatible file structures

## License

MIT License

## Author

Joshua Marcus

## Citation

If you use this software in your research, please cite:
```
Marcus, J. (2024). Microscopy Intensity Analysis Pipeline.
https://github.com/marcusjoshm/pb_cap_analysis
```

## Contributing

Contributions are welcome! Please submit issues and pull requests to the GitHub repository.

## Changelog

### v2.0.0 (2025-01)
- Removed background multiplication factor (`bg_factor`) parameter
- Added dot file filtering for macOS compatibility
- Improved handling of missing files (graceful skipping)
- Updated to match PerCell plugin implementation
- Enhanced documentation and error messages

### v1.0.0 (2024-10)
- Initial release with per-ROI background subtraction
- Gaussian peak detection for background estimation
- Multi-configuration analysis support
