#!/usr/bin/env python3
"""
P-body Cap Enrichment Analysis
Analyzes Cap enrichment in P-bodies using per-particle background subtraction.
"""

import numpy as np
from pathlib import Path
from PIL import Image
import tifffile
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from roifile import ImagejRoi
import zipfile


def load_image(filepath):
    """Load an image file and convert to numpy array."""
    img = Image.open(filepath)
    return np.array(img)


def load_imagej_rois(roi_zip_path):
    """
    Load ROIs from ImageJ ROI zip file.

    Args:
        roi_zip_path: Path to ImageJ ROI zip file

    Returns:
        List of ROI dictionaries with coordinates and properties
    """
    rois = []
    with zipfile.ZipFile(roi_zip_path, 'r') as zf:
        for name in sorted(zf.namelist()):  # Sort to ensure consistent ordering
            with zf.open(name) as roi_file:
                roi_bytes = roi_file.read()
                roi = ImagejRoi.frombytes(roi_bytes)
                rois.append({
                    'name': name,
                    'roi': roi,
                    'coordinates': roi.coordinates()
                })
    return rois


def roi_to_mask(roi_coords, image_shape):
    """
    Convert ROI coordinates to a binary mask.

    Args:
        roi_coords: ROI coordinates array (N x 2)
        image_shape: Shape of the output mask

    Returns:
        Binary mask with ROI pixels set to 1
    """
    from skimage.draw import polygon

    mask = np.zeros(image_shape, dtype=np.uint8)
    if len(roi_coords) > 0:
        # ROI coordinates are (x, y), but polygon needs (row, col) = (y, x)
        rr, cc = polygon(roi_coords[:, 1], roi_coords[:, 0], image_shape)
        mask[rr, cc] = 1
    return mask


def gaussian(x, amp, mean, std):
    """Gaussian function for peak fitting."""
    return amp * np.exp(-((x - mean) ** 2) / (2 * std ** 2))


def find_gaussian_peaks(data, n_bins=50):
    """
    Find Gaussian peaks in a histogram of the data.

    Args:
        data: Array of intensity values
        n_bins: Number of bins for histogram

    Returns:
        Dictionary with peak information and histogram data
    """
    if len(data) == 0:
        return None

    # Create histogram
    hist, bin_edges = np.histogram(data, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Smooth the histogram slightly to reduce noise
    from scipy.ndimage import gaussian_filter1d
    hist_smooth = gaussian_filter1d(hist.astype(float), sigma=1)

    # Find peaks in the histogram
    peaks, properties = find_peaks(hist_smooth, prominence=np.max(hist_smooth) * 0.1)

    if len(peaks) == 0:
        # No clear peaks found, use mean as fallback
        return {
            'n_peaks': 0,
            'peaks': [],
            'hist': hist,
            'bin_centers': bin_centers,
            'background_value': np.mean(data)
        }

    # Sort peaks by prominence (highest first)
    peak_prominences = properties['prominences']
    sorted_indices = np.argsort(peak_prominences)[::-1]
    peaks = peaks[sorted_indices]

    # Get the peak positions (intensity values)
    peak_positions = bin_centers[peaks]

    # Determine background value
    if len(peaks) >= 2:
        # Two or more peaks: use the lower peak value
        background_value = min(peak_positions[0], peak_positions[1])
    elif len(peaks) == 1:
        # Single peak: use that peak value
        background_value = peak_positions[0]
    else:
        background_value = np.mean(data)

    return {
        'n_peaks': len(peaks),
        'peaks': peak_positions[:2] if len(peaks) >= 2 else peak_positions,
        'hist': hist,
        'bin_centers': bin_centers,
        'hist_smooth': hist_smooth,
        'peak_indices': peaks[:2] if len(peaks) >= 2 else peaks,
        'background_value': background_value
    }


def analyze_particle_enrichment_from_rois(mask_rois, perimeter_rois, cap_intensity,
                                          g3bp1_intensity, image_shape, method='gaussian_peaks'):
    """
    Analyze Cap and G3BP1 enrichment for each particle using ROI definitions.

    Args:
        mask_rois: List of ROI dictionaries for particles
        perimeter_rois: List of ROI dictionaries for perimeters
        cap_intensity: Cap intensity image
        g3bp1_intensity: G3BP1 intensity image
        image_shape: Shape of the images
        method: 'minimum' or 'gaussian_peaks'

    Returns:
        Dictionary with per-particle analysis results
    """
    n_particles = len(mask_rois)
    print(f"  Found {n_particles} particles from ROI files")

    results = []

    for i in range(n_particles):
        # Create masks from ROIs
        particle_mask = roi_to_mask(mask_rois[i]['coordinates'], image_shape)
        perimeter_mask = roi_to_mask(perimeter_rois[i]['coordinates'], image_shape)

        # Extract intensities
        cap_particle_intensities = cap_intensity[particle_mask > 0]
        cap_perimeter_intensities = cap_intensity[perimeter_mask > 0]

        g3bp1_particle_intensities = g3bp1_intensity[particle_mask > 0]
        g3bp1_perimeter_intensities = g3bp1_intensity[perimeter_mask > 0]

        # Skip if no pixels found
        if len(cap_particle_intensities) == 0 or len(cap_perimeter_intensities) == 0:
            continue

        # Determine background value based on method
        if method == 'minimum':
            cap_background = np.min(cap_perimeter_intensities)
            g3bp1_background = np.min(g3bp1_perimeter_intensities)
            cap_peak_info = None
            g3bp1_peak_info = None
        elif method == 'gaussian_peaks':
            cap_peak_info = find_gaussian_peaks(cap_perimeter_intensities)
            g3bp1_peak_info = find_gaussian_peaks(g3bp1_perimeter_intensities)
            cap_background = cap_peak_info['background_value'] if cap_peak_info else np.mean(cap_perimeter_intensities)
            g3bp1_background = g3bp1_peak_info['background_value'] if g3bp1_peak_info else np.mean(g3bp1_perimeter_intensities)
        else:
            raise ValueError(f"Unknown method: {method}")

        # Background-subtracted intensities
        cap_bg_subtracted = cap_particle_intensities - cap_background
        g3bp1_bg_subtracted = g3bp1_particle_intensities - g3bp1_background

        # Calculate statistics
        particle_result = {
            'particle_id': i + 1,  # 1-indexed particle ID
            'roi_name': mask_rois[i]['name'],
            'n_pixels': len(cap_particle_intensities),
            'n_perimeter_pixels': len(cap_perimeter_intensities),

            # Cap statistics
            'cap_mean_raw': np.mean(cap_particle_intensities),
            'cap_median_raw': np.median(cap_particle_intensities),
            'cap_background': cap_background,
            'cap_mean_bg_subtracted': np.mean(cap_bg_subtracted),
            'cap_median_bg_subtracted': np.median(cap_bg_subtracted),
            'cap_perimeter_mean': np.mean(cap_perimeter_intensities),
            'cap_perimeter_std': np.std(cap_perimeter_intensities),
            'cap_peak_info': cap_peak_info,

            # G3BP1 statistics
            'g3bp1_mean_raw': np.mean(g3bp1_particle_intensities),
            'g3bp1_median_raw': np.median(g3bp1_particle_intensities),
            'g3bp1_background': g3bp1_background,
            'g3bp1_mean_bg_subtracted': np.mean(g3bp1_bg_subtracted),
            'g3bp1_median_bg_subtracted': np.median(g3bp1_bg_subtracted),
            'g3bp1_perimeter_mean': np.mean(g3bp1_perimeter_intensities),
            'g3bp1_perimeter_std': np.std(g3bp1_perimeter_intensities),
            'g3bp1_peak_info': g3bp1_peak_info,
        }

        results.append(particle_result)

    return {
        'particles': results,
        'n_particles': n_particles,
        'method': method
    }


def visualize_particle_histograms(analysis_results, data_dir, dataset_name, max_plots=9):
    """
    Create histogram visualizations for the first few particles showing
    perimeter intensity distributions and detected peaks.
    """
    particles = analysis_results['particles']
    n_to_plot = min(max_plots, len(particles))

    if n_to_plot == 0:
        return

    # Create subplots
    n_cols = 3
    n_rows = (n_to_plot + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5 * n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    for idx in range(n_to_plot):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        particle = particles[idx]
        peak_info = particle['cap_peak_info']

        if peak_info is None:
            ax.text(0.5, 0.5, f"Particle {particle['particle_id']}\nNo data",
                   ha='center', va='center')
            ax.axis('off')
            continue

        # Plot histogram
        ax.bar(peak_info['bin_centers'], peak_info['hist'],
               width=np.diff(peak_info['bin_centers'])[0] * 0.8,
               alpha=0.6, label='Histogram')

        # Plot smoothed histogram
        if 'hist_smooth' in peak_info:
            ax.plot(peak_info['bin_centers'], peak_info['hist_smooth'],
                   'r-', linewidth=2, label='Smoothed')

        # Mark detected peaks
        if len(peak_info['peaks']) > 0:
            for peak_val in peak_info['peaks']:
                ax.axvline(peak_val, color='green', linestyle='--', linewidth=2)

        # Mark background value
        ax.axvline(peak_info['background_value'], color='blue',
                  linestyle='-', linewidth=2, label=f"BG={peak_info['background_value']:.1f}")

        ax.set_title(f"Particle {particle['particle_id']}\n"
                    f"{peak_info['n_peaks']} peak(s), {particle['n_perimeter_pixels']} px")
        ax.set_xlabel('Cap Intensity')
        ax.set_ylabel('Count')
        ax.legend(fontsize=8)

    # Hide unused subplots
    for idx in range(n_to_plot, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout()
    output_path = Path(data_dir) / f"{dataset_name}_perimeter_histograms.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path.name}")


def save_summary_statistics(analysis_results, data_dir, dataset_name):
    """Save summary statistics to a CSV file."""
    import csv

    particles = analysis_results['particles']
    output_path = Path(data_dir) / f"{dataset_name}_enrichment_analysis.csv"

    with open(output_path, 'w', newline='') as f:
        if len(particles) == 0:
            return

        # Get fieldnames from first particle (exclude peak_info objects)
        fieldnames = [k for k in particles[0].keys() if 'peak_info' not in k]

        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for particle in particles:
            # Filter out peak_info
            row = {k: v for k, v in particle.items() if 'peak_info' not in k}
            writer.writerow(row)

    print(f"  Saved: {output_path.name}")


def main():
    """Main analysis function."""
    import sys

    # Get base directory from command line argument or use default
    if len(sys.argv) > 1:
        base_dir = sys.argv[1]
    else:
        base_dir = "/Volumes/NX-01-A/2025-10-08_test_data"

    base_path = Path(base_dir)

    if not base_path.exists():
        print(f"Error: Base directory '{base_dir}' does not exist")
        return

    # Find all subdirectories containing .tif files
    subdirs = [d for d in base_path.iterdir() if d.is_dir()]

    if len(subdirs) == 0:
        print(f"Error: No subdirectories found in '{base_dir}'")
        return

    print(f"Found {len(subdirs)} subdirectories to analyze")
    print()

    # Process each subdirectory
    for data_dir in subdirs:
        dataset_name = data_dir.name
        method = "gaussian_peaks"  # Default method

        print("=" * 60)
        print(f"Analyzing {dataset_name} Dataset - Method: {method}")
        print("=" * 60)

        # Load intensity images
        cap_files = [f for f in data_dir.glob("*.tif") if "cap" in f.name.lower() and "intensity" in f.name.lower()]
        g3bp1_files = [f for f in data_dir.glob("*.tif") if "g3bp1" in f.name.lower() and "intensity" in f.name.lower()]

        if len(cap_files) == 0 or len(g3bp1_files) == 0:
            print(f"  Skipping {dataset_name}: Missing intensity images")
            print()
            continue

        cap_intensity = load_image(str(cap_files[0]))
        g3bp1_intensity = load_image(str(g3bp1_files[0]))
        image_shape = cap_intensity.shape

        # Load ROI files
        mask_roi_files = [f for f in data_dir.glob("*.zip") if "mask" in f.name.lower() and "dilated" not in f.name.lower() and "perimeter" not in f.name.lower()]
        perimeter_roi_files = [f for f in data_dir.glob("*.zip") if "perimeter" in f.name.lower()]

        # If no perimeter ROI files, look for dilated mask ROI files
        if len(perimeter_roi_files) == 0:
            perimeter_roi_files = [f for f in data_dir.glob("*.zip") if "dilated" in f.name.lower()]

        if len(mask_roi_files) == 0 or len(perimeter_roi_files) == 0:
            print(f"  Skipping {dataset_name}: Missing ROI files")
            print()
            continue

        print(f"  Loading ROIs from: {mask_roi_files[0].name}")
        print(f"  Loading perimeter ROIs from: {perimeter_roi_files[0].name}")

        mask_rois = load_imagej_rois(str(mask_roi_files[0]))
        perimeter_rois = load_imagej_rois(str(perimeter_roi_files[0]))

        # Analyze enrichment using ROIs
        results = analyze_particle_enrichment_from_rois(
            mask_rois, perimeter_rois, cap_intensity, g3bp1_intensity,
            image_shape, method=method
        )

        # Print summary
        print(f"\n  Summary Statistics (Cap enrichment):")
        cap_bg_subtracted = [p['cap_mean_bg_subtracted'] for p in results['particles']]
        if len(cap_bg_subtracted) > 0:
            print(f"    Mean bg-subtracted: {np.mean(cap_bg_subtracted):.2f}")
            print(f"    Median bg-subtracted: {np.median(cap_bg_subtracted):.2f}")
            print(f"    Particles with enrichment (>0): {np.sum(np.array(cap_bg_subtracted) > 0)}/{len(cap_bg_subtracted)}")

        # Save results
        save_summary_statistics(results, data_dir, dataset_name)
        visualize_particle_histograms(results, data_dir, dataset_name)

        print()


if __name__ == "__main__":
    main()
