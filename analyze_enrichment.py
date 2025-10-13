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
                    'coordinates': roi.coordinates(),
                    'bytes': roi_bytes  # Store original bytes for ImageJ processing
                })
    return rois


def enlarge_rois_with_imagej(rois, pixels):
    """
    Enlarge ROIs using binary dilation (equivalent to ImageJ's RoiEnlarger).

    Args:
        rois: List of ROI dictionaries from load_imagej_rois
        pixels: Number of pixels to enlarge (positive) or shrink (negative)

    Returns:
        List of enlarged ROI dictionaries with updated coordinates

    Note: This uses scipy's binary dilation which produces similar results to
    ImageJ's RoiEnlarger but is faster and doesn't require Java/ImageJ initialization.
    """
    print(f"  Enlarging ROIs by {pixels} pixels using binary dilation...")

    # We need image_shape to create masks. Get it from the first successful coordinate conversion
    # For now, we'll determine it from the max coordinates
    max_x = max_y = 0
    for roi_dict in rois:
        coords = roi_dict['coordinates']
        if len(coords) > 0:
            max_x = max(max_x, int(np.max(coords[:, 0])))
            max_y = max(max_y, int(np.max(coords[:, 1])))

    # Add some padding to ensure ROIs fit
    image_shape = (max_y + 100, max_x + 100)

    enlarged_rois = []

    for roi_dict in rois:
        # Convert ROI to mask
        mask = roi_to_mask(roi_dict['coordinates'], image_shape)

        # Create structuring element for dilation
        structure = ndimage.generate_binary_structure(2, 1)  # 4-connected

        # Enlarge (dilate) or shrink (erode) the mask
        if pixels > 0:
            enlarged_mask = ndimage.binary_dilation(mask, structure=structure,
                                                   iterations=pixels).astype(np.uint8)
        elif pixels < 0:
            enlarged_mask = ndimage.binary_erosion(mask, structure=structure,
                                                  iterations=abs(pixels)).astype(np.uint8)
        else:
            enlarged_mask = mask

        # Extract new coordinates from the enlarged mask boundary
        from skimage.measure import find_contours
        contours = find_contours(enlarged_mask, 0.5)

        if len(contours) > 0:
            # Use the largest contour
            largest_contour = max(contours, key=len)
            # Convert from (row, col) to (x, y)
            new_coords = np.column_stack((largest_contour[:, 1], largest_contour[:, 0]))
        else:
            # If no contour found, use original coordinates
            new_coords = roi_dict['coordinates']

        enlarged_rois.append({
            'name': roi_dict['name'],
            'roi': roi_dict['roi'],  # Keep original roi object
            'coordinates': new_coords,
            'bytes': roi_dict.get('bytes', b'')
        })

    print(f"  Successfully enlarged {len(enlarged_rois)} ROIs by {pixels} pixels")
    return enlarged_rois


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


def find_gaussian_peaks(data, n_bins=50, max_background=None):
    """
    Find Gaussian peaks in a histogram of the data.

    Args:
        data: Array of intensity values
        n_bins: Number of bins for histogram
        max_background: If specified, only consider peaks below this value for background

    Returns:
        Dictionary with peak information and histogram data
    """
    if len(data) == 0:
        return None

    # Create histogram with range starting at 0 to capture near-zero peaks
    data_max = float(np.max(data))
    hist, bin_edges = np.histogram(data, bins=n_bins, range=(0, data_max))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Smooth the histogram slightly to reduce noise
    from scipy.ndimage import gaussian_filter1d
    hist_smooth = gaussian_filter1d(hist.astype(float), sigma=2)

    # Find peaks in the histogram with lower prominence to catch near-zero peaks
    peaks, properties = find_peaks(hist_smooth, prominence=np.max(hist_smooth) * 0.05)

    if len(peaks) == 0:
        # Try even lower prominence threshold
        peaks, properties = find_peaks(hist_smooth, prominence=1)

    if len(peaks) == 0:
        # Still no peaks found, use the bin with maximum count as the peak
        peak_idx = np.argmax(hist_smooth)
        background_value = bin_centers[peak_idx]
        return {
            'n_peaks': 1,
            'peaks': [background_value],
            'hist': hist,
            'bin_centers': bin_centers,
            'hist_smooth': hist_smooth,
            'background_value': background_value
        }

    # Get the peak positions (intensity values)
    peak_positions = bin_centers[peaks]
    peak_prominences = properties['prominences']

    # Determine background value
    if max_background is not None:
        # Constrained mode: Only consider peaks below max_background threshold
        peaks_below_threshold = peak_positions < max_background

        if np.any(peaks_below_threshold):
            # Select the most prominent peak below threshold
            prominences_below = peak_prominences[peaks_below_threshold]
            positions_below = peak_positions[peaks_below_threshold]
            most_prominent_idx = np.argmax(prominences_below)
            background_value = positions_below[most_prominent_idx]
        else:
            # No peaks below threshold, find bin with highest count below threshold
            bins_below_threshold = bin_centers < max_background
            if np.any(bins_below_threshold):
                # Get histogram counts for bins below threshold
                hist_below = hist[bins_below_threshold]
                bin_centers_below = bin_centers[bins_below_threshold]
                max_count_idx = np.argmax(hist_below)
                background_value = bin_centers_below[max_count_idx]
            else:
                # Fallback: use the leftmost peak
                sorted_by_position = np.argsort(peak_positions)
                background_value = peak_positions[sorted_by_position[0]]
    else:
        # Default mode: Use the leftmost (lowest intensity) peak
        sorted_by_position = np.argsort(peak_positions)
        background_value = peak_positions[sorted_by_position[0]]

    # For reporting, also sort by prominence to return the most prominent peaks
    peak_prominences = properties['prominences']
    sorted_indices = np.argsort(peak_prominences)[::-1]
    peaks_sorted_by_prominence = peaks[sorted_indices]
    peak_positions_for_report = bin_centers[peaks_sorted_by_prominence]

    return {
        'n_peaks': len(peaks),
        'peaks': peak_positions_for_report[:2] if len(peaks) >= 2 else peak_positions_for_report,
        'hist': hist,
        'bin_centers': bin_centers,
        'hist_smooth': hist_smooth,
        'peak_indices': peaks_sorted_by_prominence[:2] if len(peaks) >= 2 else peaks_sorted_by_prominence,
        'background_value': background_value
    }


def analyze_particle_enrichment_from_rois(mask_rois, perimeter_rois, cap_intensity,
                                          image_shape, method='gaussian_peaks',
                                          bg_multiplication_factor=1.0, max_background=None):
    """
    Analyze Cap enrichment for each particle using ROI definitions.

    Args:
        mask_rois: List of ROI dictionaries for particles
        perimeter_rois: List of ROI dictionaries for perimeters
        cap_intensity: Cap intensity image
        image_shape: Shape of the images
        method: 'minimum' or 'gaussian_peaks'
        bg_multiplication_factor: Multiplication factor for background value (default 1.0)

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

        # Skip if no pixels found
        if len(cap_particle_intensities) == 0 or len(cap_perimeter_intensities) == 0:
            continue

        # Determine background value based on method
        if method == 'minimum':
            cap_background_raw = np.min(cap_perimeter_intensities)
            cap_peak_info = None
        elif method == 'gaussian_peaks':
            cap_peak_info = find_gaussian_peaks(cap_perimeter_intensities, max_background=max_background)
            cap_background_raw = cap_peak_info['background_value'] if cap_peak_info else np.mean(cap_perimeter_intensities)
        else:
            raise ValueError(f"Unknown method: {method}")

        # Apply multiplication factor to background value
        cap_background = cap_background_raw * bg_multiplication_factor

        # Background-subtracted intensities
        cap_bg_subtracted = cap_particle_intensities - cap_background

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
            # Filter out peak_info and convert negative values to empty strings
            row = {}
            for k, v in particle.items():
                if 'peak_info' not in k:
                    # If the value is numeric and negative, use empty string
                    if isinstance(v, (int, float, np.integer, np.floating)) and v < 0:
                        row[k] = ''
                    else:
                        row[k] = v
            writer.writerow(row)

    print(f"  Saved: {output_path.name}")


def main():
    """Main analysis function."""
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(description='P-body Cap Enrichment Analysis')
    parser.add_argument('base_dir', nargs='?', default="/Volumes/NX-01-A/2025-10-08_test_data",
                       help='Base directory containing subdirectories with microscopy data')
    parser.add_argument('--bg-factor', type=float, default=1.0,
                       help='Background subtraction multiplication factor (default: 1.0)')
    parser.add_argument('--enlarge-rois', type=int, default=0,
                       help='Number of pixels to enlarge dilated ROIs (default: 0, no enlargement)')
    parser.add_argument('--max-background', type=float, default=None,
                       help='Maximum background value threshold. Only peaks below this value will be considered for background (default: None, no constraint)')

    args = parser.parse_args()

    base_dir = args.base_dir
    bg_multiplication_factor = args.bg_factor
    roi_enlargement_pixels = args.enlarge_rois
    max_background = args.max_background

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
    print(f"Background multiplication factor: {bg_multiplication_factor}")
    if roi_enlargement_pixels > 0:
        print(f"ROI enlargement: {roi_enlargement_pixels} pixels")
    else:
        print("ROI enlargement: disabled")
    if max_background is not None:
        print(f"Maximum background constraint: {max_background}")
    else:
        print("Maximum background constraint: disabled")
    print()

    # Process each subdirectory
    for data_dir in subdirs:
        dataset_name = data_dir.name
        method = "gaussian_peaks"  # Default method

        print("=" * 60)
        print(f"Analyzing {dataset_name} Dataset - Method: {method}")
        print(f"Background multiplication factor: {bg_multiplication_factor}")
        if roi_enlargement_pixels > 0:
            print(f"ROI enlargement: {roi_enlargement_pixels} pixels")
        if max_background is not None:
            print(f"Maximum background constraint: {max_background}")
        print("=" * 60)

        # Load intensity images
        cap_files = [f for f in data_dir.glob("*.tif") if "cap" in f.name.lower() and "intensity" in f.name.lower()]

        if len(cap_files) == 0:
            print(f"  Skipping {dataset_name}: Missing Cap intensity image")
            print()
            continue

        cap_intensity = load_image(str(cap_files[0]))
        image_shape = cap_intensity.shape

        # Load ROI files
        mask_roi_files = [f for f in data_dir.glob("*.zip") if "mask" in f.name.lower() and "dilated" not in f.name.lower()]
        dilated_roi_files = [f for f in data_dir.glob("*.zip") if "dilated" in f.name.lower() and "mask" in f.name.lower()]

        if len(mask_roi_files) == 0 or len(dilated_roi_files) == 0:
            print(f"  Skipping {dataset_name}: Missing ROI files")
            print()
            continue

        print(f"  Loading ROIs from: {mask_roi_files[0].name}")
        print(f"  Loading dilated ROIs from: {dilated_roi_files[0].name}")

        mask_rois = load_imagej_rois(str(mask_roi_files[0]))
        dilated_rois = load_imagej_rois(str(dilated_roi_files[0]))

        # Optionally further enlarge the dilated ROIs using ImageJ's RoiEnlarger
        if roi_enlargement_pixels > 0:
            perimeter_rois = enlarge_rois_with_imagej(dilated_rois, pixels=roi_enlargement_pixels)
        else:
            perimeter_rois = dilated_rois

        # Analyze enrichment using ROIs
        results = analyze_particle_enrichment_from_rois(
            mask_rois, perimeter_rois, cap_intensity,
            image_shape, method=method, bg_multiplication_factor=bg_multiplication_factor,
            max_background=max_background
        )

        # Print summary
        print(f"\n  Summary Statistics (Cap enrichment):")
        cap_backgrounds = [p['cap_background'] for p in results['particles']]
        cap_bg_subtracted = [p['cap_mean_bg_subtracted'] for p in results['particles']]
        if len(cap_backgrounds) > 0:
            print(f"    Mean background: {np.mean(cap_backgrounds):.2f}")
            print(f"    Median background: {np.median(cap_backgrounds):.2f}")
            print(f"    Particles with enrichment (>0): {np.sum(np.array(cap_bg_subtracted) > 0)}/{len(cap_bg_subtracted)}")

        # Save results
        save_summary_statistics(results, data_dir, dataset_name)
        visualize_particle_histograms(results, data_dir, dataset_name)

        print()


if __name__ == "__main__":
    main()
