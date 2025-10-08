#!/usr/bin/env python3
"""
Microscopy Data Extraction Script
Extracts Cap and G3BP1 intensity values from mask regions and perimeters.
"""

import numpy as np
from pathlib import Path
from PIL import Image
import tifffile
import matplotlib.pyplot as plt
import glob


def load_image(filepath):
    """Load an image file and convert to numpy array."""
    img = Image.open(filepath)
    return np.array(img)


def find_file(directory, keywords):
    """
    Find a file in directory that contains all specified keywords.

    Args:
        directory: Path to search in
        keywords: List of keywords that must all be in filename

    Returns:
        Path to matching file
    """
    files = list(Path(directory).glob("*.tif"))
    for file in files:
        if all(keyword.lower() in file.name.lower() for keyword in keywords):
            return str(file)
    raise FileNotFoundError(f"No file found with keywords: {keywords}")


def extract_intensity_data(data_dir, save_perimeter_mask=True, visualize=True):
    """
    Extract microscopy intensity data from a dataset directory.

    Args:
        data_dir: Path to directory containing mask and intensity images
        save_perimeter_mask: If True, saves the perimeter mask to the data directory

    Returns:
        Dictionary containing:
            - cap_mask_intensities: Cap intensity values within mask
            - g3bp1_mask_intensities: G3BP1 intensity values within mask
            - cap_perimeter_intensities: Cap intensity values at perimeter
            - g3bp1_perimeter_intensities: G3BP1 intensity values at perimeter
            - perimeter_mask: Binary mask of perimeter regions
    """
    # Find files
    mask_file = find_file(data_dir, ["Mask"])
    # Make sure we get the non-dilated mask (exclude "Dilated")
    files = list(Path(data_dir).glob("*.tif"))
    mask_candidates = [f for f in files if "mask" in f.name.lower() and "dilated" not in f.name.lower()]
    if mask_candidates:
        mask_file = str(mask_candidates[0])

    dilated_mask_file = find_file(data_dir, ["Dilated", "Mask"])
    cap_file = find_file(data_dir, ["Cap", "Intensity"])
    g3bp1_file = find_file(data_dir, ["G3BP1", "Intensity"])

    print(f"Loading files from {data_dir}:")
    print(f"  Mask: {Path(mask_file).name}")
    print(f"  Dilated Mask: {Path(dilated_mask_file).name}")
    print(f"  Cap Intensity: {Path(cap_file).name}")
    print(f"  G3BP1 Intensity: {Path(g3bp1_file).name}")

    # Load images
    mask = load_image(mask_file)
    dilated_mask = load_image(dilated_mask_file)
    cap_intensity = load_image(cap_file)
    g3bp1_intensity = load_image(g3bp1_file)

    # Convert masks to binary (0 or 1)
    mask_binary = (mask > 0).astype(np.uint8)
    dilated_mask_binary = (dilated_mask > 0).astype(np.uint8)

    # Create perimeter mask by subtracting dilated from mask (or mask from dilated)
    # The perimeter is the absolute difference, giving the ring/donut shape
    perimeter_mask = np.abs(mask_binary.astype(np.int16) - dilated_mask_binary.astype(np.int16)).astype(np.uint8)

    # Save perimeter mask if requested
    if save_perimeter_mask:
        # Convert to 8-bit (0 or 255 for visibility)
        perimeter_mask_8bit = (perimeter_mask * 255).astype(np.uint8)

        # Create output filename based on the original mask filename
        mask_basename = Path(mask_file).stem
        output_filename = f"{mask_basename.replace('Mask', 'Perimeter Mask')}.tif"
        output_path = Path(data_dir) / output_filename

        # Save the perimeter mask using tifffile to preserve exact values
        tifffile.imwrite(output_path, perimeter_mask_8bit)
        print(f"  Saved perimeter mask: {output_filename}")

    # Extract intensity values
    # For mask regions: get all non-zero values
    cap_mask_intensities = cap_intensity[mask_binary > 0]
    g3bp1_mask_intensities = g3bp1_intensity[mask_binary > 0]

    # For perimeter regions: get all values where perimeter mask is non-zero
    cap_perimeter_intensities = cap_intensity[perimeter_mask > 0]
    g3bp1_perimeter_intensities = g3bp1_intensity[perimeter_mask > 0]

    print(f"\nExtracted data summary:")
    print(f"  Cap mask pixels: {len(cap_mask_intensities)}")
    print(f"  G3BP1 mask pixels: {len(g3bp1_mask_intensities)}")
    print(f"  Cap perimeter pixels: {len(cap_perimeter_intensities)}")
    print(f"  G3BP1 perimeter pixels: {len(g3bp1_perimeter_intensities)}")

    # Visualize perimeter intensities if requested
    if visualize and len(cap_perimeter_intensities) > 0:
        dataset_name = Path(data_dir).name

        # Create masked arrays where only perimeter pixels are shown
        cap_perimeter_img = np.zeros_like(cap_intensity, dtype=float)
        cap_perimeter_img[perimeter_mask > 0] = cap_intensity[perimeter_mask > 0]

        g3bp1_perimeter_img = np.zeros_like(g3bp1_intensity, dtype=float)
        g3bp1_perimeter_img[perimeter_mask > 0] = g3bp1_intensity[perimeter_mask > 0]

        # Set zeros to NaN for better visualization
        cap_perimeter_img[cap_perimeter_img == 0] = np.nan
        g3bp1_perimeter_img[g3bp1_perimeter_img == 0] = np.nan

        # Plot and save Cap intensity - full image
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(cap_intensity, cmap='viridis')
        ax.set_title('Cap Intensity - Full Image')
        ax.axis('off')
        plt.colorbar(im, ax=ax)
        output_path = Path(data_dir) / f"{dataset_name}_Cap_Full.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path.name}")

        # Plot and save Cap perimeter intensities
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(cap_perimeter_img, cmap='viridis')
        ax.set_title('Cap Intensity - Perimeter Only')
        ax.axis('off')
        plt.colorbar(im, ax=ax)
        output_path = Path(data_dir) / f"{dataset_name}_Cap_Perimeter.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path.name}")

        # Plot and save G3BP1 intensity - full image
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(g3bp1_intensity, cmap='viridis')
        ax.set_title('G3BP1 Intensity - Full Image')
        ax.axis('off')
        plt.colorbar(im, ax=ax)
        output_path = Path(data_dir) / f"{dataset_name}_G3BP1_Full.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path.name}")

        # Plot and save G3BP1 perimeter intensities
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(g3bp1_perimeter_img, cmap='viridis')
        ax.set_title('G3BP1 Intensity - Perimeter Only')
        ax.axis('off')
        plt.colorbar(im, ax=ax)
        output_path = Path(data_dir) / f"{dataset_name}_G3BP1_Perimeter.png"
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_path.name}")

    return {
        'cap_mask_intensities': cap_mask_intensities,
        'g3bp1_mask_intensities': g3bp1_mask_intensities,
        'cap_perimeter_intensities': cap_perimeter_intensities,
        'g3bp1_perimeter_intensities': g3bp1_perimeter_intensities,
        'perimeter_mask': perimeter_mask,
        'mask': mask_binary,
        'dilated_mask': dilated_mask_binary
    }


def main():
    """Main function to process both datasets."""
    base_dir = "/Volumes/NX-01-A/2025-10-08_test_data"

    # Process Untreated dataset
    print("=" * 60)
    print("Processing Untreated Dataset")
    print("=" * 60)
    untreated_dir = Path(base_dir) / "Untreated"
    untreated_data = extract_intensity_data(untreated_dir)

    # Process Arsenite Treated dataset
    print("\n" + "=" * 60)
    print("Processing Arsenite Treated Dataset")
    print("=" * 60)
    as_treated_dir = Path(base_dir) / "As Treated"
    as_treated_data = extract_intensity_data(as_treated_dir)

    print("\n" + "=" * 60)
    print("Extraction Complete!")
    print("=" * 60)

    return {
        'untreated': untreated_data,
        'as_treated': as_treated_data
    }


if __name__ == "__main__":
    data = main()
