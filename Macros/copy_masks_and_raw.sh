#!/bin/bash
# Base directory
BASE_DIR="/Volumes/NX-01-A/2025-12-16_cap_analysis"

# Find all directories that DON'T end with _BS (these are the source directories)
for base_name in "${BASE_DIR}"/*/; do
    # Remove trailing slash
    base_name="${base_name%/}"
    
    # Skip if this is a _BS directory
    [[ "$base_name" == *_BS ]] && continue
    
    echo "Processing: $(basename "$base_name")"
    
    # Define the _BS directory path
    bs_dir="${base_name}_BS"
    
    # Create _BS directory and subdirectories if they don't exist
    mkdir -p "$bs_dir"
    mkdir -p "${bs_dir}/Masks"
    mkdir -p "${bs_dir}/Raw Data"
    echo "  Ensured directories exist: $(basename "$bs_dir"), Masks, Raw Data"
    
    # Copy combined_masks .tif files to Masks
    masks_source="${base_name}/combined_masks"
    masks_dest="${bs_dir}/Masks"
    
    if [ -d "$masks_source" ]; then
        echo "  Copying masks from $masks_source to $masks_dest"
        find "$masks_source" -name "*.tif" -exec cp -v {} "$masks_dest/" \;
    else
        echo "  Warning: $masks_source not found"
    fi
    
    # Copy ch1 and ch2 .tif files from raw_data to Raw Data
    raw_source="${base_name}/raw_data"
    raw_dest="${bs_dir}/Raw Data"
    
    if [ -d "$raw_source" ]; then
        echo "  Copying ch1 and ch2 files from $raw_source to $raw_dest"
        find "$raw_source" -name "*ch1*.tif" -exec cp -v {} "$raw_dest/" \;
        find "$raw_source" -name "*ch2*.tif" -exec cp -v {} "$raw_dest/" \;
    else
        echo "  Warning: $raw_source not found"
    fi
    
    echo ""
done

# Clean up macOS dot files
echo "Running dot_clean on ${BASE_DIR}..."
dot_clean "$BASE_DIR"

echo "Copy operation completed!"