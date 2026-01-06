#!/bin/bash
# Base directory
BASE_DIR="/Volumes/NX-01-A/2025-12-17_cap_analysis"

# Clean up macOS dot files before processing
echo "Running dot_clean on ${BASE_DIR}..."
dot_clean "$BASE_DIR"
echo ""

# Find all directories that end with _BS
for bs_dir in "${BASE_DIR}"/*_BS; do
    # Get the base name by removing _BS suffix
    base_name="${bs_dir%_BS}"
    
    # Check if the corresponding base directory exists
    if [ -d "$base_name" ]; then
        echo "Processing: $(basename "$base_name")"
        
        # Find all ch0 .tif files in raw_data
        raw_data_dir="${base_name}/raw_data"
        processed_dir="${bs_dir}/Processed"
        
        if [ -d "$raw_data_dir" ] && [ -d "$processed_dir" ]; then
            # Find all ch0 files
            while IFS= read -r -d '' ch0_file; do
                # Extract the filename
                filename=$(basename "$ch0_file")
                
                # Determine the condition (1Hr_NaAsO2 or Untreated)
                if [[ "$filename" == *"1Hr_NaAsO2"* ]]; then
                    condition="1Hr_NaAsO2"
                elif [[ "$filename" == *"Untreated"* ]]; then
                    condition="Untreated"
                else
                    echo "  Warning: Could not determine condition for $filename"
                    continue
                fi
                
                # Check if the subdirectory exists in Processed
                dest_subdir="${processed_dir}/${condition}"
                if [ -d "$dest_subdir" ]; then
                    # Create new filename: Condition_Cap_Intensity.tif
                    new_filename="${condition}_Cap_Intensity.tif"
                    dest_file="${dest_subdir}/${new_filename}"
                    
                    echo "  Copying: $filename -> ${condition}/${new_filename}"
                    cp -p "$ch0_file" "$dest_file"
                    
                    # Verify the copy
                    if [ -f "$dest_file" ]; then
                        src_size=$(stat -f%z "$ch0_file" 2>/dev/null || stat -c%s "$ch0_file" 2>/dev/null)
                        dest_size=$(stat -f%z "$dest_file" 2>/dev/null || stat -c%s "$dest_file" 2>/dev/null)
                        if [ "$src_size" != "$dest_size" ]; then
                            echo "  ERROR: File size mismatch for $new_filename"
                        fi
                    fi
                else
                    echo "  Warning: Destination directory $dest_subdir does not exist"
                fi
            done < <(find "$raw_data_dir" -name "*ch0*.tif" -print0)
        else
            echo "  Warning: Required directories not found"
            [ ! -d "$raw_data_dir" ] && echo "    Missing: $raw_data_dir"
            [ ! -d "$processed_dir" ] && echo "    Missing: $processed_dir"
        fi
        
        echo ""
    else
        echo "Warning: Base directory $base_name not found for $bs_dir"
    fi
done
echo "Copy and rename operation completed!"