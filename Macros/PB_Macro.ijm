// ImageJ Macro for P-body Analysis with Stress Granule Subtraction
// This macro processes microscopy data to subtract stress granule regions from p-bodies

// Set batch mode for faster processing
setBatchMode(true);

// Get the main directory containing Masks and Raw Data folders
mainDir = getDirectory("Choose the main directory (containing Masks and Raw Data folders)");
masksDir = mainDir + "Masks" + File.separator;
rawDataDir = mainDir + "Raw Data" + File.separator;
outputBaseDir = mainDir + "Processed" + File.separator;

// Create output directory if it doesn't exist
File.makeDirectory(outputBaseDir);

// Get list of all mask files
maskFiles = getFileList(masksDir);

// Create an array to store unique condition strings
conditionList = newArray();

// Extract unique conditions from ch2 files (p-bodies)
for (i = 0; i < maskFiles.length; i++) {
    if (indexOf(maskFiles[i], "_ch2_") > 0 && endsWith(maskFiles[i], ".tif")) {
        // Extract the unique string part
        filename = maskFiles[i];
        // Remove "MASK_MAX_z-stack_" prefix
        temp = replace(filename, "MASK_MAX_z-stack_", "");
        // Remove "_Merged_ch2_t00.tif" suffix
        condition = replace(temp, "_Merged_ch2_t00.tif", "");
        conditionList = Array.concat(conditionList, condition);
    }
}

// Process each condition
for (c = 0; c < conditionList.length; c++) {
    condition = conditionList[c];
    print("Processing condition: " + condition);
    
    // Define file paths
    ch1MaskPath = masksDir + "MASK_MAX_z-stack_" + condition + "_Merged_ch1_t00.tif";
    ch2MaskPath = masksDir + "MASK_MAX_z-stack_" + condition + "_Merged_ch2_t00.tif";
    ch2RawPath = rawDataDir + "MAX_z-stack_" + condition + "_Merged_ch2_t00.tif";
    
    // Check if all required files exist
    if (!File.exists(ch1MaskPath) || !File.exists(ch2MaskPath) || !File.exists(ch2RawPath)) {
        print("Skipping " + condition + " - missing files");
        continue;
    }
    
    // Create output directory for this condition
    conditionOutputDir = outputBaseDir + condition + File.separator;
    File.makeDirectory(conditionOutputDir);
    
    // Clear ROI Manager
    roiManager("reset");
    
    // Step 1: Open ch1 mask (stress granules) and get ROIs
    open(ch1MaskPath);
    ch1Title = getTitle();
    run("Analyze Particles...", "add");
    numSGRois = roiManager("count");
    print("  Found " + numSGRois + " stress granule ROIs");
    
    // Close ch1
    close(ch1Title);
    
	// Step 2: Open ch2 mask (p-bodies)
	open(ch2MaskPath);
	ch2Title = getTitle();
	getStatistics(area, mean, min, max);
	print("  Ch2 image type: " + bitDepth() + "-bit");
	print("  Ch2 min/max: " + min + "/" + max);
	getThreshold(lower, upper);
	print("  Ch2 initial threshold: " + lower + "-" + upper);
	
	// Set threshold for binary mask (in case it's not already set)
	setThreshold(1, 255);
	setOption("BlackBackground", true);
	
	// Step 3: Subtract stress granule regions from p-bodies
	if (numSGRois > 0) {
	    roiManager("Combine");
	    run("Clear", "slice");
	    run("Select None");
	}
	
	// Clear ROI manager and analyze the subtracted p-body mask
	roiManager("reset");
	run("Analyze Particles...", "add");
    numPBodyRois = roiManager("count");
    print("  Found " + numPBodyRois + " p-body ROIs after subtraction");
    
    // Step 4: Save the p-body ROI list (with "Mask" in name)
    if (numPBodyRois > 0) {
        roiOutputPath = conditionOutputDir + condition + "_PB_Mask.zip";
        roiManager("Save", roiOutputPath);
        print("  Saved ROI list: " + roiOutputPath);
        
        // Step 5: Enlarge all ROIs by 2 pixels
        for (r = 0; r < numPBodyRois; r++) {
            roiManager("select", r);
            run("Enlarge...", "enlarge=2 pixel");
            roiManager("update");
        }
        
        // Step 6: Save the dilated ROI list (with "Dilated" and "Mask" in name)
        dilatedRoiOutputPath = conditionOutputDir + condition + "_PB_Dilated_Mask.zip";
        roiManager("Save", dilatedRoiOutputPath);
        print("  Saved dilated ROI list: " + dilatedRoiOutputPath);
    }
    
    // Close ch2 mask
    close(ch2Title);
    
    // Step 7: Copy the ch2 raw data image to output directory (with "Intensity" in name)
    rawOutputPath = conditionOutputDir + condition + "_DDX6_Intensity.tif";
    File.copy(ch2RawPath, rawOutputPath);
    print("  Copied raw data: " + rawOutputPath);
    
    // Clear ROI manager for next iteration
    roiManager("reset");
    print("Completed: " + condition);
    print("---");
}

// Exit batch mode
setBatchMode(false);

print("=================================");
print("Processing complete!");
print("Output saved to: " + outputBaseDir);
print("=================================");