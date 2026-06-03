//@ File    (label="Input directory", style="directory", description="Folder containing the raw images to process.") inputDir

// @ String (label="File suffix", choices={".nd2", ".czi", ".tif"}, style="listBox") fileExtension

//@ Integer (label="Crop size in pixels", value=1024, min=1, stepSize=1, description="Square crop size in pixels, e.g. 1024 creates 1024 x 1024 crops.") cropSize
//@ Integer (label="Reference channel", value=3, min=1, stepSize=1, description="Channel used for ROI definition, filtering, and exported crop. Use 0 to use all channels") refChannel
//@ Integer (label="Number of crops per image", value=1, min=1, stepSize=1, description="Number of crops generated from each input image.") nCrops

//@ String  (label="Crop strategy", choices={"random", "random+filtered", "intensity_based"}, value="random+filtered", description="Strategy used to select crop positions.") strategy

//@ Boolean (label="Use MIP for filtering", value=true, description="If enabled, filtering intensity is measured on a maximum intensity projection of the selected channel. If disabled, filtering is measured on the current slice.") useMIPforFiltering

//@ Boolean (label="Reject near-empty crops", value=true, description="If enabled, crops with too little signal variation are rejected using the minimum standard deviation threshold.") rejectNearEmpty
//@ Double  (label="Minimum standard deviation", value=5, min=0, stepSize=0.1, description="Minimum standard deviation required to keep a crop when near-empty rejection is enabled.") minStdDev

//@ Double  (label="Minimum mean intensity", value=100, min=0, stepSize=1, description="Used by the random+filtered strategy. Crops below this mean intensity are rejected.") minMeanIntensity

//@ Double  (label="Top fraction", value=0.30, min=0, max=1, stepSize=0.01, style="format:0.00", description="Used by the intensity_based strategy. Fraction of brightest candidate regions considered.") topFraction

//@ Boolean (label="Normalize before cropping", value=false, description="If enabled, the selected channel stack is normalized before crop extraction.") normalizeBeforeCropping
//@ Double  (label="Saturation", value=0.01, min=0, max=100, stepSize=0.01, style="format:0.00", description="Saturation value used during normalization.") saturation

// @ Float(label="XY scaling factor", value=1.0) scaleXY
// @ Float(label="Z scaling factor", value=1.0) scaleZ
// @ Boolean(label="Use interpolation", value=false) interpolation

//@ String  (label="Output format", choices={"TIFF", "HDF5"}, value="TIFF", description="Format used to save the extracted crops.") outputFormat

//@ Boolean (label="Use HDF5-ready naming", value=true, description="If enabled, crop names are padded: basename_001, basename_002, etc.") hdf5ReadyNaming

//@ Integer (label="HDF5 compression level", value=0, min=0, max=9, stepSize=1, description="Compression level used for HDF5 output. 0 means no compression; 9 is maximum compression.") compressionLevel


/* 25/05/2026
Nicolas Peredo, PhD
Image Analysis Expert
VIB BioImaging Core Leuven, VIB Technologies, Center for Neuroscience Leuven, Leuven, Belgium
VIB BioImaging Core Leuven, Department of Neurosciences, KU Leuven, Leuven, Belgium
Nikon Center of Excellence
Campus Gasthuisberg - ON5 - room 04.367
Herestraat 49 - box 62
3000 Leuven
Belgium
*/


// ------------------------------------------------------------
// Batch random crops for ilastik training
// Single-channel crop export version
//
// Features:
// - loops through folder
// - stored image IDs
// - strategies: random / random+filtered / intensity_based
// - optional MIP-based filtering
// - optional rejection of near-empty crops
// - exports ONLY the selected channel
// - optional normalization of the selected channel before cropping
// - output as TIFF or HDF5
// ------------------------------------------------------------


// ------------------------------
// GET FILE LIST
// ------------------------------
outputDir = inputDir + "/training_data";
File.makeDirectory(outputDir);

fileList = getFilesList(inputDir, fileExtension);

if (fileList.length == 0)
    exit("No files found with extension: " + fileExtension);

setBatchMode(true);



// ------------------------------
// MAIN LOOP
// ------------------------------

for (f=0; f<fileList.length; f++) {

    filename = fileList[f];
    fullPath = inputDir + File.separator + filename;

    print("\\Clear");
    print("Processing: " + filename);

    // Open image
    run("Bio-Formats Importer", "open=[" + fullPath + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");

    sourceTitle = getTitle();
    sourceID = getImageID();
    basename = getBasename(filename, ".");

    getDimensions(w, h, c, z, t);

    if (refChannel < 1 || refChannel > c) {
        selectImage(sourceID);
        close();
        exit("Selected channel exceeds number of channels.");
    }

    // ----------------------------------------
    // CREATE SINGLE-CHANNEL WORKING STACK
    // This stack is the selected channel only.
    // It is the one optionally normalized and used for cropping/export.
    // ----------------------------------------
    selectImage(sourceID);
    Stack.setChannel(refChannel);

    workTitle = basename + "_ch" + refChannel + "_work";
    run("Duplicate...", "title=[" + workTitle + "] duplicate channels=" + refChannel + "-" + refChannel + " slices=1-" + z + " frames=1-" + t);
    workID = getImageID();
    workTitle = getTitle();

    if (normalizeBeforeCropping) {
        selectImage(workID);
        run("Enhance Contrast...", "saturated=" + saturation + " normalize process_all use");
    }
    
    if (scaleXY != 1.0 || scaleZ != 1.0) {
		originalDimensions = scalingXYZ_multichannel(image_name, scaleXY, scaleZ, interpolation);
		close(workTitle);
		workID = getImageID();
    	workTitle = getTitle();
	}

    // ----------------------------------------
    // PREPARE FILTERING IMAGE
    // Either the working stack itself or its MIP
    // ----------------------------------------
    filterImageID = workID;
    filterImageTitle = workTitle;

    if (useMIPforFiltering) {
        selectImage(workID);

        mipTitle = basename + "_filterMIP";
        run("Duplicate...", "title=[" + mipTitle + "] duplicate");

        mipSourceID = getImageID();

        getDimensions(w2, h2, c2, z2, t2);
        if (z2 > 1) {
            run("Z Project...", "projection=[Max Intensity]");
            filterImageID = getImageID();
            filterImageTitle = getTitle();

            selectImage(mipSourceID);
            close();
        } else {
            filterImageID = mipSourceID;
            filterImageTitle = getTitle();
        }
    }

    // ----------------------------------------
    // CREATE MAXIMUM GRID ROIs
    // ----------------------------------------
    if (!isOpen("ROI Manager"))
        run("ROI Manager...");
    roiManager("Reset");

    cols = floor(w / cropSize);
    rows = floor(h / cropSize);
    maxROIs = cols * rows;

    if (maxROIs < 1) {
        print("Skipping " + filename + ": crop size too large.");

        if (useMIPforFiltering && filterImageID != workID) {
            selectImage(filterImageID);
            close();
        }

        selectImage(workID);
        close();

        selectImage(sourceID);
        close();
        continue;
    }

    candidateX = newArray(maxROIs);
    candidateY = newArray(maxROIs);

    idx = 0;
    for (yy=0; yy<rows; yy++) {
        for (xx=0; xx<cols; xx++) {
            x = xx * cropSize;
            y = yy * cropSize;

            candidateX[idx] = x;
            candidateY[idx] = y;

            makeRectangle(x, y, cropSize, cropSize);
            roiManager("Add");

            idx++;
        }
    }

    // ----------------------------------------
    // SELECT ROI INDICES ACCORDING TO STRATEGY
    // ----------------------------------------
    if (strategy == "random") {
        selectedIndices = strategy_random_with_optional_empty_rejection(
            candidateX, candidateY, maxROIs, nCrops, cropSize,
            filterImageID, rejectNearEmpty, minStdDev
        );
    }
    else if (strategy == "random+filtered") {
        selectedIndices = strategy_random_filtered(
            candidateX, candidateY, maxROIs, nCrops, cropSize,
            filterImageID, minMeanIntensity, rejectNearEmpty, minStdDev
        );
    }
    else if (strategy == "intensity_based") {
        selectedIndices = strategy_intensity_based(
            candidateX, candidateY, maxROIs, nCrops, cropSize,
            filterImageID, topFraction, rejectNearEmpty, minStdDev
        );
    }
    else {
        if (useMIPforFiltering && filterImageID != workID) {
            selectImage(filterImageID);
            close();
        }

        selectImage(workID);
        close();

        selectImage(sourceID);
        close();

        exit("Unknown strategy: " + strategy);
    }

    if (selectedIndices.length == 0) {
        print("No valid crops found for: " + filename);
        roiManager("Reset");

        if (useMIPforFiltering && filterImageID != workID) {
            selectImage(filterImageID);
            close();
        }

        selectImage(workID);
        close();

        selectImage(sourceID);
        close();
        continue;
    }

    // ----------------------------------------
    // CROP + SAVE
    // Export ONLY the selected channel working stack
    // ----------------------------------------
    for (i=0; i<selectedIndices.length; i++) {

        selectImage(workID);

        roiManager("Deselect");
        roiManager("Select", selectedIndices[i]);

        roiIndexForName = i + 1;

        if (hdf5ReadyNaming)
            cropName = basename + "_" + IJ.pad(roiIndexForName, 3);
        else
            cropName = basename + "_" + roiIndexForName;

        run("Duplicate...", "title=[" + cropName + "] duplicate");
        cropID = getImageID();
        cropTitle = getTitle();

        if (outputFormat == "TIFF") {
            saveAs("Tiff", outputDir + "/" + cropName + ".tif");
        }
        else if (outputFormat == "HDF5") {
            exporth5(cropTitle, cropName, outputDir, compressionLevel);
        }
        else {
            selectImage(cropID);
            close();

            if (useMIPforFiltering && filterImageID != workID) {
                selectImage(filterImageID);
                close();
            }

            selectImage(workID);
            close();

            selectImage(sourceID);
            close();

            exit("Unknown output format: " + outputFormat);
        }

        selectImage(cropID);
        close();
    }

    // ----------------------------------------
    // CLEANUP
    // ----------------------------------------
    roiManager("Reset");

    if (useMIPforFiltering && filterImageID != workID) {
        selectImage(filterImageID);
        close();
    }

    selectImage(workID);
    close();

    selectImage(sourceID);
    close();

    print("Saved " + selectedIndices.length + " crops for " + filename);
}

setBatchMode(false);
print("Finished processing all images.");



// ------------------------------------------------------------
// STRATEGY 1: RANDOM
// ------------------------------------------------------------
function strategy_random_with_optional_empty_rejection(candidateX, candidateY, total, n, cropSize, measureImageID, rejectNearEmpty, minStdDev) {

    valid = collectValidIndices(candidateX, candidateY, total, cropSize, measureImageID, -1, rejectNearEmpty, minStdDev);
    validCount = valid.length;

    if (validCount == 0)
        return newArray(0);

    if (n > validCount)
        n = validCount;

    randOrder = pickNUniqueRandomIndices(validCount, n);

    out = newArray(n);
    for (k=0; k<n; k++)
        out[k] = valid[randOrder[k]];

    return out;
}



// ------------------------------------------------------------
// STRATEGY 2: RANDOM + FILTERED
// ------------------------------------------------------------
function strategy_random_filtered(candidateX, candidateY, total, n, cropSize, measureImageID, minMean, rejectNearEmpty, minStdDev) {

    valid = collectValidIndices(candidateX, candidateY, total, cropSize, measureImageID, minMean, rejectNearEmpty, minStdDev);
    validCount = valid.length;

    if (validCount == 0)
        return newArray(0);

    if (n > validCount)
        n = validCount;

    randOrder = pickNUniqueRandomIndices(validCount, n);

    out = newArray(n);
    for (k=0; k<n; k++)
        out[k] = valid[randOrder[k]];

    return out;
}



// ------------------------------------------------------------
// STRATEGY 3: INTENSITY BASED
// ------------------------------------------------------------
function strategy_intensity_based(candidateX, candidateY, total, n, cropSize, measureImageID, topFraction, rejectNearEmpty, minStdDev) {

    if (topFraction <= 0 || topFraction > 1)
        topFraction = 0.30;

    valid = newArray(0);
    means = newArray(0);

    count = 0;
    for (i=0; i<total; i++) {
        stats = getROIstats(candidateX[i], candidateY[i], cropSize, measureImageID);
        meanVal = stats[0];
        stdVal  = stats[1];

        passEmpty = true;
        if (rejectNearEmpty && stdVal < minStdDev)
            passEmpty = false;

        if (passEmpty) {
            valid = Array.concat(valid, i);
            means = Array.concat(means, meanVal);
            count++;
        }
    }

    if (count == 0)
        return newArray(0);

    sortIndicesByValuesDescending(valid, means);

    keepN = floor(count * topFraction);
    if (keepN < 1)
        keepN = 1;

    if (n > keepN)
        n = keepN;

    brightest = newArray(keepN);
    for (i=0; i<keepN; i++)
        brightest[i] = valid[i];

    randOrder = pickNUniqueRandomIndices(keepN, n);

    out = newArray(n);
    for (k=0; k<n; k++)
        out[k] = brightest[randOrder[k]];

    return out;
}



// ------------------------------------------------------------
// Collect valid candidate indices
// If minMean < 0, mean filtering is ignored
// ------------------------------------------------------------
function collectValidIndices(candidateX, candidateY, total, cropSize, measureImageID, minMean, rejectNearEmpty, minStdDev) {

    out = newArray(0);

    for (i=0; i<total; i++) {
        stats = getROIstats(candidateX[i], candidateY[i], cropSize, measureImageID);
        meanVal = stats[0];
        stdVal  = stats[1];

        passMean = true;
        passEmpty = true;

        if (minMean >= 0 && meanVal < minMean)
            passMean = false;

        if (rejectNearEmpty && stdVal < minStdDev)
            passEmpty = false;

        if (passMean && passEmpty)
            out = Array.concat(out, i);
    }

    return out;
}



// ------------------------------------------------------------
// Return [mean, stdDev] for ROI on measurement image
// ------------------------------------------------------------
function getROIstats(x, y, size, measureImageID) {

    selectImage(measureImageID);
    makeRectangle(x, y, size, size);
    getStatistics(area, mean, min, max, std);

    stats = newArray(2);
    stats[0] = mean;
    stats[1] = std;
    return stats;
}



// ------------------------------------------------------------
// Sort indices by values descending
// ------------------------------------------------------------
function sortIndicesByValuesDescending(indices, values) {

    n = values.length;

    for (i=0; i<n-1; i++) {
        for (j=0; j<n-1-i; j++) {
            if (values[j] < values[j+1]) {

                tmpVal = values[j];
                values[j] = values[j+1];
                values[j+1] = tmpVal;

                tmpIdx = indices[j];
                indices[j] = indices[j+1];
                indices[j+1] = tmpIdx;
            }
        }
    }
}



// ------------------------------------------------------------
// Pick n unique random integers from [0 .. total-1]
// ------------------------------------------------------------
function pickNUniqueRandomIndices(total, n) {

    indices = newArray(total);

    for (i=0; i<total; i++)
        indices[i] = i;

    for (i=0; i<n; i++) {
        j = i + floor(random() * (total - i));

        tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }

    out = newArray(n);
    for (i=0; i<n; i++)
        out[i] = indices[i];

    return out;
}



// ------------------------------------------------------------
// USER PROVIDED FUNCTIONS
// ------------------------------------------------------------

// Extract a string from another string at the given input smaller string (eg ".")
function getBasename(filename, SubString){
  dotIndex = indexOf(filename, SubString);
  basename = substring(filename, 0, dotIndex);
  return basename;
}

// Return a file list contained in directory dir filtered by extension.
function getFilesList(dir, fileExtension) {  
  tmplist=getFileList(dir);
  list = newArray(0);
  imageNr=0;
  for (i=0; i<tmplist.length; i++)
  {
    if (endsWith(tmplist[i], fileExtension)==true || endsWith(tmplist[i], fileExtension + ".roi")==true)
    {
      list[imageNr]=tmplist[i];
      imageNr=imageNr+1;
    }
  }
  Array.sort(list);
  return list;
}

// === Function: scalingXYZ_multichannel ===
// Scales a multichannel 3D stack with different factors for XY and Z.
// Works by splitting channels, scaling individually, and merging back.
// scaleXY : scaling factor for X and Y
// scaleZ  : scaling factor for Z
// interp  : true = Bicubic interpolation, false = None
function scalingXYZ_multichannel(image_name, scaleXY, scaleZ, interp) {
	selectWindow(image_name);
	getDimensions(width, height, channels, slices, frames);
	originalDimensions = newArray(width, height, channels, slices, frames);

	// Determine interpolation mode (ImageJ macro syntax)
	if (interp) {
		interpMode = "Bicubic";
	} else {
		interpMode = "None";
	}

	// Compute target dimensions for Scale...
	numSlices = nSlices / channels;
	scaledX = width * scaleXY;
	scaledY = height * scaleXY;
	scaledZ = numSlices * scaleZ;

	// --- Single channel case ---
	if (channels == 1) {
		selectWindow(image_name);
		run("Scale...", 
			"x=" + scaleXY +
			" y=" + scaleXY +
			" z=" + scaleZ +
			" width=" + scaledX +
			" height=" + scaledY +
			" depth=" + scaledZ +
			" interpolation=" + interpMode +
			" average process create");
		return originalDimensions;
	}

	// --- Multi-channel case ---
	scaledChannelIDs = newArray(channels);

	for (c = 1; c <= channels; c++) {
		// Select original image and duplicate current channel only
		selectWindow(image_name);
		run("Duplicate...", "duplicate channels=" + c + "-" + c);
		rename("temp_channel_" + c);

		// Scale current channel
		run("Scale...", 
			"x=" + scaleXY +
			" y=" + scaleXY +
			" z=" + scaleZ +
			" width=" + scaledX +
			" height=" + scaledY +
			" depth=" + scaledZ +
			" interpolation=" + interpMode +
			" average process create");

		// Rename scaled channel image to its channel index for merging
		rename("scaled_channel_" + c);
		scaledChannelIDs[c - 1] = getImageID();

		// Close temporary duplicate
		close("temp_channel_" + c);
	}

	// --- Merge all scaled channels ---
	mergeCommand = "";
	for (c = 1; c <= channels; c++) {
		mergeCommand += " c" + c + "=" + "scaled_channel_" + c;
	}
	mergeCommand += " create";
	run("Merge Channels...", mergeCommand);

	return originalDimensions;
}

// ------------------------------------------------------------
// HDF5 EXPORT
// ------------------------------------------------------------
function exporth5(image, output_name, output_dir, compression) { 
    export_path = output_dir + File.separator + output_name + ".h5";
    dataset = "/data";
    run("Export HDF5", "input=[" + image + "] exportpath=[" + export_path + "] datasetname=[" + dataset + "] compressionlevel=" + compression);
}