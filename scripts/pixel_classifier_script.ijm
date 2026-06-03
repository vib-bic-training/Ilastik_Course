// @ String(visibility=MESSAGE, value="File location and format", required=false) msg1
// @ File(label="File directory", style="Directory") dir
// @ String (label="File suffix", choices={".czi", ".nd2", ".tif"}, style="listBox") suffix

// @ String(visibility=MESSAGE, value="Pixel Classifier Parameters", required=false) msg2
// @ File(label="Pixel classifier path", style="file") pixel_class
// @ Integer (label="Channel to segment. Use 0 for all channels", value=2, min=0) channel
// @ String (label="File suffix", choices={"Probabilities", "Segmentation"}, style="listBox") output
// @ Float(label="Probability threshold for segmentation (if probabilities output was selected)", value=0.5) proba_thresh

// @ String(visibility=MESSAGE, value="Pre-processing parameters prior classification", required=false) msg3
// @ Boolean(label="Normalize intensity?", value=false) normalize
// @ String(label="Normalization parameters. Paste from recorder.", value="saturated=0.1 normalize process_all use") normalize_param

// @ Float(label="XY scaling factor", value=1.0) scaleXY
// @ Float(label="Z scaling factor", value=1.0) scaleZ
// @ Boolean(label="Use interpolation", value=false) interpolation

// @ String(visibility=MESSAGE, value="Pixel classifier validation", required=false) msg4
// @ Boolean(label="Test pixel classifier? (Ground truth necessary)", value=false) do_validation
// @ File(label="File directory", style="Directory", required=false) groundTruth_dir
// @ String(label="Ground truth files suffix", value="_groundTruth.tif") groundTruth_suffix


/* 29/05/2026
Nicolas Peredo, PhD
Image Analysis Expert
VIB BioImaging Core Leuven, VIB Technologies, Center for Brain and Disease Research, Leuven, Belgium
VIB BioImaging Core Leuven, Department of Neurosciences, KU Leuven, Leuven, Belgium
Nikon Center of Excellence
Campus Gasthuisberg - ON5 - room 04.367
Herestraat 49 - box 62
3000 Leuven
Belgium
phone +32 (0)16/37.70.03

When you publish data analyzed with this script please add the references of the used plug-ins:
Ilastik plug-in
Ilastik itself
Morpholibj
*/


//File directory
fileList = getFilesList(dir, suffix);
Array.sort(fileList);

//Create the different folders with results
File.makeDirectory(dir + File.separator + "processed");

//Arrays for table
jaccard_array = newArray(fileList.length);

for (files = 0; files < fileList.length; files++) {

	//Getting file names
	file = fileList[files];
	basename = getBasename(file, suffix);
	
	//Open image
	run("Bio-Formats Importer", "open=[" + dir + File.separator + file + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
	rename("raw_image");
	
	//Create individual channels
	run("Duplicate...", "duplicate channels=" + channel + "-" + channel);
	rename("channel_forSegmentation");
	
	if (normalize) {
		normalizeChannels("channel_forSegmentation",normalize_param);
		rename("normalized");
	}
	
	if (scaleXY != 1.0 || scaleZ != 1.0) {
		originalDimensions = scalingXYZ_multichannel(image_name, scaleXY, scaleZ, interpolation);
		rename("scaled");
	}

	rename("to_process");
	
	// Apply classifier to segment the signal
	runIlastikClassification("to_process", pixel_class, output);
	rename("processed");
	
	//// Generate binary from probabilities image
	if (output == "Probabilities") {
		run("Duplicate...", "duplicate channels=1-1");
		rename("temp");
		setThreshold(proba_thresh, 1000000000000000000000000000000.0000);
		run("Convert to Mask", "background=Dark black create");
		rename("segmented");
		close("temp");
	}
	
	// Generate binary from segmentation image
	else {
		setThreshold(1, 1, "raw");
		run("Convert to Mask", "background=Dark black create");
		rename("segmented");
	}
	
	if (do_validation) {
		// Open ground truth image
		run("Bio-Formats Importer", "open=[" + groundTruth_dir + File.separator + basename  + groundTruth_suffix + "] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1");
		rename("ground_truth");
	    
		run("Label Overlap Measures", "source=segmented target=ground_truth jaccard");
		
		jaccard_array[files] = getResult("JaccardIndex", 0);
	}
	
	saveAs("Tiff", dir + File.separator + "processed" + File.separator + basename + ".tif");
	
	close("*");
	
}

if (do_validation) {
	Table.create("Validation_table");
	Table.setColumn("Filename", fileList);
	Table.setColumn("Jaccard Index", jaccard_array);
	
	File.makeDirectory(dir + File.separator + "validation");
	saveAs("Results", dir + File.separator + "validation" + File.separator + "segmentation_validation.csv");
}

//Extract a string from another string at the given input smaller string (eg ".")
function getBasename(filename, SubString){
  dotIndex = indexOf(filename, SubString);
  basename = substring(filename, 0, dotIndex);
  return basename;
}

//Return a file list contain in the directory dir filtered by extension.
function getFilesList(dir, fileExtension) {  
  tmplist=getFileList(dir);
  list = newArray(0);
  imageNr=0;
  for (i=0; i<tmplist.length; i++)
  {
    if (endsWith(tmplist[i], fileExtension)==true)
    {
      list[imageNr]=tmplist[i];
      imageNr=imageNr+1;
      //print(tmplist[i]);
    }
  }
  Array.sort(list);
  return list;
}

// === Function: normalizeChannels_split ===
// Normalizes intensity of each channel using Enhance Contrast
// by splitting channels (more memory efficient than duplicating).
// Works for both single- and multi-channel stacks.
function normalizeChannels(image_name,parameters) {
	selectWindow(image_name);
	name = getTitle();
	getDimensions(width, height, channels, slices, frames);

	// --- Single channel case ---
	if (channels == 1) {
		selectWindow(image_name);
		run("Enhance Contrast...", parameters);
		return;
	}

	// --- Multi-channel case ---
	selectWindow(image_name);
	run("Split Channels");  // Creates separate images for each channel

	// Build array of split channel titles (ImageJ names them automatically)
	channelTitles = newArray(channels);
	for (c = 1; c <= channels; c++) {
		channelTitles[c - 1] = "C" + c + "-" + name;
	}

	// Apply normalization to each split channel
	for (c = 0; c < channels; c++) {
		selectWindow(channelTitles[c]);
		run("Enhance Contrast...", parameters);
	}

	// Merge all channels back
	mergeCommand = "";
	for (c = 1; c <= channels; c++) {
		mergeCommand += " c" + c + "=" + channelTitles[c - 1];
	}
	mergeCommand += " create";
	run("Merge Channels...", mergeCommand);

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

// function description
function process_file(image_name) { 
	selectWindow(image_name);
	originalDimensions = scalingXYZ_multichannel(image_name, 0.6, 0.6, true);
	rename("processed");
	
	return originalDimensions;
}

// === Function: runIlastikPrediction ===
// Runs Ilastik Pixel Classification Prediction on the current image
// using the given .ilp project and input .h5 file path.
// After prediction, restores original voxel calibration.
//output (string): can be either "Probabilities" or "Segmentation"
function runIlastikClassification(image_name, pixel_class,output) {
	// --- Get original calibration ---
	getVoxelSize(width, height, depth, unit);
	xy_size = width;
	z_size = depth;

	// --- Run Ilastik prediction ---
	run("Run Pixel Classification Prediction", 
		"projectfilename=" + pixel_class + 
		" inputimage=" + image_name + 
		" pixelclassificationtype=" + output);
		
	getDimensions(width, height, channels, slices, frames);

	// --- Restore voxel calibration ---
	Stack.setXUnit(unit);
	run("Properties...", 
		"channels=" + channels + " slices=" + slices + " frames=1 " +
		"pixel_width=" + xy_size + " pixel_height=" + xy_size + " voxel_depth=" + z_size);
	
}


// 3D rescalling with or without interpolation. 
// This functions rescales scaled images to their original size. 
// This function takes as an argument the originaldimensions variable from the scaling3D() function.
// Interp =  interpolation, write either true or false
function rescaling3D(orig_dim, interp){
	if (interp) {
		run("Scale...", "x=- y=- z=- width=" + orig_dim[0] + " height=" + orig_dim[1] + " depth=" + orig_dim[3] + " interpolation=Bicubic average process create");
	}
	else {
		run("Scale...", "x=- y=- z=- width=" + orig_dim[0] + " height=" + orig_dim[1] + " depth=" + orig_dim[3] + " interpolation=None process create");
	}
}

// === Function: addChannelToStack ===
// Adds a new single-channel image as an extra channel to the current multichannel stack.
// multichannel_stack : the title (window name) of the multichannel image stack.
// new_image_title : the title (window name) of the image to add.
function addChannelToStack(multichannel_stack, new_channel) {
	// Get current multichannel stack info
	selectWindow(multichannel_stack);
	getDimensions(width, height, channels, slices, frames);

	// --- Single channel case ---
	if (channels == 1) {
		// If original is single channel, just merge with the new image directly
		mergeCommand = "c1=" + getTitle() + " c2=" + new_channel + " create";
		run("Merge Channels...", mergeCommand);
		return;
	}

	// --- Multi-channel case ---
	// Split original multichannel stack
	selectWindow(multichannel_stack);
	run("Split Channels");

	// Build array of split channel titles
	channelTitles = newArray(channels);
	for (c = 1; c <= channels; c++) {
		channelTitles[c - 1] = "C" + c + "-" + multichannel_stack;
	}

	// Build the Merge Channels command including the new image as the last channel
	mergeCommand = "";
	for (c = 1; c <= channels; c++) {
		mergeCommand += " c" + c + "=" + channelTitles[c - 1];
	}
	mergeCommand += " c" + (channels + 1) + "=" + new_channel;
	//mergeCommand += " create";

	run("Merge Channels...", mergeCommand);

	// Optional: close split channel windows to clean up
	/*
	for (c = 0; c < channels; c++) {
		selectWindow(channelTitles[c]);
		close();
	}
	*/
}

/*
 * Add a constant-value column to a table.
 *
 * Arguments:
 * - tableName  : name of the table to modify
 * - columnName : name of the new column
 * - value      : value to add to every row
 */
function addConstantColumnToTable(tableName, columnName, value) {

    nRows = Table.size(tableName);

    for (row = 0; row < nRows; row++) {
        Table.set(columnName, row, value, tableName);
    }

    Table.update(tableName);
}