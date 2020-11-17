/*
This macro is used to assist manual annotation of images in ImageJ.
ROIs are saved to a zip file and an overlay on the original image is stored.
Image filenames are maintained in the generated files.
When executing this macro on an image that has been processed, the already generated ROIs will be imported.

First open a single image, then run the macro.
*/

title = getTitle();
path = getDirectory("Image");
exportpath = path + "/export/";
roipath = path + "/export/" + title + "_roi.zip";

if (!File.exists(exportpath)) {
	File.makeDirectory(exportpath);
}
if (File.exists(roipath)) {
	roiManager("open", roipath);
	run("From ROI Manager");
	roiManager("reset");
}
roiManager("reset");
setTool("brush");

showMessage("Use key B to store ROI");
waitForUser("Press OK when done");

// Save ROIs
run("To ROI Manager");
roiManager("Show All with labels");
roiManager("save", roipath);

// Make Overlay
run("From ROI Manager");
roiManager("reset");
run("Overlay Options...", "stroke=yellow width=0.5 fill=none apply show");
Stack.setDisplayMode("composite");
run("Flatten");
saveAs("jpg", exportpath + title + "_overlay.jpg");
close();
