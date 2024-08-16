target = "10Nov_TAK1"; // CHANGE for new group
z_plane = 1; // CHANGE each time

run("Clear Results");
run("Crop");
run("8-bit");
rename("FirstStack");
run("Duplicate...", "duplicate");
rename("ASubstack");
selectImage("ASubstack");
run("Invert LUT"); // Disable if using pERK
run("Subtract Background...", "rolling=50");
run("Smooth");
run("Gaussian Blur...", "sigma=2");
run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None*");
run("Enhance Contrast...", "saturated=0.4");
run("Find Maxima...", "prominence=5 output=[Maxima Within Tolerance]");
setOption("BlackBackground", true);
run("Erode");
run("Watershed");
run("Set Measurements...", "area mean standard centroid stack nan redirect=None decimal=3");
run("Analyze Particles...", "  circularity=0.05-1.00 show=Outlines exclude clear summarize add");
selectImage("Drawing of ASubstack Maxima");
selectImage("FirstStack");
npoints = RoiManager.size;
k = 1;
for (i=0; i<npoints; i +=k) {;
roiManager("Select", i);
roiManager("Measure");
}

dir = "/Users/kowteckfong/Desktop/MP InVivoCaImaging/pERK data/TAK/" + target + "/"; // CHANGE for new group
saveAs("Results",  dir + target + "_" + z_plane + ".csv") // Add pERK.csv if pERK

close("ASubstack");
close("ROI Manager");
close("Summary");
close("ASubstack Maxima");
close("FirstStack");
close("Drawing of ASubstack Maxima");