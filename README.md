# mmqt
**Microglia morphology quantification tool (MMQT)**

Toolbox for the extraction of morphological features of microglia cells from confocal
microscopy image stacks containing two channels:
* **anti-Iba1** staining of microglia
* **DAPI** staining of nuclei

The software will segment foreground from background, create a skeleton of the microglia,
segregate individual microglia cells and extract morphological features.

For details see the method paper _[Reference pending]_.  
Please read the LICENSE file before using MMQT.

### Disclaimer

**This is a preliminary upload of the MMQT scripts for peer review purposes only.
More detailed information and instructions will follow after publication of the manuscript.**

### Content

* Hardware and software requirements
* Getting started
* Accepted data formats
* Output files

### Hardware and software requirements

* Modern CPU, e.g. Intel Core i7
* at least 16 GB RAM
* Obligatory:
  * MATLAB® R2016B ([The Mathworks](https://www.mathworks.com), 
    newer versions untested but might work) with the toolboxes:
  * Image Processing Toolbox 9.5
  * Statistics and Machine Learning Toolbox 11.0
* When using Carl Zeiss Image Data File (CZI-file; file-extension .czi) as input::  
  Bio-Formats software tool for MATLAB available from the 
  [Open Microscopy Environment consortium](http://www.openmicroscopy.org).

### Getting started

* Add the MMQT folder to the Matlab path
* Try one of the example scripts in the subfolder `./example`
* Example CZI data can be downloaded from the 
  [ISD Research Server](http://download.isd-muc.de/mmqt/sample_czi.zip)
	
### Accepted data formats

The MMQT toolbox needs confocal microscopy Z-stack image data with two channels:
* Channel **1**: DAPI staining of nuclei
* Channel **2**: anti-Iba1 staining of microglia

Images can be provided in two file formats:

* **Carl Zeiss Image Data File** (CZI-file; file-extension .czi).  
  For details on how to import them see example script `mmqt_example_CZI.m`.

* **MATLAB® formatted binary file** (MAT-file; file-extension: .mat).  
  MAT-files have to contain two variables:
  * `img`: 4D array containing the raw image data with 8 bit integer encoding; dimensions
	should correspond to X/Y/Z directions and the color channels.
  * `hdr`: structure with two fields, describing dimensionality and resolution (given as 
	scaling factors) of the raw data in the array "img":
    * `hdr.dim`: vector of length 4, with values specifying the number of dimensions
	  in array "img" and the size of each dimension.  
	  For example, `hdr.dim=[4,1024,1024,132,2]` specifies that the array "img" has 
	  4 dimensions with the size of dimensions being respectively 1024, 1024, 132 and 2.
    * `hdr.pixdim`: vector of length 4, with values specifying the number of dimensions in 
      array "img" and the factors to scale each dimension to micrometer.	The forth 
      dimension corresponds to the color channels and therefore has no scaling applied to
      it. Therefore, set its scaling factor to 1.  
      For example, `hdr.pixdim=[4,0.2076,0.2076,0.4,1]` specifies that the array "img" has 
      4 dimensions and the X/Y dimension has to be scaled by 0.2076 and the Z dimension 
      by 0.4.

Please note that scaling information for the spatial dimensions is crucial for correct 
processing of the image stacks. This scaling factors sould result in a micrometer scale.
The MMQT software will try to automatically extract this information from CZI-files.
If you work with MAT-files, make sure that this information is set correctly in 
this file (see above).

### Output files

The output files creates by MMQT are documented in the following files:
`mmqt_output_files.xlsx` (MS Excel file) and `mmqt_output_files.txt` (tab-separated text)

The scores for the extracted shape-features are stored in a tab-separated text-file 
with the suffix `_features_summarized.txt`  
The ID of a cell, as given in this text-file, corresponds to the value of the voxels, 
which belong to the same cell, in the volumetric dataset being saved in the file with 
suffix `_stack5_segregated_cells.mat`.

