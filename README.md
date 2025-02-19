# Dynamic B0 Shimming for ASL Matlab Code

Matlab code for calculating B0 shim settings within the labelling plane to optimise the labelling efficiency of pseudocontinuous arterial spin labelling, based on this paper:

- *Ji Y, Woods JG, Li H, Okell TW. Dynamic  B0  field shimming for improving pseudo‐continuous arterial spin labeling at 7 T. Magnetic Resonance in Med. 2025; 93: 1674–1689. [[DOI link](https://doi.org/10.1002/mrm.30387)]*

## Brief explanation
Pseudocontinuous arterial spin labeling (PCASL) is a non-invasive MRI method to generate perfusion images or angiograms. However, PCASL is sensitive to main magnetic field ($B_0$) inhomogeneities, which can compromise labelling efficiency, particularly at higher field strengths (e.g. 7T). The $B_0$ field is often only optimised for the brain, giving poor field homogeneity at the labelling plane in the neck. If we attempt to shim for both the brain and the labelling plane, the result is often poor homogeneity across both regions, compromising labelling efficiency and image quality. 

Instead, here we propose optimising the static shim for the brain but then dynamically modify the linear shim terms (gradients) between the labelling and imaging periods, to achieve good $B_0$ homogeneity in the neck during the PCASL pulse train, giving good labelling efficiency, and reverting to the brain-optimised shim values during imaging to maintain high image quality. Since we only care about the $B_0$ homogeneity within the vessels during labelling, we can use small ROIs to optimise the shim terms in the neck. 

This Matlab code allows the reading in and reconstruction of raw k-space double-echo fieldmap data, interactive selection of vessel ROIs and outputs the required linear shim and global frequency offset terms required for optimal correction during the PCASL labelling period.

## Code structure
The main dynamic shim calculation is performed in `run.m`, which calls a number of other functions in the `functions` folder.  A few options are available that can be set in `run.m`:
- `roi_exist`: If set to true, the code will look for a file named `mask_shim.mat` and load this to use as the vessel mask.
- `autoroi`: If set to true, the user is prompted to interactively select the central point of each vessel from the central slice. A cylindrical ROI is then created around this point with radius defined by `r_mm`. To avoid include noisy voxels (e.g. air), the mask is further thresholded at `Frac` * 95th percentile signal intensity within the initial mask.
- If both of the above are set to false, the user is prompted to draw a manual polygon ROI over each vessel within each slice.

## External dependencies
Tested using MATLAB 2024b. We rely on the excellent `mapVBVD` by Philipp Ehses and colleagues. The original version can be found here: https://github.com/pehses/mapVBVD. 

## Contributors
Original code: Yang Ji: https://orcid.org/0000-0003-4134-374X

Updates: Thomas W. Okell: https://orcid.org/0000-0001-8258-0659


## How to cite

Please cite this paper:

- *Ji Y, Woods JG, Li H, Okell TW. Dynamic  B0  field shimming for improving pseudo‐continuous arterial spin labeling at 7 T. Magnetic Resonance in Med. 2025; 93: 1674–1689. [[DOI link](https://doi.org/10.1002/mrm.30387)]*

## Copyright

Copyright © 2025, University of Oxford

