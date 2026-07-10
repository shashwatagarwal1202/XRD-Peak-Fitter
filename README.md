# XRD-Peak-Fitter
Interactive Python-based XRD peak-fitting tool featuring automatic background estimation, Voigt-profile fitting, residual analysis, manual refinement, and export of fitted peak parameters. Developed to automate diffraction data analysis during my internship at BARC.

## Features

- Automatic peak detection
- Manual peak addition
- Manual peak exclusion
- Voigt fitting
- Background estimation
- Batch processing
- Multiple export formats

## Example Fit
Step 1:

![img_1.png](img_1.png)
Load the experimental XRD dataset: The software imports the diffraction pattern and displays the raw experimental data. If junk regions are present, the user can optionally exclude them by clicking on them.

Step 2:

![img_4.png](img_4.png)
Review automatically detected peaks: The algorithm estimates the background, detects candidate diffraction peaks, performs a Voigt fit, and labels every fitted peak. If an incorrect peak is detected, the user can click it to remove it.

Step 3:

![img_5.png](img_5.png)
Add missing peaks (optional): If the algorithm misses a peak, the user can click its approximate position. The optimizer then refits all peaks simultaneously.

Step 4:

![img_6.png](img_6.png)
Final optimized fit

Black dots: Experimental data

Red line: Total fitted model

Green dashed line: Voigt peak model

Blue dashed line: Estimated background

Lower panel: Residual (difference between experiment and model)

A small, randomly distributed residual indicates that the mathematical model closely represents the experimental data.

Step 5:

![img_7.png](img_7.png)
Export results

## Copyright

Copyright © 2026 Shashwat Agarwal.

This repository is provided for viewing and educational purposes only.
No license is granted for copying, modifying, redistributing, or commercial use without prior written permission.
