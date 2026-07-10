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
Step 1: Load the experimental XRD dataset: The software imports the diffraction pattern and displays the raw experimental data. If junk regions are present, the user can optionally exclude them by clicking on them.

<img width="892" height="485" alt="img_1" src="https://github.com/user-attachments/assets/f504faab-69f4-44ad-9f16-dc1ed1dbe9cc" />


Step 2: Review automatically detected peaks: The algorithm estimates the background, detects candidate diffraction peaks, performs a Voigt fit, and labels every fitted peak. If an incorrect peak is detected, the user can click it to remove it.


<img width="855" height="465" alt="img_4" src="https://github.com/user-attachments/assets/ebb1aaa5-9c32-4240-ab52-76c90e5e4328" />

Step 3: Add missing peaks (optional): If the algorithm misses a peak, the user can click its approximate position. The optimizer then refits all peaks simultaneously.

<img width="858" height="467" alt="img_5" src="https://github.com/user-attachments/assets/074a8b89-1848-4dc2-a26f-9e4418a224ef" />


Step 4: Final optimized fit

<img width="867" height="471" alt="img_6" src="https://github.com/user-attachments/assets/945b033e-0c14-47d2-957d-a977176e3fc2" />

Black dots: Experimental data

Red line: Total fitted model

Green dashed line: Voigt peak model

Blue dashed line: Estimated background

Lower panel: Residual (difference between experiment and model)

A small, randomly distributed residual indicates that the mathematical model closely represents the experimental data.


Step 5: Export results

<img width="702" height="530" alt="img_7" src="https://github.com/user-attachments/assets/95fdf2be-2e2f-4053-9e97-577b70836c87" />

## Copyright

Copyright © 2026 Shashwat Agarwal.

This repository is provided for viewing and educational purposes only.
No license is granted for copying, modifying, redistributing, or commercial use without prior written permission.
