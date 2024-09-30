# ANTEMA

Automated Nanoparticle Transmission Electron Micrograph Analysis

This program is designed to automatically analyze high resolution transmission electron microscopic images (.dm3 files) of nanoparticles for shape and size related parameters by using a neural network for image segmentation. 
This program can be run in MATLAB by running ANTEMA_Main.m. Necessary inputs will be requested in the command window.
Or you can install a version with a graphical user interface as a MATLAB app.

You can test the program on the dm3 files in the testimages folder.

MATLAB Version R2021b
Required Toolboxes:
- Deep Learning Toolbox
- Image Processing Toolbox
- Computer Vision Toolbox
- Statistics and Machine Learning Toolbox

More detailed explanations can be found in the corresponding paper: N. Gumbiowski, K. Loza, M. Heggen, M. Epple, Nanoscale advances 2023, 2318 – 2326, DOI: 10.1039/d2na00781a

The version of the program referred to in this paper is v1.0.
If you use this program in your work, please cite the corresponding paper.

Since then the program has been updated in the following ways:
- a lower boundary for the convexity was implemented as well as a maximum number of holes a particle area is allowed to have, so that particles with too much overlap are excluded from evaluation and not fed into the separation routine
- tif images can now also be used as an input
- the option of using watershed instead of UECS for separation was implemented which can be useful for certain types of images
- a python version of the trained network has been uploaded
- a python script for the segmentation has been added (this python script is still very much incomplete, the only part that is working and fully implemented is the preprocessing of the image and the subsequent segmentation using the model, all post processing steps have not been implemented yet)

For further information contact: matthias.epple@uni-due.de
