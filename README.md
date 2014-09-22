Implementation of Chan-Vese like segmentation algorithm for 3D medical images
=============================================================================

Authors:
-------
* Ahsanul Karim <karim.ahsanul@gmail.com>
* Sourabh Daptardar <saurabh.daptardar@gmail.com>

License:
-------
BSD


About the code:
--------------
This software implements Chan-Vese like segmentation algorithm based upon level
sets. Major differences are : 1. It uses 'gaussian' curvature, instead of mean
curvature 2. Targeted for 3D objects rather than 2D images.

This code was created for course project in medical imaging course 
by Prof. Allen Tanenbaum at Stony Brook University.
https://www.cs.stonybrook.edu/people/faculty/AllenTannenbaum

Requires:
--------
MATLAB

Steps to run
------------
On Matlab prompt:
cd src/
main

main.m contains demos to run the Slicer3D sample data files as well as
synthetic 3D dataset like 3D blobs, 3D solid sphere.


Sample data
-----------
Slicer3D sample data (MR images)  is placed in the data/ directory
This data was obtained from Slicer 3D app :  http://www.slicer.org/
Refer to the website for terms of usage.

Demos
-----
* Blob: https://www.youtube.com/watch?v=lxOBL_5ZKCA
* MR Head: https://www.youtube.com/watch?v=tsz1UCLJ4qw

Also shipped along with the following third party codes which the software uses:
-------------------------------------------------------------------------------
1. nrrdread.m : Read NRRD images into matlab
http://www.mathworks.com/matlabcentral/fileexchange/34653-nrrd-format-file-reader/content/nrrdread.m 
BSD license

2. vol3d.m : Render 3D graphics in matlab
http://www.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2
BSD license

References:
----------
1. Active Contours Without Edges, IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 10, NO. 2, FEBRUARY 2001
   Tony F. Chan, and Luminita A. Vese
   http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=902291

2. Summer Research School on Medical Imaging, 2012, The Fields Institute for Research in Mathematical Sciences.
   Todd Wittman
http://www.math.ucla.edu/~wittman/Fields/

3. For general overview on level set methods
http://www.museth.org/Ken/Publications_files/Breen-etal_SIG04.pdf
