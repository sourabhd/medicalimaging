
clear all; close all; clc; imtool close all;
dbstop if error;

% filename = 'mrtumor.jpg';
% I = read_image(filename);
% %curve_evolution(I);
% KG = gaussian_curvature(I);
% size(I)

filename = '../data/slicer3d/MRHead.nrrd';
[X, meta] = nrrdread(filename);
h = vol3d('cdata',X,'texture','3D');
view(3); 
axis tight;  daspect([1 1 1]);
alphamap('rampup');
alphamap(.06 .* alphamap);
%alphamap('default');
%alphamap(.01 .* alphamap);
%set(gcf,'DoubleBuffer','on')
%hold on; contour3(F); hold off;
%f2 = figure;
cv_3d(X,0.2);