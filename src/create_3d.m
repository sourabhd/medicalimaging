
clear all; close all; clc;

filename = '../data/MRHead.nrrd';
[D, meta] = nrrdread(filename);


%load mri
%D = squeeze(D);
%D(:,1:60,:) = [];
% p = patch(isosurface(D, 5), 'FaceColor', 'none', 'EdgeColor', 'none');
% p2 = patch(isocaps(D, 5), 'FaceColor', 'interp', 'EdgeColor', 'none');
% view(3); axis tight;  daspect([1 1 1])
% colormap(gray(100))
% camlight; lighting gouraud
% isonormals(D, p);

%load mri.mat
%D = squeeze(D);
h = vol3d('cdata',D,'texture','3D');
view(3);  
axis tight;  daspect([1 1 1]);
%alphamap('rampup');
%alphamap(.06 .* alphamap);
alphamap('default');
alphamap(.01 .* alphamap);



% for i = 1:130
% imagesc(X(:,:,i))
% pause(1);
% end