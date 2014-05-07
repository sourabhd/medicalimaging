

filename = '../data/slicer3d/MRHead.nrrd';
[X, meta] = nrrdread(filename);

for i = 1:130
imagesc(X(:,:,i))
pause(1);
end