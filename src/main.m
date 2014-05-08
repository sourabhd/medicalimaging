
clear all; close all; clc; imtool close all;
dbstop if error;


T = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run on 3D slicer data

%filename = '../data/MRHead.nrrd';
%filename = '../data/MRBrainTumor1.nrrd';
%filename = '../data/MRBrainTumor2.nrrd';
%cv_3d(filename,'yes',T);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Run on synthetic data

% Test Case 1 : Solid Sphere 
%blob('solid',T);

% Test Case 2 : Blob
blob('blob',T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%