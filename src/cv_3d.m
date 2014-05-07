%function [F] = cv_3d(f, lambda)

% Chan-Vese Active Contours Segmentation using Level Sets
% Input: Original grayscale image f
%        Parameter lambda (Try lambda=0.1)
% Output: Level Set Function Phi

%Set default value for lambda if user does not provide lambda.
clear all; close all; clc; 
imtool close all;
dbstop if error;

% filename = 'mrtumor.jpg';
% I = read_image(filename);
% %curve_evolution(I);
% KG = gaussian_curvature(I);
% size(I)

%filename = '../data/MRHead.nrrd';
filename = '../data/MRBrainTumor1.nrrd';
[X, meta] = nrrdread(filename);

X([X>255]) = 255;
h = vol3d('cdata',X,'texture','3D');
view(3); 
axis tight;  daspect([1 1 1]);
alphamap('default');
alphamap(.1 .* alphamap);

%X(1:256,1:256,1:130) = 0;
%X(25:75,25:75,25:75) = 240;
%[xx yy zz] = meshgrid(1:100,1:100,1:100);
%X = sqrt((xx-50).^2+(yy-50).^2+(zz-50).^2)<=25;
%alphamap(.1 .* alphamap);

%if ~exist('lambda')
lambda = 0.1;
%end;

%Parameters
T = 10;    %Stopping time
dt = 0.2;  %Time step
a = 0.01;  %Fudge factor to avoid division by zero.
epsilon = 0.1;  %Epsilon in delta approximation
%Fidelity terms lambda_in & lambda_out could be different.
%Both set = lambda here for simplicity.
lambda_in = lambda;
lambda_out = lambda;  

%I = im2double(X);
I = double(X);
I = I - min(I(:));
I = I / max(I(:));
X = double(X);
%F = im2double(f);
[m,n,p,k] = size(X);
%I = double(f);

%Initialize with large box around perimeter of image.
F(1:m,1:n,1:p) = -1;  
%F(2:m-1,2:n-1,2:p-1) = 0;  
F(100:200,100:200,50:100) = 1;
 

%F(1,:,:) = 0;
%F(:,1,:) = 0;
%F(:,:,1) = 0;
%F(m,:,:) = 0;
%F(:,n,:) = 0;
%F(:,:,p) = 0;

figure;
for t = 0:dt:T
    %Approximate derivatives
    F_x = (F([2:m,m],:,:) - F(1:m,:,:))/2;
    F_y = (F(:,[2:n,n],:) - F(:,1:n,:))/2;
    F_z = (F(:,:,[2:p,p]) - F(:,:,1:p))/2;
    
    F_xx = (F_x([2:m,m],:,:) - F_x(1:m,:,:))/2;
    F_yy = (F_y(:,[2:n,n],:) - F_y(:,1:n,:))/2;
    F_zz = (F_z(:,:,[2:p,p]) - F_z(:,:,1:p))/2;
    
    F_xy = (F_x(:,[2:n,n],:) - F_x(:,1:n,:))/2; 
    F_yz = (F_y(:,:,[2:p,p]) - F_y(:,:,1:p))/2;
    F_xz = (F_x(:,:,[2:p,p]) - F_x(:,:,1:p))/2;
    
    %TV term = Num/Den
    Num = (F_z.*(F_xx.*F_z- 2.*F_x.*F_xz) +(F_x.^2).* F_zz).*(F_z.*(F_yy.*F_z- 2.*F_y.*F_yz) +(F_y.^2) .* F_zz)-((F_z.*(-1.*F_x.*F_yz- F_xy.*F_z - F_xz.*F_y) +F_x.* F_y.* F_zz).^2);
    Den = F_z.^2.*(F_x.^2+F_y.^2+F_z.^2) + a;
    
    c_in = sum(sum(sum([F>0].*X)))./(a+sum(sum(sum([F>0]))));
    c_out = sum(sum(sum([F<0].*X)))./(a+sum(sum(sum([F<0]))));
    
    F = F + dt*epsilon./(pi*(epsilon^2+F.^2)).*( Num./Den  - lambda_in*(X-c_in).^2 + lambda_out*(X-c_out).^2);
    
    %Plot results.  Note drawing every iteration slows down the program.
    
    hold on;
    %montage(F);
    %subplot(1,3,1); 
    %
    tic;pause(0.01);toc;
    %axis equal;
    subplot(3,3,1); imshow(I(:,:,30)); title('Level Set');
    hold on; contour(F(:,:,30),[0,0],'r'); hold off;
    subplot(3,3,2); imshow(I(:,:,40)); title('Level Set');
    hold on; contour(F(:,:,65),[0,0],'r'); hold off;
    subplot(3,3,3); imshow(I(:,:,50)); title('Level Set');
    hold on; contour(F(:,:,100),[0,0],'r'); hold off;
    subplot(3,3,4); imshow(I(:,:,30)); title('Level Set');
    hold on; contour(F(:,:,30),[0,0],'r'); hold off;
    subplot(3,3,5); imshow(I(:,:,60)); title('Level Set');
    hold on; contour(F(:,:,65),[0,0],'r'); hold off;
    subplot(3,3,6); imshow(I(:,:,70)); title('Level Set');
    hold on; contour(F(:,:,100),[0,0],'r'); hold off;
    subplot(3,3,7); imshow(I(:,:,80)); title('Level Set');
    hold on; contour(F(:,:,30),[0,0],'r'); hold off;
    subplot(3,3,8); imshow(I(:,:,90)); title('Level Set');
    hold on; contour(F(:,:,65),[0,0],'r'); hold off;
    subplot(3,3,9); imshow(I(:,:,100)); title('Level Set');
    hold on; contour(F(:,:,100),[0,0],'r'); hold off;
    
    %subplot(2,3,4); imshow(F(:,:,30)); title('Level Set');
    
    %subplot(2,3,5); imshow(F(:,:,65)); title('Level Set');
    
    %subplot(2,3,6); imshow(F(:,:,100)); title('Level Set');
    
    %hold on; contour(F(:,:,50),[0,0],'r'); hold off;
    %hold on;
    %A = 255.*[F(:,:,50)> 0];
    %imagesc(A);
    %hold off;
    %hold on; contour(F(:,:,50)>0,[0,0],'r'); hold off;
    %subplot(1,3,2); imagesc(F(:,:,50)); title('Level Set'); 
    %subplot(1,3,3); imagesc(F(:,:,99)); title('Level Set'); 
    
    %subplot(122); imagesc(f); title(['CV \lambda=',num2str(lambda)]);  xlabel(['t=',num2str(t)]);
    %hold on; contour(F(),[0,0],'r'); hold off;
    %colormap gray;
   
    hold off;
    %hold on; plot3(F==0)
    %drawnow;
   
end;
figure;
h = vol3d('cdata',F,'texture','3D'); 
view(3);
daspect([1 1 1]);
alphamap('default');
alphamap(.1 .* alphamap);
