function [F] = cv_3d(f, lambda)

% Chan-Vese Active Contours Segmentation using Level Sets
% Input: Original grayscale image f
%        Parameter lambda (Try lambda=0.1)
% Output: Level Set Function Phi

%Set default value for lambda if user does not provide lambda.
if ~exist('lambda')
    lambda = 0.1;
end;

%Parameters
T = 100;    %Stopping time
dt = 0.2;  %Time step
a = 0.01;  %Fudge factor to avoid division by zero.
epsilon = 0.1;  %Epsilon in delta approximation
%Fidelity terms lambda_in & lambda_out could be different.
%Both set = lambda here for simplicity.
lambda_in = lambda;
lambda_out = lambda;  

I = im2double(f);
F = im2double(f);
[m,n,p,k] = size(f);

%Initialize with large box around perimeter of image.
%F(1:m,1:n,1:p) = -1;  F(2:m-1,2:n-1,2:p-1)=1;
F(1,:,:) = 0;
F(:,1,:) = 0;
F(:,:,1) = 0;
F(m,:,:) = 0;
F(:,n,:) = 0;
F(:,:,p) = 0;

for t = 0:dt:T
    %Approximate derivatives
    F_x = (F([2:m,m],:,:) - F(1:m,:,:));
    F_y = (F(:,[2:n,n],:) - F(:,1:n,:));
    F_z = (F(:,:,[2:p,p]) - F(:,:,1:p));
    
    F_xx = (F_x([2:m,m],:,:) - F_x(1:m,:,:));
    F_yy = (F_y(:,[2:n,n],:) - F_y(:,1:n,:));
    F_zz = (F_z(:,:,[2:p,p]) - F_z(:,:,1:p));
    
    F_xy = (F_x(:,[2:n,n],:) - F_x(:,1:n,:)); 
    F_yz = (F_y(:,:,[2:p,p]) - F_y(:,:,1:p));
    F_xz = (F_x(:,:,[2:p,p]) - F_x(:,:,1:p));
    
    %TV term = Num/Den
    Num = (F_z.*(F_xx.*F_z- 2.*F_x.*F_xz) +(F_x.^2).* F_zz).*(F_z.*(F_yy.*F_z- 2.*F_y.*F_yz) +(F_y.^2) .* F_zz)-((F_z.*(-1.*F_x.*F_yz- F_xy.*F_z - F_xz.*F_y) +F_x.* F_y.* F_zz).^2);
    Den = F_z.^2.*(F_x.^2+F_y.^2+F_z.^2);
    
    c_in = sum(sum(sum([F>0].*I)))./(a+sum(sum(sum([F>0]))));
    c_out = sum(sum(sum([F<0].*I)))./(a+sum(sum(sum([F<0]))));
    
    %Add to previous iteration of u.
    %size(F)
    %size(dt*epsilon./(pi*(epsilon^2+F.^2)))
    %size(Num ./ Den )
    %size(lambda_in)
    %fprintf('f :\n');
    %size(f)
    %size(c_in)
    %size(lambda_in*(f-c_in).^2)
    %size(lambda_out*(f-c_out).^2)
    
    F = F + dt*epsilon./(pi*(epsilon^2+F.^2)).*( Num./Den  - lambda_in*(I-c_in).^2 + lambda_out*(I-c_out).^2);
    


    
    %Plot results.  Note drawing every iteration slows down the program.
    %subplot(121); imagesc(Phi); title('Level Set'); 
    %subplot(122); imagesc(f); title(['CV \lambda=',num2str(lambda)]);  xlabel(['t=',num2str(t)]);
    %hold on; contour(Phi,[0,0],'r'); hold off;
    %colormap gray;
    h = vol3d('cdata',f,'texture','3D'); 
    view(3);
    daspect([1 1 1]);
    alphamap('rampup');
    alphamap(.06 .* alphamap);
    % hold off;
    %hold on; plot3(F==0)
    %drawnow;
    t
end;
