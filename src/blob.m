function blob(issolid,T)

bx = 100;
by = 100;
bz = 100;


cx = 50;
cy = 50;
cz = 50;

r = 25;


[x y z] = meshgrid(1:bx,1:by,1:bz);
%X = [x(:), y(:), z(:)]';
mu = [cx, cy,cz]';
%sigma = ones(3);
sigma = 0.25;
n = 3;



Bsolid = sqrt((x-cx).^2 + (y-cy).^2 + (z-cz).^2) <= r;
Bblob = (1.0/sqrt(2*pi)) * exp(- 0.5 * sigma * ( sqrt((x-cx).^2 + (y-cy).^2 + (z-cz).^2)));

if(strcmp(issolid,'solid'))
    B = Bsolid;
else
    B = Bblob;
end

Bnorm = (B - min(min(min(B)))) ./(max(max(max(B))) - min(min(min(B))));

%BVol = vol3d('cdata', 255*Bnorm, 'texture', '3D');
%view(3);

cv_3d(255*Bnorm,'no',T);

end