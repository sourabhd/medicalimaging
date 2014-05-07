function [I] = read_image(imgfile)

    dbstop if error;
    J = imread(imgfile);
    if size(J,3) > 1
        J = rgb2gray(J);
    end
    I = im2double(J);
end