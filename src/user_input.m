function curve_evolution(I)

    f = figure;
    imshow(I);
    %h = imfreehand(gca);
    %h = imellipse(gca);
    h = imrect(gca);
    position = wait(h);
    disp(position);
    mask = h.createMask();
    imshow(mask);
    title('Binary mask of the region', 'FontSize', 10);
end
