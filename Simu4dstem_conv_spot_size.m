function img_conv = Simu4dstem_conv_spot_size(img, spot_size_pxl)
%SIMU4DSTEM_CONV_SPOT_SIZE convolute a spot size on an image
%   

hsize = ceil(spot_size_pxl * 4);
if rem(hsize,2) == 0
    hsize = hsize + 1;
end
h = fspecial('gaussian',hsize,spot_size_pxl);
img_conv = imfilter(img,h,'same','conv');

end

