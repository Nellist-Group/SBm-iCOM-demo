function imgOut = Simu4dstem_spot_size_conv(img, dx, spotSize)
%SIMU4DSTEM_SPOT_SIZE_CONV a gaussian spotsize function to convolute with
%image
%   input:
%       img --> the image to convolute
%       dx --> dx of the image, in A
%       spotSize --> spotSize of beam source in A, as well as the FWHM of
%       the gaussion convolution kernal

% imgOut = imgaussfilt(img, spotSize/2.35482/dx);
kernelSize = ceil(spotSize/2.35482/dx) .* 2 + 1;

k = fspecial('gaussian', kernelSize, spotSize/2.35482/dx);
imgp = padarray(img,[(kernelSize-1)/2,(kernelSize-1)/2],'symmetric');

imgOut = conv2(imgp,k,'valid');

end

