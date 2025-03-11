function xf = Simu4dstem_fft2_4d34(x)
%SIMU4DSTEM_FFT2_4D12 Summary of this function goes here
%   Detailed explanation goes here

xf = fftshift(fftn(ifftshift(x)))/sqrt(numel(x(:)));
xf = fftshift(ifft2(ifftshift(xf))).*sqrt(size(xf,1).* size(xf,2));

end

