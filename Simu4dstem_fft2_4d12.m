function xf = Simu4dstem_fft2_4d12(x)
%SIMU4DSTEM_FFT2_4D12 Summary of this function goes here
%   Detailed explanation goes here

xf = fftshift(fft2(ifftshift(x)))/sqrt(size(x,1).* size(x,2));
end

