function xf = Simu4dstem_fft4(x)
%SIMU4DSTEM_FFT4 Summary of this function goes here
%   Detailed explanation goes here

xf = fftshift(fftn(ifftshift(x)))/sqrt(numel(x(:)));

end

