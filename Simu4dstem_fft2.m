function xf = Simu4dstem_fft2(x)
%EPIE_FFT2
xf = fftshift(fft2(ifftshift(x)))/sqrt(numel(x(:)));
end

