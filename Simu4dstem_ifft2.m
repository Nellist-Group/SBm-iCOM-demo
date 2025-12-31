function xif = Simu4dstem_ifft2(x)
%EPIE_IFFT2 

xif = fftshift(ifft2(ifftshift(x)))*sqrt(numel(x(:)));

end

