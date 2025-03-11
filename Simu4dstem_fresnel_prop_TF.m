function [OW] = Simu4dstem_fresnel_prop_TF(IW,lambda,z,dxy)
%FRESNELPROPTF transfer function approach.
% input:
%   IW --> incident wave
%   lambda --> wave length in m
%   z --> propogation distance in m
%   dxy --> sampling, the size of each pixel in m
%
% Zhiyuan Ding
% Ref: Computational Fourier Optics, a MATLAB tutorial, P67

[M,N] = size(IW);
% k = 2*pi/lambda;
fx = -1/(2*dxy(1)):1/(M*dxy(1)):1/(2*dxy(1))-1/(M*dxy(1));
fy = -1/(2*dxy(2)):1/(N*dxy(2)):1/(2*dxy(2))-1/(N*dxy(2));
[FX,FY] = meshgrid(fy,fx);
H = exp(-1j*pi*lambda*z*(FX.^2 + FY.^2));
H = fftshift(H);
IW = fft2(fftshift(IW));
OW = H .* IW;
OW = ifftshift(ifft2(OW));


end

