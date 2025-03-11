function [ptycho] = FT_from_M_to_G(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To Fourier transform M with respect to probe position in order to build G.
% The purpose of this step is found in J.M. Rodenburg, et al. 
% Ultramicroscopy 48 (1993) 304-314.
% Function requires structured array ptycho. 
% Input variables required are:
% ptycho.m_wp: 4 D m-matrix with truncated ronchigrams
% ptycho.ObjSize: vector with the pixel of the scanned area + zero-padding
% Output variables added to ptycho are:
% ptycho.mean_m_wp: matrix containing the pacbed of the truncated
% ronchigram
% ptycho.sumI: sum of the total intensity in all ronchigrams and probe positions
% ptycho.G_wp: 4D array that contains the ronchigram for each spatial frecuency of probe position   
% ptycho.size4D: vector containing the dimensions of ptycho.G_wp
% ptycho.varfunctions.FT_from_M_to_G: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size4D = [size(ptycho.m_wp,1), size(ptycho.m_wp,2), ptycho.ObjSize(1), ptycho.ObjSize(2)];

display('FT{M} w.r.t Probe Positions...')
% allocate memory for G matrix
G_wp = zeros(size4D,'like',ptycho.m_wp);
%Now take FFT of each detector pixel with respect to probe position
%ImagePadMean pads the average intensity of the detector over the whole
%scan around the probe possitions actually scanned if a pad is set. It does
%not check if ptycho.use_padding is set so if no pad is set, it just 
%catches errors in your data loading..
tmp = zeros(1,1,ptycho.ObjSize(1),ptycho.ObjSize(2),'like',G_wp);
for i = 1:size4D(1)
    for j=1:size4D(2)
        %disp(['FT{M}:',num2str(i), '/', num2str(j)]);
        a = ImagePadMean(squeeze(ptycho.m_wp(i,j,:,:)),ptycho.ObjSize);
        b = fftshift(fft2(ifftshift(a)));
        tmp(1,1,:,:) =  b;
        G_wp(i,j,:,:) = tmp;
%         G_wp(i,j,:,:) = fftshift(fft2(ifftshift(ImagePadMean(squeeze(ptycho.m_wp(i,j,:,:)),ptycho.ObjSize))));
    end
end

% calculate the average CBED from all probe positions
if ~exist('mean_m_wp','var')
    mean_m_wp = mean(mean(ptycho.m_wp,4),3);
end
if ~exist('sumI','var')
    sumI = sum(ptycho.m_wp(:));
end

ptycho.size4D = size4D;
ptycho.mean_m_wp = mean_m_wp;
ptycho.sumI = sumI;
ptycho.G_wp = G_wp;
ptycho.varfunctions.FT_from_M_to_G = 1;

end