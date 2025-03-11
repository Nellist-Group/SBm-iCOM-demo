function [ptycho] = parallel_FT_from_M_to_G(ptycho)
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

size4D = [size(ptycho.m_wp,1),size(ptycho.m_wp,2), ptycho.ObjSize(1), ptycho.ObjSize(2)];
size3D = [size(ptycho.m_wp,1)*size(ptycho.m_wp,2), size(ptycho.m_wp,3),size(ptycho.m_wp,4)];
size3Dpadd = [size(ptycho.m_wp,1)*size(ptycho.m_wp,2), ptycho.ObjSize(1), ptycho.ObjSize(2)];

mres = reshape(ptycho.m_wp,size3D);
ObjSize = ptycho.ObjSize;

display('FT{M} w.r.t Probe Positions...')
% allocate memory for G matrix
G_wp = zeros(size3Dpadd);

parfor i = 1:size3Dpadd(1)
        %disp(['FT{M}:',num2str(i)]);
        G_wp(i,:,:) = fftshift(fft2(ifftshift( ImagePadMean(squeeze(mres(i,:,:)),ObjSize))));
end

G_wp = reshape(G_wp,size4D);

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