function [mnoise] = add_noise_to_CBEDs(m,dose,pixel_size)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function adds Poisson distributed noise to the simulated CBED
% patterns.
% Input variables required are:
% m: 4D array that contains the ronchigram for each probe position. This
% variable is usually stored in ptycho.m
% dose: integer value indicating the number of electrons per Angstrom
% squared illuminating the sample. This variable is defined in the parameter
% file as ptycho.dose
% pixel_size: Image pixel size in [Angstrom]. This corresponds to the sampling of the
% probe positions and it should be defined as [pixel_size_in_x
% pixel_size_in_y]. This variable is defined in the parameter file as
% ptycho.pix_size
% Output variables are:
% mnoise: 4D array that contains the modified ronchigram for each probe
% position. Ronchigrams will be affected by the noise corresponding to the
% indicated dose. 
%
% Intended use of this function:
% Adding noise to the CBED patterns is intented for more realistic
% simulations analyses. The intended use of this function under the
% 'ptycho' framework is:
% ptycho.m = add_noise_to_CBEDs(ptycho.m,ptycho.dose,ptycho.pix_size);
% IMPORTANT: Note that the variable 'ptycho.m' will be modified and the original
% simulated CBED patterns are lost. If you would like to keep the original
% CBED patterns, then ptycho.m should be stored in a different variable,
% since the processing functions will work on the variable ptycho.m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mnoise = zeros(size(m));

doseperpixel = dose*(pixel_size(1)*pixel_size(2));

for i = 1:1:size(m,3)
    for j =1:1:size(m,4)
       i
      j 

        %mnoise(:,:,i,j) =  poissrnd(doseperpixel*m(:,:,i,j)/sum(sum(m(:,:,i,j))));
        mnoise(:,:,i,j) =  poissrnd(doseperpixel*m(:,:,i,j));

    end
end


end
