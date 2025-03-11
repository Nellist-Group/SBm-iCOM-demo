function [ptycho] = calculate_G_power_spectrum_wrt_Kf(ptycho)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To calculate the power spectrum of G with respect to Kf
% Function requires structured array ptycho.
% Input variables required in ptycho are:
% ptycho.varfunctions.FT_from_M_to_G: Flag to indicate the function has
% been executed
% ptycho.ObjSize: vector with the pixel of the scanned area + zero-padding (if
% stated in parameter file). Padding is useful to avoid artifacts when
% doing FFTs.
% ptycho.G_wp: 4D array that contains the FFT of of matrix M with respect
% to probe position
% Output variables added to ptycho are:
% ptycho.psG: matrix containing the power spectrum of G with respect to Kf
% ptycho.varfunctions.calculate_G_power_spectrum_wrt_Kf: Flag to indicate that 
% the function has been executed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if ptycho.varfunctions.FT_from_M_to_G
        psG = zeros(ptycho.ObjSize);
        for yy=1:ptycho.ObjSize(1)   
            for xx=1:ptycho.ObjSize(2) 
                g =  squeeze(ptycho.G_wp(:,:,yy,xx));
                psG(yy,xx) = sum(sum(g.*conj(g)));
            end
        end
 else
     display('ptycho.G_wp variable is not in ptycho array');
     display('calculate it using FT_from_M_to_G.m function');
 end
    
 ptycho.psG = psG;
 ptycho.varfunctions.calculate_G_power_spectrum_wrt_Kf = 1;
end