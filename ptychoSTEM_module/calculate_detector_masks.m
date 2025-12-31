function [ptycho] = calculate_detector_masks(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To calculate and define the different detector masks applied to the
% ronchigram in order to make synthetic images, such as BF, IBF, ABF, ADF
% and DPC
% Function requires structured array ptycho.
% Input variables required are:
% ptycho.ObjApt_angle: Probe convergence angle
% ptycho.rot_angle: Rotation angle from trotters?  
% ptycho.pacbed: matrix containing the Position Averaged Convergent Beam
% Electron Diffraction (PACBED) pattern.
% ptycho.tx: matrix of scattering angle in x from meshgrid
% ptycho.ty: matrix of scattering angle in y from meshgrid
% ptycho.theta: matrix of scattering angles
% Output variables added to ptycho are:
% ptycho.maskBF: matrix of mask for Bright Field detector
% ptycho.maskABF: matrix of mask for Annular Bright Field detector
% ptycho.maskIBF: matrix of mask for Incoherent Bright Field detector
% ptycho.maskDF: matrix of mask for Annular Dark Field detector
% ptycho.maskWP: matrix of mask for Weak Phase Ptychography
% ptycho.maskPIE: matrix of mask for PIE
% ptycho.R11: matrix of mask for DPC detector
% ptycho.R12: matrix of mask for DPC detector
% ptycho.R13: matrix of mask for DPC detector
% ptycho.R14: matrix of mask for DPC detector
% ptycho.R21: matrix of mask for DPC detector
% ptycho.R22: matrix of mask for DPC detector
% ptycho.R23: matrix of mask for DPC detector
% ptycho.R24: matrix of mask for DPC detector
% ptycho.R31: matrix of mask for DPC detector
% ptycho.R32: matrix of mask for DPC detector
% ptycho.R33: matrix of mask for DPC detector
% ptycho.R34: matrix of mask for DPC detector
% ptycho.varfunctions.calculate_detector_masks: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a mask for BF disc, DF ring, ABF ring.
mask0 = ones(size(ptycho.theta));
% mask for weak phase ptychography using only BF discs
% for thick samples, there's significant intensities outside Bf,
% therefore
maskWP = mask0; 
if isfield(ptycho,'mask_mode') && strcmp(ptycho.mask_mode,'none')
    % do not change maskWP
else % default: no field ptycho.mask_mode, or mask_mode = 'semiangle'
    maskWP(ptycho.theta>ptycho.ObjApt_angle*1.05)=0;
end
%maskWP(ptycho.theta>ptycho.ObjApt_angle*1.3)=0;
% mask for PIE using both BF and DF signals.
maskPIE = mask0; 
maskPIE(ptycho.theta>ptycho.ObjApt_angle*1.3)=0;

% IBF
maskIBF = mask0; 
maskIBF(ptycho.theta>ptycho.ObjApt_angle*1.00)=0;
% BF
maskBF  = mask0; 
maskBF(ptycho.theta>ptycho.ObjApt_angle*0.2)=0;
% ABF
maskABF = mask0; 
maskABF(ptycho.theta>ptycho.ObjApt_angle*1.0)=0;
maskABF(ptycho.theta<ptycho.ObjApt_angle*0.5)=0;
% ADF
maskDF  = mask0; 
maskDF(ptycho.theta<ptycho.ObjApt_angle*1.05)=0;
%maskDF(ptycho.theta>ptycho.ObjApt_angle*5.52)=0;
% modify the above line(s) if you want to use an outer angle
% or modify inner  angle for the ADF

% DPC
segment_angle = ptycho.rot_angle;

tx_rot = ptycho.tx.*cos(segment_angle) - ptycho.ty.*sin(segment_angle);
ty_rot = ptycho.tx.*sin(segment_angle) + ptycho.ty.*cos(segment_angle);
R1 = mask0; 
R1(ptycho.theta>ptycho.ObjApt_angle*0.5)=0;
R11 = R1; R11(tx_rot<0)=0; R11(ty_rot<0)=0;
R12 = R1; R12(tx_rot<0)=0; R12(ty_rot>0)=0;
R13 = R1; R13(tx_rot>0)=0; R13(ty_rot>0)=0;
R14 = R1; R14(tx_rot>0)=0; R14(ty_rot<0)=0;
R2 = mask0; 
R2(ptycho.theta>ptycho.ObjApt_angle*1.01)=0;
R2(ptycho.theta<ptycho.ObjApt_angle*0.5)=0;
R21 = R2; R21(tx_rot<0)=0; R21(ty_rot<0)=0;
R22 = R2; R22(tx_rot<0)=0; R22(ty_rot>0)=0;
R23 = R2; R23(tx_rot>0)=0; R23(ty_rot>0)=0;
R24 = R2; R24(tx_rot>0)=0; R24(ty_rot<0)=0;
R3 = mask0; 
R3(ptycho.theta<=ptycho.ObjApt_angle*1.02)=0;
R3(ptycho.theta>ptycho.ObjApt_angle*1.25)=0;
R31 = R3; R31(tx_rot<0)=0; R31(ty_rot<0)=0;
R32 = R3; R32(tx_rot<0)=0; R32(ty_rot>0)=0;
R33 = R3; R33(tx_rot>0)=0; R33(ty_rot>0)=0;
R34 = R3; R34(tx_rot>0)=0; R34(ty_rot<0)=0;

if ptycho.plot_figures
figure; imagesc(R11+R12*3+R13*1+R14*2 + R21*3+R22*4+R23*3+R24*4 + R31*5+R32*6+R33*5+R34*6);
axis image
colormap jet

figure;
subplot(2,3,1); imagesc(ptycho.pacbed.*maskBF);axis image;
title('BF detector')
subplot(2,3,2); imagesc(ptycho.pacbed.*maskABF);axis image;
title('ABF detector')
subplot(2,3,3); imagesc(ptycho.pacbed.*maskIBF);axis image;
title('IBF detector')
subplot(2,3,4); imagesc(ptycho.pacbed.*maskDF);axis image;
title('ADF detector')
subplot(2,3,5); imagesc(ptycho.pacbed.*maskWP);axis image;
title('WP')
subplot(2,3,6); imagesc(ptycho.pacbed.*maskPIE);axis image;
title('PIE')
end

ptycho.maskBF = maskBF;
ptycho.maskABF = maskABF;
ptycho.maskIBF = maskIBF;
ptycho.maskDF = maskDF;
ptycho.maskWP = maskWP;
ptycho.maskPIE = maskPIE;
ptycho.R11 = R11;
ptycho.R12 = R12;
ptycho.R13 = R13;
ptycho.R14 = R14;
ptycho.R21 = R21;
ptycho.R22 = R22;
ptycho.R23 = R23;
ptycho.R24 = R24;
ptycho.R31 = R31;
ptycho.R32 = R32;
ptycho.R33 = R33;
ptycho.R34 = R34;
ptycho.varfunctions.calculate_detector_masks = 1;

end

