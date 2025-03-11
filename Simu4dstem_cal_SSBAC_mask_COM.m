function [ptycho_0,Res] = Simu4dstem_cal_SSBAC_mask_COM(jsonData,dps_4d)
%SIMU4DSTEM_CAL_SSB_AND_WDD Summary of this function goes here
%   Detailed explanation goes here
addpath('ptychoSTEM_module/');
ptycho = a_loop_control_parameters();
% parameters change
ptycho.identifier = 'Simu4dstem_cal_SSBAC_mask_COM';
% Acceleration voltage in [kV]
ptycho.voltage_kV  = jsonData.voltage;
% Wavelength in [Angtrom]
ptycho.wavelength = Wavelength(ptycho.voltage_kV);
% Probe convergence angle in [Radian]
ptycho.ObjApt_angle = jsonData.alpha * 1e-3;
% Rotation angle in [Radian]
ptycho.rot_angle = 0;
% Image pixel size in [Angstrom]. This corresponds to the sampling of the
% probe positions and it should be defined as [pixel_size_in_x pixel_size_in_y];
ptycho.pix_size = [jsonData.step_size,jsonData.step_size]; % in A, zyd added

ptycho.aberr.C1 = -jsonData.defocus * 1e-10  ; % Defocus      Unit: meters. (common value = 1nm)
ptycho.aberr.C12a = -jsonData.cond_lens_c_12_amp_A * 1e-10  ; % A1, 2-fold Stig, a direction Unit: meters. (Common value = 5nm)
ptycho.aberr.C12b = jsonData.cond_lens_c_12_angle ; % A1, 2-fold Stig, b direction Unit: meters. (Common value = 5nm)
ptycho.aberr.C21a = 0  ; % B2, Axial Coma   Unit: meters. (Common value = 250nm)
ptycho.aberr.C21b = 0  ; % B2, Axial Coma   Unit: meters. (Common value = 250nm)
ptycho.aberr.C23a = jsonData.cond_lens_c_23_amp_A * 1e-10  ; % A2, 3-fold Stig  Unit: meters. (Common value = 250nm)
if jsonData.cond_lens_c_23_amp_A == 0
    ptycho.aberr.C23b = 0;
else
    ptycho.aberr.C23b = jsonData.cond_lens_c_23_angle; % A2, 3-fold Stig  Unit: meters. (Common value = 250nm)
end
ptycho.aberr.C30a = -jsonData.cond_lens_c_30_mm * 1e-3  ; % C3, Cs           Unit: meters. (Common value = 50um)
ptycho.aberr.C32a = 0  ; % S3, 3rd-order 2-fold Stig  Unit: meters. (Common value = 20um)
ptycho.aberr.C32b = 0  ; % S3, 3rd-order 2-fold Stig  Unit: meters. (Common value = 20um)
ptycho.aberr.C34a = 0  ; % A3, 4-fold Stig  Unit: meters. (Common value = 25um)
ptycho.aberr.C34b = 0  ; % A3, 4-fold Stig  Unit: meters. (Common value = 25um)
ptycho.aberr.C50a = 0  ; % C5, 5th ord. Cs  Unit: meters. (Common value = 50mm)
ptycho.aberr.C56a = 0  ; % A5, 6-fold Stig  Unit: meters. (Common value = 50mm)

% other parameter change

% get ptycho_0
ptycho_0 = ptycho;

% Read dataset
ptycho.m = dps_4d;
ptycho.probe_range = [0 0 size(ptycho.m,4) size(ptycho.m,3)];
ptycho.varfunctions.load_data = 1;

% Calculate WDD
fprintf('%s[LOG]: calculation SSBAC mask COM\n',datestr(now,31));
[ptycho] = define_real_space(ptycho);
[ptycho] = detector_binning(ptycho);
[ptycho] = define_reciprocal_space(ptycho);
[ptycho] = calculate_detector_masks(ptycho);
[ptycho] = truncate_ronchigrams_within_collection_angle(ptycho);
[ptycho] = define_reciprocal_space_for_truncated_detector_plane(ptycho);
[ptycho] = FT_from_M_to_G(ptycho);
[ptycho] = define_reciprocal_space_resampled_ronchigram(ptycho);
[ptycho] = define_probe_function_Modified2(ptycho); %define_probe_function_Modified2 use JEOL and CEOS notation.
[ptycho] = FT_from_G_to_H(ptycho);
% [ptycho] = wigner_distribution_deconvolution_reconstruction_Modified2(ptycho);
[ptycho]= Aberration_corrected_single_side_band_reconstruction_COM(ptycho);
% Res.ssb_r = ptycho.trotterimgRC;
% Res.ssb_l = ptycho.trotterimgLC;
Res.iCOM = ptycho.iCOM;
Res.dCOM = ptycho.dCOM;
Res.COMx = ptycho.COMx;
Res.COMy = ptycho.COMy;

% spot size convolution
% if jsonData.spot_size_A > 0
%    Res.ssb_r = Simu4dstem_spot_size_conv(Res.ssb_r,jsonData.step_size,jsonData.spot_size_A);
%    Res.ssb_l = Simu4dstem_spot_size_conv(Res.ssb_l,jsonData.step_size,jsonData.spot_size_A);
% end
% 
% if isfield(ptycho,'trotterimgL_ss')
%     Res.ssb_l_ss = ptycho.trotterimgL_ss;
% end
% if isfield(ptycho,'trotterimgR_ss')
%     Res.ssb_r_ss = ptycho.trotterimgR_ss;
% end
% if isfield(ptycho,'func_probe')
%     Res.probe = ptycho.func_probe;
% end

fprintf('%s[LOG]: SSBAC mask COM finished\n',datestr(now,31));

% save result
% save(sprintf('%sRes.mat',ptycho.savepath),'Res');
% imwrite(mat2gray(abs(Res.ssb_r)),sprintf('%sssb_amp.png',ptycho.savepath));
% imwrite(mat2gray(angle(Res.ssb_r)),sprintf('%sssb_phs.png',ptycho.savepath));
% imwrite(mat2gray(abs(Res.wdd)),sprintf('%swdd_amp.png',ptycho.savepath));
% imwrite(mat2gray(angle(Res.wdd)),sprintf('%swdd_phs.png',ptycho.savepath));
% imwrite(mat2gray(abs(Res.wdd_probe)),sprintf('%swdd_probe_amp.png',ptycho.savepath));
% imwrite(mat2gray(angle(Res.wdd_probe)),sprintf('%swdd_probe_phs.png',ptycho.savepath));
% fprintf('%s[LOG]: ssb and wdd results saved\n',datestr(now,31));
end

