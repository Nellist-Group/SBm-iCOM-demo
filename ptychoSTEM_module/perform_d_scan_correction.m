function [ptycho] = perform_d_scan_correction(ptycho)

%%%%%% DscanCoor_v1 not available, check with Hao


% Centre of Mass Map and D-scan correction

% input:  M, input_params.ObjApt_angle
% output:
% pacbed - average of all detector frames
% s -  structural variable that contains the centre and radius of the BF disc,
%      In case of binned detectors, it outputs the centre, major and
%      minor axis length of the ellipsis.
% Detector geometries - maskBF, maskDF, maskABF, maskIBF, maskWP, maskPIE
% Synthetic Images - BFimg, IBFimg, ABFimg, DFimg
% theta, theta_x, theta_y - Detector Plane Angle Calibration based on given ObjApt_angle

% run pacbedAnalyzer_v1.m;

% Quantifying the Centre of Mass of each Ronchigram.
run RonchiCOM_v1.m;

% based on the disc shift measured above, correct the D-scan
if IfCorrectDscan
    run DscanCorr_v1.m;
    % run PACBED again after D-scan correction;
    % check improvements in the synthetic ADF, BF images
    % check improvements in the PACBED
    run pacbedAnalyzer_v1.m;
    % run Ronchigram Centre of Mass again
    run RonchiCOM_v1.m;
    % based on the disc shift measured above, correct the D-scan
    run DscanCorr_v1.m;
end



end