% a script for controling multiple loops under different parameters
% by Zhiyuan Ding

% loop list
% LoopList = [0,113,565,2825];
LoopList = {''};
% LoopList = {'_elastic'};

% loop
for ii = 1:length(LoopList)
    Identifier = ['dps_grc_round3',LoopList{ii},'_code_test01'];
    DataPath = ['../mustem/dps_grc_padding_4dstem_round4_CA_30mrad',LoopList{ii},'.muh5']; % datapath, hdf5 or muh5 is recommended
    SavePath = sprintf('../results/%s/',Identifier);
    
    fprintf('%s[LOG]: loop ctrl --> %d/%d\n\t\t\tIdentifier --> %s\n',datestr(now,31),ii,length(LoopList),Identifier);
    
    % check SavePath
    if ~exist(SavePath,'dir')
        mkdir(SavePath);
    end
    
    % get ptycho parameters and load dataset
    ptycho = generate_ptycho_para(Identifier, DataPath, SavePath);
    % change ptycho parameters here:
    % ptycho.xx = xx;
%     ptycho.mask_mode = 'none';
    
    % save ptycho para
    save_ptycho_para(ptycho, SavePath);
    
    % get ePIEpara
%     ePIEpara = generate_ePIEpara(ptycho, [SavePath]);
    % change ePIE parameters here:
    % ePIEpara.xx = xx;
    % ePIEpara.update_function.algorithm = 'PIE';
    % ePIEpara.recon.guess_init_probe = 1; % 1==initialise probe, 2==load init_probe
    % ePIEpara.statics.init_probe = '../results/init_probe.mat'; % name of var: "probe"
%     ePIEpara.switch.iteration_mask_mode = 'none';
    
    % save ePIEpara
%     save_ePIEpara(ePIEpara, [SavePath]);
    
    % calculation
    cal_ssb_wdd(ptycho, SavePath);
%     cal_ePIE(ePIEpara, [SavePath]);
    
    % post-process
    %

end

%%


function ptycho = generate_ptycho_para(Identifier, DataPath, SavePath)
% init ptycho (parameters)
ptycho = a_loop_control_parameters();

% parameters change
ptycho.identifier = Identifier;
% Acceleration voltage in [kV]
ptycho.voltage_kV  = 60;
% Wavelength in [Angtrom]
ptycho.wavelength = Wavelength(ptycho.voltage_kV); 
% Probe convergence angle in [Radian]
ptycho.ObjApt_angle = 0.03; 
% Rotation angle in [Radian]
ptycho.rot_angle = 0;  
% Image pixel size in [Angstrom]. This corresponds to the sampling of the
% probe positions and it should be defined as [pixel_size_in_x pixel_size_in_y];
ptycho.pix_size = [0.1,0.1]; % in A, zyd added

ptycho.aberr.C1   = 0  ; % Defocus      Unit: meters. (common value = 1nm)
ptycho.aberr.C12a = 0  ; % A1, 2-fold Stig, a direction Unit: meters. (Common value = 5nm)
ptycho.aberr.C12b = 0  ; % A1, 2-fold Stig, b direction Unit: meters. (Common value = 5nm)
ptycho.aberr.C21a = 0  ; % B2, Axial Coma   Unit: meters. (Common value = 250nm)
ptycho.aberr.C21b = 0  ; % B2, Axial Coma   Unit: meters. (Common value = 250nm)
ptycho.aberr.C23a = 0  ; % A2, 3-fold Stig  Unit: meters. (Common value = 250nm)
ptycho.aberr.C23b = 0  ; % A2, 3-fold Stig  Unit: meters. (Common value = 250nm)
ptycho.aberr.C30a = 0  ; % C3, Cs           Unit: meters. (Common value = 50um)
ptycho.aberr.C32a = 0  ; % S3, 3rd-order 2-fold Stig  Unit: meters. (Common value = 20um)
ptycho.aberr.C32b = 0  ; % S3, 3rd-order 2-fold Stig  Unit: meters. (Common value = 20um)
ptycho.aberr.C34a = 0  ; % A3, 4-fold Stig  Unit: meters. (Common value = 25um)
ptycho.aberr.C34b = 0  ; % A3, 4-fold Stig  Unit: meters. (Common value = 25um)
ptycho.aberr.C50a = 0  ; % C5, 5th ord. Cs  Unit: meters. (Common value = 50mm)
ptycho.aberr.C56a = 0  ; % A5, 6-fold Stig  Unit: meters. (Common value = 50mm)

% other parameter change

% Read dataset
fprintf('%s[LOG]: reading 4d dataset\n',datestr(now,31));
ptycho = a_loop_control_read4ddataset(DataPath,ptycho);
end

function ePIEpara = generate_ePIEpara(ptycho, SavePath)
% ePIE parameter
ePIEpara = ePIE_parameters();
ePIEpara.dp.dk = 1.9; % mrad
ePIEpara.iter = 5; % iteration number
ePIEpara.recon.dp_order = 1; % 0=by index.1=random order.
ePIEpara.update_function.algorithm = 'ePIE'; % 'ePIE'/'PIE'
ePIEpara.update_function.obj_alpha = 0.1;
ePIEpara.update_function.probe_alpha = 0.1;
ePIEpara.switch.iteration_mask_mode = 'semiangle';% 'none'/'semiangle'/'mask-file'.
% 'none' --> no mask, 'semiangle' --> generate a mask by semiangle.
% 'mask-file' --> load p.static.iteration_mask_filename, use var "mask" in the mat file.
ePIEpara.switch.gpu = false;
ePIEpara.switch.precision = 'single';
ePIEpara.coeffi(1,:) = [0 0];%nm      C10     df
ePIEpara.coeffi(2,:) = [0 0];%    nm      C12     A1
ePIEpara.coeffi(3,:) = [0 0];%    nm      C21     B2
ePIEpara.coeffi(4,:) = [0 0];%    nm      C23     A2
ePIEpara.coeffi(5,:) = [0 0];%    um      C3      C3
ePIEpara.coeffi(6,:) = [0 0];%    um      C32     S3
ePIEpara.coeffi(7,:) = [0 0];%    um      C34     A3
ePIEpara.coeffi(8,:) = [0 0];%    um      C41     A4
ePIEpara.coeffi(9,:) = [0 0];%    um      C43     B4
ePIEpara.coeffi(10,:) = [0 0];%   um      C45     D4
ePIEpara.coeffi(11,:) = [0 0];%   mm      C5      C5
ePIEpara.coeffi(12,:) = [0 0];%   mm      C56     A5

% other ePIE para change

% load parameters from SSB setting
ePIEpara.semiangle = ptycho.ObjApt_angle * 1000;
ePIEpara.beam_voltage = ptycho.voltage_kV;

% ePIE parameters from dataset
ePIEpara.dp.num_all = size(ptycho.m,3) * size(ptycho.m,4);
if size(ptycho.m,1) == size(ptycho.m,2)
    ePIEpara.dp.size_raw = size(ptycho.m,1);
else
    error('The ePIE here only support square diffraction patterns');
end
% about trans
ePIEpara.trans.scan_num = [size(ptycho.m,3),size(ptycho.m,4)];
ePIEpara.trans.scan_step = ptycho.pix_size; % scan step in A
ePIEpara.trans.theta = [ptycho.rot_angle,ptycho.rot_angle]/180*pi; % [theta1, theta2]/180*pi, theta1&2 in degree

% dps, dps are not imported here, because the parameters should
% be saved with out dataset.
ePIEpara.datapath = ptycho.m;
end

function cal_ssb_wdd(ptycho,SavePath)
% Calculate SSB
fprintf('%s[LOG]: calculating SSB\n',datestr(now,31));
[ptycho] = define_real_space(ptycho);
[ptycho] = detector_binning(ptycho);
[ptycho] = define_reciprocal_space(ptycho);
[ptycho] = calculate_detector_masks(ptycho);
[ptycho] = truncate_ronchigrams_within_collection_angle(ptycho);
[ptycho] = define_reciprocal_space_for_truncated_detector_plane(ptycho);
[ptycho] = FT_from_M_to_G(ptycho);
[ptycho] = single_side_band_reconstruction(ptycho);
Res.ssb_r = ptycho.trotterimgR;
Res.ssb_l = ptycho.trotterimgL;
fprintf('%s[LOG]: SSB finished\n',datestr(now,31));

% Calculate WDD
fprintf('%s[LOG]: calculation WDD\n',datestr(now,31));
[ptycho] = define_reciprocal_space_resampled_ronchigram(ptycho);
[ptycho] = define_probe_function_Modified(ptycho);
[ptycho] = FT_from_G_to_H(ptycho);
[ptycho] = wigner_distribution_deconvolution_reconstruction_Modified(ptycho);
Res.wdd = ptycho.psi;
Res.wdd_probe = ptycho.func_probe;
fprintf('%s[LOG]: WDD finished\n',datestr(now,31));

% save result
save(sprintf('%sRes.mat',SavePath),'Res');
imwrite(mat2gray(abs(Res.ssb_r)),sprintf('%sssb_amp.png',SavePath));
imwrite(mat2gray(angle(Res.ssb_r)),sprintf('%sssb_phs.png',SavePath));
imwrite(mat2gray(abs(Res.wdd)),sprintf('%swdd_amp.png',SavePath));
imwrite(mat2gray(angle(Res.wdd)),sprintf('%swdd_phs.png',SavePath));
imwrite(mat2gray(abs(Res.wdd_probe)),sprintf('%swdd_probe_amp.png',SavePath));
imwrite(mat2gray(angle(Res.wdd_probe)),sprintf('%swdd_probe_phs.png',SavePath));
fprintf('%s[LOG]: ssb and wdd results saved\n',datestr(now,31));

end

function cal_ePIE(ePIEpara, SavePath)
% Calculate ePIE
fprintf('%s[LOG]: calculation ePIE\n',datestr(now,31));
[probe,obj,dps,trans_exec,ePIEpara] = ePIE_preprocess(ePIEpara);
ePIEres = ePIE_main(probe,obj,dps,trans_exec,ePIEpara);

% save ePIE result
save(sprintf('%sePIEres.mat',SavePath),'ePIEres');
ePIE_show_wf(ePIEres.obj,'Save',true,...
    'SaveNameStartWith',sprintf('%sePIE',SavePath),...
    'Trans',ePIEres.trans_exec,'DrawFigure',false);
ePIE_show_wf(ePIEres.probe,'Save',true,...
    'SaveNameStartWith',sprintf('%sePIE_probe',SavePath),...
    'DrawFigure',false);
fprintf('%s[LOG]: ePIE results saved\n',datestr(now,31));

end

function save_ptycho_para(ptycho, SavePath)
ptycho.m = [];
save(sprintf('%spara_raw.mat',SavePath),'ptycho');
end

function save_ePIEpara(ePIEpara, SavePath)
ePIEpara.datapath = [];
save(sprintf('%sePIEpara_raw.mat',SavePath),'ePIEpara');
end

