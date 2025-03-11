% set experimental parameters
dataPath = '../graphene_data/df0/binned_diff_20230819_163033.hdf5'; % datapath for data file
dataFormat = 'h5-frms6'; % 'mat-3d'/'mat-4d'/'h5-frms6'
varName = 'm_31p5mrad'; % name of var in the dataPath, 3d mat
savePath = '../test_df0/';

save_img_path_4_epdp = [savePath,'result_img_epdp','/'];
save_mat_path_4_epdp = [savePath,'result_mat_epdp','/'];
%% parameters
voltage = 80; % kV
convAngle = 24.8; % mrad
dk = 0.6613; % mrad, pixel size on camera
scanNum = [254,255]; % number of scan positions
rotAngle = 129; % 48 degree, the rot angle of dp
scanStep = [0.222,0.222]; % A

STEM_detector_angle_mrad = [60,200;
                            25,50];
DPC_detector_angle_mrad = [10,30];

%% load dataset and set save path
switch dataFormat
    case 'mat-3d'
        load(dataPath,varName);
        eval(['dps = ',varName,';']);
        if ~strcmp('dps',varName)
            clear(varName);
        end
        dps_4d = reshape(dps,[size(dps,1),size(dps,2),scanNum(1),scanNum(2)]);
    case 'mat-4d'
        load(dataPath,varName);
        eval(['dps_4d = ',varName,';']);
        if ~strcmp('dps_4d',varName)
            clear(varName);
        end
        dps_4d = permute(dps_4d,[1,2,3,4]);
        dps = reshape(dps_4d,[size(dps_4d,1),size(dps_4d,2),size(dps_4d,3).*size(dps_4d,4)]);
    case 'h5-frms6'
        dps_4d = h5read(dataPath,'/Experiments/__unnamed__/data');
%         dps_4d = permute(dps_4d,[1,2,4,3]);
        dps_4d = permute(dps_4d,[2,1,4,3]);
        dps = reshape(dps_4d,[size(dps_4d,1),size(dps_4d,2),size(dps_4d,3) * size(dps_4d,4)]);
end

dps_4d = single(dps_4d);
dps = single(dps);

%% subarea
x_range = [127,254];
y_range = [1,128];
dps_4d = dps_4d(:,:,x_range(1):x_range(2),y_range(1):y_range(2));
dps = reshape(dps_4d,[size(dps_4d,1),size(dps_4d,2),size(dps_4d,3) * size(dps_4d,4)]);
scanNum = [x_range(2)-x_range(1)+1,y_range(2)-y_range(1)+1];

%% ADF

if ~exist(save_img_path_4_epdp,'dir')
    mkdir(save_img_path_4_epdp);
end
if ~ exist(save_mat_path_4_epdp,'dir')
    mkdir(save_mat_path_4_epdp);
end

stem_img = cell(size(STEM_detector_angle_mrad,1),1);
for ii_stem_det = 1:size(STEM_detector_angle_mrad,1)
    stem_img{ii_stem_det} = Simu4dstem_cal_stem(dps_4d,dk,STEM_detector_angle_mrad(ii_stem_det,:));
    stem_img_file_name = [save_img_path_4_epdp,'img_STEM_',...
        num2str(STEM_detector_angle_mrad(ii_stem_det,1)),'_',...
        num2str(STEM_detector_angle_mrad(ii_stem_det,2)),'.png'];
    imwrite(mat2gray(stem_img{ii_stem_det}),stem_img_file_name);
end
save([save_mat_path_4_epdp,'STEM.mat'],'stem_img');

%% COM

if ~exist(save_img_path_4_epdp,'dir')
    mkdir(save_img_path_4_epdp);
end
if ~ exist(save_mat_path_4_epdp,'dir')
    mkdir(save_mat_path_4_epdp);
end

COM_result = Simu4dstem_cal_COM(dps_4d,dk, rotAngle);
imwrite(mat2gray(COM_result.COMx),[save_img_path_4_epdp,'img_COMx.png']);
imwrite(mat2gray(COM_result.COMy),[save_img_path_4_epdp,'img_COMy.png']);
imwrite(mat2gray(COM_result.dCOM),[save_img_path_4_epdp,'img_dCOM.png']);
imwrite(mat2gray(COM_result.iCOM),[save_img_path_4_epdp,'img_iCOM.png']);
imwrite(mat2gray(COM_result.iCOMf),[save_img_path_4_epdp,'img_iCOMf.png']);
save([save_mat_path_4_epdp,'COM.mat'],'COM_result');
% best rotate angle: on COMx png, left side of atom: white, right side of atom: black. on COMy png, upside: white, downside: black.
% use images with axis image should be the same direction.

%% prepare G for trotters
addpath('ptychoSTEM_module/');

if ~exist(save_img_path_4_epdp,'dir')
    mkdir(save_img_path_4_epdp);
end
if ~exist(save_mat_path_4_epdp,'dir')
    mkdir(save_mat_path_4_epdp);
end


ptycho = a_loop_control_parameters();
ptycho.voltage_kV  = voltage;
% Wavelength in [Angtrom]
ptycho.wavelength = Wavelength(ptycho.voltage_kV); 
% Probe convergence angle in [Radian]
ptycho.ObjApt_angle = convAngle*1e-3; 
% Rotation angle in [Radian]
ptycho.rot_angle = rotAngle*pi/180;
ptycho.pix_size = scanStep;
ptycho_0 = ptycho;

ptycho.m = dps_4d;
ptycho.probe_range = [0 0 size(ptycho.m,4) size(ptycho.m,3)];
ptycho.varfunctions.load_data = 1;
ptycho.mask_mode = 'semiangle'; % 'none' or 'semiangle'
ptycho.set_dk_rad = ''; % ''(empty) or dk value in rad

% prepare G
fprintf('%s[LOG]: calculating SSB\n',datestr(now,31));
[ptycho] = define_real_space(ptycho);
[ptycho] = detector_binning(ptycho);
[ptycho] = define_reciprocal_space(ptycho);
[ptycho] = calculate_detector_masks(ptycho);
[ptycho] = truncate_ronchigrams_within_collection_angle(ptycho);
[ptycho] = define_reciprocal_space_for_truncated_detector_plane(ptycho);
[ptycho] = FT_from_M_to_G(ptycho);


%% manul aberration:
aberr_input = [0;   % C1, in m
                0;       % C12a, in m
                0;       % C12b, in m
                0;       % C23a, in m
                0;       % C23b, in m
                0;       % C21a, in m
                0;       % C21b, in m
                0;       % C30a, in m
                0;       % C34a, in m
                0;       % C34b, in m
                0;       % C32a, in m
                0];      % C32b, in m
                
ptycho.aberr_input = aberr_input;
ptycho.aberr.C1 = aberr_input(1);
ptycho.aberr.C12a = aberr_input(2);
ptycho.aberr.C12b = aberr_input(3);
ptycho.aberr.C23a = aberr_input(4);
ptycho.aberr.C23b = aberr_input(5);
ptycho.aberr.C21a = aberr_input(6);
ptycho.aberr.C21b = aberr_input(7);
ptycho.aberr.C30a = aberr_input(8);
ptycho.aberr.C34a = aberr_input(9);
ptycho.aberr.C34b = aberr_input(10);
ptycho.aberr.C32a = aberr_input(11);
ptycho.aberr.C32b = aberr_input(12);
%% ACSSB
[ptycho]= Aberration_corrected_single_side_band_reconstruction(ptycho);
Res.ssb_r = ptycho.trotterimgRC;
Res.ssb_l = ptycho.trotterimgLC;
if isfield(ptycho,'trotterimgL_ss')
    Res.ssb_l_ss = ptycho.trotterimgL_ss;
end
if isfield(ptycho,'trotterimgR_ss')
    Res.ssb_r_ss = ptycho.trotterimgR_ss;
end
if isfield(ptycho,'func_probe')
    Res.probe = ptycho.func_probe;
end
fprintf('%s[LOG]: SSBAC finished\n',datestr(now,31));
ptycho_res = Res;

save(sprintf('%sptycho_ssbac.mat',save_mat_path_4_epdp),'ptycho_0','ptycho_res');
imwrite(mat2gray(abs(ptycho_res.ssb_r)),sprintf('%simg_SSBAC_amp_r.png',save_img_path_4_epdp));
imwrite(mat2gray(angle(ptycho_res.ssb_r)),sprintf('%simg_SSBAC_phs_r.png',save_img_path_4_epdp));
imwrite(mat2gray(abs(ptycho_res.ssb_l)),sprintf('%simg_SSBAC_amp_l.png',save_img_path_4_epdp));
imwrite(mat2gray(angle(ptycho_res.ssb_l)),sprintf('%simg_SSBAC_phs_l.png',save_img_path_4_epdp));
if isfield(ptycho_res,'ssb_l_ss')
    imwrite(mat2gray(abs(ptycho_res.ssb_l_ss)),sprintf('%simg_SSBAC_amp_l_ss.png',save_img_path_4_epdp));
    imwrite(mat2gray(angle(ptycho_res.ssb_l_ss)),sprintf('%simg_SSBAC_phs_l_ss.png',save_img_path_4_epdp));
end
if isfield(ptycho_res,'ssb_r_ss')
    imwrite(mat2gray(abs(ptycho_res.ssb_r_ss)),sprintf('%simg_SSBAC_amp_r_ss.png',save_img_path_4_epdp));
    imwrite(mat2gray(angle(ptycho_res.ssb_r_ss)),sprintf('%simg_SSBAC_phs_r_ss.png',save_img_path_4_epdp));
end

%% SSBAC masked COM
[ptycho]= Aberration_corrected_single_side_band_reconstruction_COM(ptycho);
% Res.ssb_r = ptycho.trotterimgRC;
% Res.ssb_l = ptycho.trotterimgLC;
Res.iCOM = ptycho.iCOM;
Res.dCOM = ptycho.dCOM;
Res.COMx = ptycho.COMx;
Res.COMy = ptycho.COMy;
imwrite(mat2gray(Res.COMx),[save_img_path_4_epdp,'img_SSBmCOMx.png']);
imwrite(mat2gray(Res.COMy),[save_img_path_4_epdp,'img_SSBmCOMy.png']);
imwrite(mat2gray(Res.dCOM),[save_img_path_4_epdp,'img_SSBmdCOM.png']);
imwrite(mat2gray(Res.iCOM),[save_img_path_4_epdp,'img_SSBmiCOM.png']);

% iCOM filter
dpNum = [size(dps_4d,3),size(dps_4d,4)];
fm = Simu4dstem_generate_mask('center-dark-dot',dpNum,min(dpNum).*0.05);
iCOMf = Simu4dstem_fft2(Res.iCOM);
iCOMf = abs(iCOMf) .* fm .* exp(1j .* angle(iCOMf));
Res.iCOMf = real(Simu4dstem_ifft2(iCOMf));
Res.iCOMfm = -Res.iCOMf;
imwrite(mat2gray(Res.iCOMf),[save_img_path_4_epdp,'img_SSBmiCOMf.png']);
imwrite(mat2gray(Res.iCOMfm),[save_img_path_4_epdp,'img_SSBmiCOMfm.png']);
COM_result = Res;
save(sprintf('%sCOM_ac.mat',save_mat_path_4_epdp),'COM_result');

