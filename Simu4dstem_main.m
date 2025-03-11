function Simu4dstem_main(varargin)
% Simu4dstem_main
% input instructions:
%   if no input parameter, then open a ui to get json file
%   if 1 input parameter, 
%       if it is a string read as the file path of json parameters
%       if it is a struct, read as jsonData
%   if 2 input parameter, read 1st input as json file path, 2nd input as
%       ePIE save prefix
%% check input
if nargin == 1
    if isa(varargin{1},'struct')
        jsonData = varargin{1};
        json_name = 'parameters.json';
    elseif isa(varargin{1},'char')
        jsonPath_input = varargin{1};
        [json_path,json_name,json_ext] = fileparts(jsonPath_input);
        json_name = [json_name,json_ext];
        jsonData = jsondecode(fileread(jsonPath_input));
    end
    ePIE_save_prefix = '';
elseif nargin == 2
    if isa(varargin{1},'struct')
        jsonData = varargin{1};
        json_name = 'parameters.json';
    elseif isa(varargin{1},'char')
        jsonPath_input = varargin{1};
        [json_path,json_name,json_ext] = fileparts(jsonPath_input);
        json_name = [json_name,json_ext];
        jsonData = jsondecode(fileread(jsonPath_input));
    end
    ePIE_save_prefix = varargin{2};
elseif nargin == 0
    [json_name, json_path] = uigetfile('*.json');
    jsonPath_input = [json_path,json_name];
    ePIE_save_prefix = '';
    jsonData = jsondecode(fileread(jsonPath_input));
end

%%
% parameter check
if jsonData.nx ~= jsonData.ny
    error('jsonData: nx ~= ny');
end
if ~exist(jsonData.savepath_root,'dir')
    mkdir(jsonData.savepath_root);
end
json_str = jsonencode(jsonData);
json_save_id = fopen([jsonData.savepath_root,filesep,json_name],'w');
fprintf(json_save_id,json_str);
fclose(json_save_id);
%copyfile(jsonPath_input,[jsonData.savepath_root,filesep,json_name]);

% parameters calculation
lambda = 1.23984244/sqrt(jsonData.voltage*(2*510.99906+jsonData.voltage))*1e-9;% m
dxy = (lambda*1e10) / jsonData.nx / (jsonData.dk*1e-3); % in A

% generate probe
if jsonData.process_ctrl.do_generate_probe == 1
    probe = Simu4dstem_gen_probe_for_multem(jsonData);
    save([jsonData.savepath_root,'probe.mat'],'probe');
elseif jsonData.process_ctrl.do_generate_probe == 2
    load([jsonData.savepath_root,'probe.mat'],'probe');
end

% disp
fprintf('%s[info]Simu4dstem: probe calculated\n',datestr(now,31));

% calculate object transmission function (tf)
if jsonData.process_ctrl.do_cal_obj_tf == 1
    % process atoms
    [atoms,lx,ly,lz] = Simu4dstem_preprocess_atoms(jsonData);
    % lx, ly, lz in A
    
    if jsonData.use_gpu == 0
        device_flag_multem = 1;
    elseif jsonData.use_gpu == 1
        device_flag_multem = 2;
    end
    % calculate object function (transmission function) by multislice using multem
    [obj_multem_strc,Slice] = Simu4dstem_cal_obj_fun(atoms,lx,ly,lz,...
        round(lx/dxy),round(ly/dxy),...
        jsonData.model_slice_mode,jsonData.dz,...
        jsonData.voltage,device_flag_multem);
    save([jsonData.savepath_root,'tf.mat'],'obj_multem_strc','Slice','-v7.3');
    
    % disp
    fprintf('%s[info]Simu4dstem: object transmission function calculated (%d layers)\n',datestr(now,31),size(Slice,1));
elseif jsonData.process_ctrl.do_cal_obj_tf == 2
    load([jsonData.savepath_root,'tf.mat'],'obj_multem_strc','Slice');
    fprintf('%s[info]Simu4dstem: object transmission function loaded (%d layers)\n',datestr(now,31),size(Slice,1));
end

% scan 4dstem dataset
if jsonData.process_ctrl.do_scan_4dstem == 1
    % convert multem result to scanning tf
    obj_tf = cell(1,length(obj_multem_strc));
    for ii = 1:length(obj_multem_strc)
        obj_tf{ii} = obj_multem_strc{ii}.trans;
    end
    clear('obj_multem_strc');
    sliceloc = (Slice(:,1) + Slice(:,2))./2;
    
    % get scanning trans_A (trans in A)
    trans_A = Simu4dstem_generate_trans(...
        [jsonData.array_size_x,jsonData.array_size_y],...
        [jsonData.step_size,jsonData.step_size],...
        [jsonData.extra_pxl_num_on_atom_pot*dxy + jsonData.nx/2*dxy,jsonData.extra_pxl_num_on_atom_pot*dxy + jsonData.ny/2*dxy],...
        [0,0]);
    trans_pxl = trans_A ./ dxy;
    save([jsonData.savepath_root,'trans_pos.mat'],'trans_A','trans_pxl','dxy','-v7.3');
    
    % init dps
    dps = zeros(size(probe,1),size(probe,2),size(trans_A,1));
    if jsonData.use_gpu == 1
        dps = gpuArray(dps);
        probe = gpuArray(probe);
        for count_tf_mv = 1:length(obj_tf)
            obj_tf{count_tf_mv} = gpuArray(obj_tf{count_tf_mv});
        end
    end
    
    for count_scan = 1:size(trans_A)
        % get obj_cut
        obj_cut = cell(1,length(obj_tf));
        cut_range_x = [round(trans_pxl(count_scan,1))-floor(size(probe,1)/2),...
            round(trans_pxl(count_scan,1))+ceil(size(probe,1)/2)-1];
        cut_range_y = [round(trans_pxl(count_scan,2))-floor(size(probe,2)/2),...
            round(trans_pxl(count_scan,2))+ceil(size(probe,2)/2)-1];
        for ii = 1:length(obj_tf)
            obj_cut{ii} = obj_tf{ii}(cut_range_x(1):cut_range_x(2),cut_range_y(1):cut_range_y(2));
        end
        
        % calculate exit wave
        [exit_wave, incident_wave] = Simu4dstem_mul_slice_prop_forward(probe,obj_cut,sliceloc,lambda,dxy);
        
        exit_wave_final = exit_wave{end};
        % select-area aperture
        if jsonData.select_area_aperture_radius_in_nm ~= 0
            sa_aperture = Simu4dstem_generate_mask('center-dot',size(exit_wave{end}),...
                round(jsonData.select_area_aperture_radius_in_nm./(dxy/10)));
            exit_wave_final = exit_wave_final .* sa_aperture;
        end
        
        % calculate diffraction patterns
        dps(:,:,count_scan) = abs(Simu4dstem_fft2(exit_wave_final)).^2;
    end
    
    % gather from gpu
    if jsonData.use_gpu == 1
        dps = gather(dps);
    end
    dps0 = dps;
    
    % disp
    fprintf('%s[info]Simu4dstem: scan finished\n',datestr(now,31));
    
    % set inelastic simulation
    if jsonData.inelastic_distribution == 1
        fprintf('%s[info]Simu4dstem: setting energy window\n',datestr(now,31));
        load(jsonData.EELS_path,'eels');
        save([jsonData.savepath_root,'EELS_spec.mat'],'eels');
        eels_win = Simu4dstem_generate_energy_window(eels,jsonData.dE,jsonData.energy_window);
        for iiew = 1:size(jsonData.energy_window,1)
            fprintf('%s[info]Simu4dstem: calculating inelastic distribution for energy window [%d,%d]eV\n',datestr(now,31),jsonData.energy_window(iiew,1),jsonData.energy_window(iiew,2));
            dps0 = Simu4dstem_add_inelastic_distribution(...
                dps0, jsonData.dk, eels_win{1,iiew}, eels_win{2,iiew}, jsonData.dE, jsonData.voltage);
            % save
            fprintf('%s[info]Simu4dstem: saving dps for energy window [%d,%d]eV\n',datestr(now,31),jsonData.energy_window(iiew,1),jsonData.energy_window(iiew,2));
            save([jsonData.savepath_root,'dps_ew[',num2str(jsonData.energy_window(iiew,1)),',',num2str(jsonData.energy_window(iiew,2)),']_epdp_0.mat'],'dps','-v7.3');
            % set dose by poisson distribution
            for ii = 1:length(jsonData.electron_num_per_dp)
                if jsonData.electron_num_per_dp(ii) ~= 0
                    dps = Simu4dstem_add_poisson_noise(dps0, jsonData.electron_num_per_dp(ii) .* eels_win{3,iiew});
                    save([jsonData.savepath_root,'dps_ew[',num2str(jsonData.energy_window(iiew,1)),',',num2str(jsonData.energy_window(iiew,2)),']_epdp_',num2str(jsonData.electron_num_per_dp(ii)),'.mat'],'dps','-v7.3');
                end
            end
        end
	save([jsonData.savepath_root,'eels_win.mat'],'eels_win','-v7.3');
    else % not inelastic mode, only poisson noise
        % save
        fprintf('%s[info]Simu4dstem: saving dps\n',datestr(now,31));
        save([jsonData.savepath_root,'dps_epdp_0.mat'],'dps','-v7.3');
        % set dose by poisson distribution
        for ii = 1:length(jsonData.electron_num_per_dp)
            if jsonData.electron_num_per_dp(ii) ~= 0
                dps = Simu4dstem_add_poisson_noise(dps0, jsonData.electron_num_per_dp(ii));
                save([jsonData.savepath_root,'dps_epdp_',num2str(jsonData.electron_num_per_dp(ii)),'.mat'],'dps','-v7.3');
            end
        end
        fprintf('%s[info]Simu4dstem: dps for epdp saved\n',datestr(now,31));
        
    end
end

%% reconstruct
if sum([jsonData.process_ctrl.do_STEM_img,jsonData.process_ctrl.do_COM_img,...
        jsonData.process_ctrl.do_SSB_img,jsonData.process_ctrl.do_SSBAC_img,...
        jsonData.process_ctrl.do_WDD_img,jsonData.process_ctrl.do_SSBAC_mask_COM,...
        jsonData.process_ctrl.do_ePIE_img,jsonData.process_ctrl.do_SSB_mask_COM,...
        jsonData.process_ctrl.do_WDD_deconv_iCOM,jsonData.process_ctrl.do_DPC_img]) > 0
if jsonData.inelastic_distribution == 1
    ewLoopNum = size(jsonData.energy_window,1);
else
    ewLoopNum = 1;
end
for iiew = 1:ewLoopNum
    if jsonData.inelastic_distribution == 1
        ew_name = ['_ew[',num2str(jsonData.energy_window(iiew,1)),',',num2str(jsonData.energy_window(iiew,2)),']'];
    else
        ew_name = '';
    end
for ii = 1:length(jsonData.electron_num_per_dp)
    % load dp
    load([jsonData.savepath_root,'dps',ew_name,'_epdp_',num2str(jsonData.electron_num_per_dp(ii)),'.mat'],'dps');
    
    % create image folder
    save_img_path_4_epdp = [jsonData.savepath_root,'result_img_epdp',num2str(jsonData.electron_num_per_dp(ii)),'/'];
    save_mat_path_4_epdp = [jsonData.savepath_root,'result_mat_epdp',num2str(jsonData.electron_num_per_dp(ii)),'/'];
    if ~exist(save_img_path_4_epdp,'dir')
        mkdir(save_img_path_4_epdp);
    end
    if ~exist(save_mat_path_4_epdp,'dir')
        mkdir(save_mat_path_4_epdp);
    end
    
    % get dps_4d
    dps_4d = reshape(dps,[size(dps,1),size(dps,2),jsonData.array_size_x,jsonData.array_size_y]);
    
    % STEM image
    if jsonData.process_ctrl.do_STEM_img
        stem_img = cell(size(jsonData.STEM_detector_angle_mrad,1),1);
        if jsonData.spot_size_A ~= 0
            stem_img_spot_size = cell(size(jsonData.STEM_detector_angle_mrad,1),1);
            spot_size_pxl = jsonData.spot_size_A ./ jsonData.step_size;
        end
        for ii_stem_det = 1:size(jsonData.STEM_detector_angle_mrad,1)
            stem_img{ii_stem_det} = Simu4dstem_cal_stem(dps_4d,jsonData.dk,jsonData.STEM_detector_angle_mrad(ii_stem_det,:));
            stem_img_file_name = [save_img_path_4_epdp,'img_STEM',ew_name,'_',...
                num2str(jsonData.STEM_detector_angle_mrad(ii_stem_det,1)),'_',...
                num2str(jsonData.STEM_detector_angle_mrad(ii_stem_det,2)),'.png'];
            imwrite(mat2gray(stem_img{ii_stem_det}),stem_img_file_name);
            if jsonData.spot_size_A ~= 0
                stem_img_spot_size{ii_stem_det} = Simu4dstem_cal_stem(dps_4d,jsonData.dk,jsonData.STEM_detector_angle_mrad(ii_stem_det,:));
                stem_img_spot_size{ii_stem_det} = Simu4dstem_conv_spot_size(stem_img_spot_size{ii_stem_det},spot_size_pxl);
                stem_img_spot_size_file_name = [save_img_path_4_epdp,'img_STEM',ew_name,'_',...
                        num2str(jsonData.STEM_detector_angle_mrad(ii_stem_det,1)),'_',...
                        num2str(jsonData.STEM_detector_angle_mrad(ii_stem_det,2)),'_ss',...
                        num2str(jsonData.spot_size_A),'A.png'];
                imwrite(mat2gray(stem_img_spot_size{ii_stem_det}),stem_img_spot_size_file_name);
            end
        end
        save([save_mat_path_4_epdp,'STEM',ew_name,'.mat'],'stem_img');
        if jsonData.spot_size_A ~= 0
            save([save_mat_path_4_epdp,'STEM_spot_size',ew_name,'.mat'],'stem_img_spot_size');
        end
    end
    
    % COM image
    if jsonData.process_ctrl.do_COM_img
        COM_result = Simu4dstem_cal_COM(dps_4d,jsonData.dk, 0);
        imwrite(mat2gray(COM_result.COMx),[save_img_path_4_epdp,'img_COMx',ew_name,'.png']);
        imwrite(mat2gray(COM_result.COMy),[save_img_path_4_epdp,'img_COMy',ew_name,'.png']);
        imwrite(mat2gray(COM_result.dCOM),[save_img_path_4_epdp,'img_dCOM',ew_name,'.png']);
        imwrite(mat2gray(COM_result.iCOM),[save_img_path_4_epdp,'img_iCOM',ew_name,'.png']);
        imwrite(mat2gray(abs(Simu4dstem_fft2(COM_result.iCOM))),[save_img_path_4_epdp,'img_iCOM_fft',ew_name,'.png']);
        save([save_mat_path_4_epdp,'COM',ew_name,'.mat'],'COM_result');
    end
    
    % DPC image
    if jsonData.process_ctrl.do_DPC_img
        COM_result = Simu4dstem_cal_DPC(dps_4d,jsonData.dk, 0, jsonData.DPC_detector_angle_mrad(1),jsonData.DPC_detector_angle_mrad(2));
        imwrite(mat2gray(COM_result.DPCx),[save_img_path_4_epdp,'img_DPCx',ew_name,'.png']);
        imwrite(mat2gray(COM_result.DPCy),[save_img_path_4_epdp,'img_DPCy',ew_name,'.png']);
        imwrite(mat2gray(COM_result.dDPC),[save_img_path_4_epdp,'img_dDPC',ew_name,'.png']);
        imwrite(mat2gray(COM_result.iDPC),[save_img_path_4_epdp,'img_iDPC',ew_name,'.png']);
        imwrite(mat2gray(abs(Simu4dstem_fft2(COM_result.iDPC))),[save_img_path_4_epdp,'img_iDPC_fft',ew_name,'.png']);
        save([save_mat_path_4_epdp,'DPC',ew_name,'.mat'],'COM_result');
    end
    
    
    % SSB mask COM
    if jsonData.process_ctrl.do_SSB_mask_COM
        [ptycho_0, ptycho_res] = Simu4dstem_cal_SSB_mask_COM(jsonData,dps_4d);
        save(sprintf('%sptycho_ssb_mask_COM%s.mat',save_mat_path_4_epdp,ew_name),'ptycho_0','ptycho_res');
%         imwrite(mat2gray(ptycho_res.ssb_r),sprintf('%simg_SSB_COM_amp_r%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(ptycho_res.ssb_r),sprintf('%simg_SSB_COM_phs_r%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(ptycho_res.ssb_r))),sprintf('%simg_SSB_COM_phs_r_fft%s.png',save_img_path_4_epdp,ew_name));
%         imwrite(mat2gray(ptycho_res.ssb_l),sprintf('%simg_SSB_COM_amp_l%s.png',save_img_path_4_epdp,ew_name));
%         imwrite(mat2gray(ptycho_res.ssb_l),sprintf('%simg_SSB_COM_phs_l%s.png',save_img_path_4_epdp,ew_name));
        if isfield(ptycho_res,'ssb_l_ss')
%             imwrite(mat2gray(ptycho_res.ssb_l_ss),sprintf('%simg_SSB_COM_amp_l_ss%s.png',save_img_path_4_epdp,ew_name));
%             imwrite(mat2gray(ptycho_res.ssb_l_ss),sprintf('%simg_SSB_COM_phs_l_ss%s.png',save_img_path_4_epdp,ew_name));
        end
        if isfield(ptycho_res,'ssb_r_ss')
%             imwrite(mat2gray(ptycho_res.ssb_r_ss),sprintf('%simg_SSB_COM_amp_r_ss%s.png',save_img_path_4_epdp,ew_name));
%             imwrite(mat2gray(ptycho_res.ssb_r_ss),sprintf('%simg_SSB_COM_phs_r_ss%s.png',save_img_path_4_epdp,ew_name));
        end
    end
    
        % SSBAC mask COM
    if jsonData.process_ctrl.do_SSBAC_mask_COM
        [ptycho_0, ptycho_res] = Simu4dstem_cal_SSBAC_mask_COM(jsonData,dps_4d);
        save(sprintf('%sptycho_ssbac_mask_COM%s.mat',save_mat_path_4_epdp,ew_name),'ptycho_0','ptycho_res');
        imwrite(mat2gray(ptycho_res.iCOM),sprintf('%simg_SSBAC_iCOM%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(ptycho_res.iCOM))),sprintf('%simg_SSBAC_iCOM_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(ptycho_res.dCOM),sprintf('%simg_SSBAC_dCOM%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(ptycho_res.dCOM))),sprintf('%simg_SSBAC_dCOM_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(ptycho_res.COMx),sprintf('%simg_SSBAC_COMx%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(ptycho_res.COMx))),sprintf('%simg_SSBAC_COMx_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(ptycho_res.COMy),sprintf('%simg_SSBAC_COMy%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(ptycho_res.COMy))),sprintf('%simg_SSBAC_COMy_fft%s.png',save_img_path_4_epdp,ew_name));
    end
    
    % WDD deconvolute iCOM
    if jsonData.process_ctrl.do_WDD_deconv_iCOM
        [ptycho_0, ptycho_res] = Simu4dstem_cal_WDD_deconv_iCOM(jsonData,dps_4d);
        save(sprintf('%sptycho_WDD_deconv_iCOM%s.mat',save_mat_path_4_epdp,ew_name),'ptycho_0','ptycho_res');
        imwrite(mat2gray(ptycho_res.wdd),sprintf('%simg_WDD_deconv_iCOM%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(ptycho_res.wdd))),sprintf('%simg_WDD_deconv_iCOM_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(ptycho_res.wdd_probe)),sprintf('%sprobe_WDD_deconv_iCOM_amp%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.wdd_probe)),sprintf('%sprobe_WDD_deconv_iCOM_phs%s.png',save_img_path_4_epdp,ew_name));
    end
    
    % SSB (automatically addpath('ptychoSTEM_module/'))
    if jsonData.process_ctrl.do_SSB_img
        [ptycho_0, ptycho_res] = Simu4dstem_cal_SSB(jsonData,dps_4d);
        save(sprintf('%sptycho_ssb%s.mat',save_mat_path_4_epdp,ew_name),'ptycho_0','ptycho_res');
        imwrite(mat2gray(abs(ptycho_res.ssb_r)),sprintf('%simg_SSB_amp_r%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.ssb_r)),sprintf('%simg_SSB_phs_r%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(angle(ptycho_res.ssb_r)))),sprintf('%simg_SSB_phs_r_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(ptycho_res.ssb_l)),sprintf('%simg_SSB_amp_l%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.ssb_l)),sprintf('%simg_SSB_phs_l%s.png',save_img_path_4_epdp,ew_name));
        if isfield(ptycho_res,'ssb_l_ss')
            imwrite(mat2gray(abs(ptycho_res.ssb_l_ss)),sprintf('%simg_SSB_amp_l_ss%s.png',save_img_path_4_epdp,ew_name));
            imwrite(mat2gray(angle(ptycho_res.ssb_l_ss)),sprintf('%simg_SSB_phs_l_ss%s.png',save_img_path_4_epdp,ew_name));
        end
        if isfield(ptycho_res,'ssb_r_ss')
            imwrite(mat2gray(abs(ptycho_res.ssb_r_ss)),sprintf('%simg_SSB_amp_r_ss%s.png',save_img_path_4_epdp,ew_name));
            imwrite(mat2gray(angle(ptycho_res.ssb_r_ss)),sprintf('%simg_SSB_phs_r_ss%s.png',save_img_path_4_epdp,ew_name));
        end
    end
    
        % SSBAC (automatically addpath('ptychoSTEM_module/'))
    if jsonData.process_ctrl.do_SSBAC_img
        [ptycho_0, ptycho_res] = Simu4dstem_cal_SSBAC(jsonData,dps_4d);
        save(sprintf('%sptycho_ssbac%s.mat',save_mat_path_4_epdp,ew_name),'ptycho_0','ptycho_res');
        imwrite(mat2gray(abs(ptycho_res.ssb_r)),sprintf('%simg_SSBAC_amp_r%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.ssb_r)),sprintf('%simg_SSBAC_phs_r%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(angle(ptycho_res.ssb_r)))),sprintf('%simg_SSBAC_phs_r_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(ptycho_res.ssb_l)),sprintf('%simg_SSBAC_amp_l%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.ssb_l)),sprintf('%simg_SSBAC_phs_l%s.png',save_img_path_4_epdp,ew_name));
        if isfield(ptycho_res,'ssb_l_ss')
            imwrite(mat2gray(abs(ptycho_res.ssb_l_ss)),sprintf('%simg_SSBAC_amp_l_ss%s.png',save_img_path_4_epdp,ew_name));
            imwrite(mat2gray(angle(ptycho_res.ssb_l_ss)),sprintf('%simg_SSBAC_phs_l_ss%s.png',save_img_path_4_epdp,ew_name));
        end
        if isfield(ptycho_res,'ssb_r_ss')
            imwrite(mat2gray(abs(ptycho_res.ssb_r_ss)),sprintf('%simg_SSBAC_amp_r_ss%s.png',save_img_path_4_epdp,ew_name));
            imwrite(mat2gray(angle(ptycho_res.ssb_r_ss)),sprintf('%simg_SSBAC_phs_r_ss%s.png',save_img_path_4_epdp,ew_name));
        end
    end
    
        % WDD (automatically addpath('ptychoSTEM_module/'))
    if jsonData.process_ctrl.do_WDD_img
        [ptycho_0, ptycho_res] = Simu4dstem_cal_WDD(jsonData,dps_4d);
        save(sprintf('%sptycho_wdd%s.mat',save_mat_path_4_epdp,ew_name),'ptycho_0','ptycho_res');
        imwrite(mat2gray(abs(ptycho_res.wdd)),sprintf('%simg_WDD_amp%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.wdd)),sprintf('%simg_WDD_phs%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(Simu4dstem_fft2(angle(ptycho_res.wdd)))),sprintf('%simg_WDD_phs_fft%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(abs(ptycho_res.wdd_probe)),sprintf('%sprobe_WDD_amp%s.png',save_img_path_4_epdp,ew_name));
        imwrite(mat2gray(angle(ptycho_res.wdd_probe)),sprintf('%sprobe_WDD_phs%s.png',save_img_path_4_epdp,ew_name));
    end
    
    % ePIE (automatically addpath('ePIE_module/'))
    if jsonData.process_ctrl.do_ePIE_img
        ePIE_result = Simu4dstem_cal_ePIE(jsonData,dps,[save_mat_path_4_epdp,'ePIE',ew_name,ePIE_save_prefix,'/']);
        save_range_ePIE_X = [floor(min(ePIE_result.trans_exec.x)+ ePIE_result.trans_exec.dp_size/2),ceil(max(ePIE_result.trans_exec.x)+ ePIE_result.trans_exec.dp_size/2)];
        save_range_ePIE_Y = [floor(min(ePIE_result.trans_exec.y)+ ePIE_result.trans_exec.dp_size/2),ceil(max(ePIE_result.trans_exec.y)+ ePIE_result.trans_exec.dp_size/2)];
        ePIE_obj_scan_area = ePIE_result.obj(save_range_ePIE_X(1):save_range_ePIE_X(2),save_range_ePIE_Y(1):save_range_ePIE_Y(2));
        imwrite(mat2gray(abs(ePIE_obj_scan_area)),[save_img_path_4_epdp,'img_ePIE_amp',ew_name,ePIE_save_prefix,'.png']);
        imwrite(mat2gray(angle(ePIE_obj_scan_area)),[save_img_path_4_epdp,'img_ePIE_phs',ew_name,ePIE_save_prefix,'.png']);
        imwrite(mat2gray(abs(ePIE_result.probe{1})),[save_img_path_4_epdp,'probe_ePIE_amp',ew_name,ePIE_save_prefix,'.png']);
        imwrite(mat2gray(angle(ePIE_result.probe{1})),[save_img_path_4_epdp,'probe_ePIE_phs',ew_name,ePIE_save_prefix,'.png']);
    end
end
end
end

end