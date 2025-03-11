function ePIE_result = Simu4dstem_cal_ePIE(jsonData,dps,mat_save_path)
%SIMU4DSTEM_CAL_EPIE Summary of this function goes here
%   Detailed explanation goes here
addpath('ePIE_module');
p = ePIE_parameters();
p.datapath = dps;
p.dataformat = 'var-3d';
p.savepath = mat_save_path;
p.beam_voltage = jsonData.voltage;
p.semiangle = jsonData.alpha;
p.dp.num_all = jsonData.array_size_x*jsonData.array_size_y;
if jsonData.nx ~= jsonData.ny
    error('Simu4dstem_cal_ePIE: jsonData.nx ~= jsonData.ny');
end
p.dp.size_raw = jsonData.nx;
p.dp.dk = jsonData.dk;
p.dp.padding_both_side_by_pxl = 0; % a factor of padding. the size of dp used in ePIE_iteration will be p.dp.size + 2*padding_both_side_by_pxl
% p.switch.gpu = true;
p.switch.show_result = false;
p.switch.save_result = true;
p.switch.iteration_mask_mode = 'threshold'; % 'none'/'semiangle'/'threshold'/'semiangle-hollow'/'mask-file'.
p.switch.save_result_every_x_loop = 0; % save result every x loops. 0 means no saving during iterations
p.switch.iteration_mask_threshold = 1;
p.trans.scan_num = [jsonData.array_size_y,jsonData.array_size_x];
p.trans.scan_step = [jsonData.step_size, jsonData.step_size]; % scan step in A
p.trans.transpose = true; % do a transpose to trans, mean x and y exchange, scan y axis first
p.trans.draw_trans = false;
p.iter = 5;
p.update_function.algorithm = 'ePIE';
p.recon.guess_init_probe = 1;
p.statics.init_probe = [jsonData.savepath_root,'init_probe.mat'];
p.recon.probe_aber_constrain = 0;
% calculate probe aberation parameters from reconsturcted probe, 
% then use the paramter to re-generate probe.
% 0 --> no
% x (x>0) --> do the calculate and probe replacement every x updates of probe
p.recon.probe_aber_fix_list = ... % only valid if probe_aber_constrain>0, nan: calculate this aberation. value: fix this aberration
[   nan; % C10, in m
    nan; % C12a, in m
    nan; % C12b, in m
    nan; % C23a, in m
    nan; % C23b, in m
    nan; % C21a, in m
    nan; % C21b, in m
    nan; % C30, in m
    nan; % C34a, in m
    nan; % C34b, in m
    nan; % C32a, in m
    nan; % C32b, in m
    nan; % phase constant.
    nan; % x shift
    nan]; % y shift

%----- Dynamic Probe -----%
p.dynaP.enable = false;
p.dynaP.mode = 'square-area';
p.dynaP.para = {10};
% mode: 'square-are'
%       divide scanning area into square areas, use the same probe in each
%       area
%       para: {x}; x--> side length in A of each small area

%----- Parallel -----%
p.para.enable = false;
p.para.mode = 'sub-array-hard-edge';
p.para.mode_parameter = {2,2};
% 'sub-array-soft-edge': divide diffraction patterns to sub-arrays with a changing edge.
%   soft edge is only suitable for dynaP=false
%   parameters: {number of subarray (x), number of subarray (y)}
% 'sub-array-hard-edge': divide diffraction patterns to sub-arrays with a fixed edge
%   suitable for dyanP=true. if dyanP mode='suqare-area', then use that as
%   the sub-array setting

p.coeffi(1,:) = [jsonData.defocus/10, 0];%nm      C10     df
p.coeffi(5,:) = [jsonData.cond_lens_c_30_mm*1e3, 0];%    um      C3      C3
if jsonData.cond_lens_c_12_amp_A == 0
    p.coeffi(2,:) = [0,0];
else
    p.coeffi(2,:) = [jsonData.cond_lens_c_12_amp_A/10, jsonData.cond_lens_c_12_angle];%    nm      C12     A1
end
if jsonData.cond_lens_c_23_amp_A == 0
    p.coeffi(4,:) = [0,0];
else
    p.coeffi(4,:) = [jsonData.cond_lens_c_23_amp_A/10, jsonData.cond_lens_c_23_angle];%    nm      C23     A2
end

% gpu
if jsonData.use_gpu == 0
    p.switch.gpu = false;
elseif jsonData.use_gpu == 1
    p.switch.gpu = true;
end

[dps,probe,obj,trans_exec,p,reconTmp]=ePIE_preprocess(p);
ePIE_result = ePIE_main(dps,probe,obj,trans_exec,p,reconTmp);

% spot size convolution
if jsonData.spot_size_A > 0
    ePIE_result.obj = Simu4dstem_spot_size_conv(ePIE_result.obj,ePIE_result.dxy,jsonData.spot_size_A);
end

% result = ePIE_result;
% pre_str = '';

% if ~exist(p.savepath,'dir')
%     mkdir(p.savepath);
% end
% save([p.savepath,'/',pre_str,'result.mat'],'result');
% if strcmp(p.dataformat,'var-3d')
%     p.datapath = [];
% end
% save([p.savepath,'/',pre_str,'p.mat'],'p');
% 
% [~] = ePIE_show_wf(result.obj,'Save',true,...
%     'SaveNameStartWith',sprintf('%s/%sobj',p.savepath,pre_str),...
%     'Trans',result.trans_exec,'DrawFigure',false);
% [~] = ePIE_show_wf(result.probe{1},'Save',true,...
%     'SaveNameStartWith',sprintf('%s/%sprobe',p.savepath,pre_str),...
%     'DrawFigure',false);

end

