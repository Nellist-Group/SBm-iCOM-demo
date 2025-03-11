%% single run
Simu4dstem_main('Graphene_rec_cell_25.6x25.6A_Si_rep_t1.json');

%% loop run
df = [-1000,-500,-250,-100,0]; % in A
C23 = [0,2500,5000,10000,20000]; % in A


jsonData = jsondecode(fileread('Graphene_rec_cell_25.6x25.6A_Si_rep_aber_probe.json'));

for iidf = 1:length(df)
    for iiC23 = 1:length(C23)
        jsonData.defocus = df(iidf);
        jsonData.cond_lens_c_23_amp_A = C23(iiC23);
        jsonData.savepath_root = ['./abtest1/df_',num2str(df(iidf)/10),'nm_C23_',num2str(C23(iiC23)/10),'nm/'];
        Simu4dstem_main(jsonData);
    end
end

%% loop run 
vol = [60,300]; % kV
ModelName = {'MoS2_sl_50x50x3A';
             'MoS2_2h_50x50x10A';
             'MoS2_3R-down_50x50x10A';
             'MoS2_3R-down-up_ABA_50x50x16A';
             'MoS2_3R-rotate-up_50x50x16A';
             'MoS2_3R-up_50x50x10A';
             'MoS2_3R-up-up_ABC_50x50x16A';
             };


jsonData = jsondecode(fileread('MoS2_sl_50x50x3A_300kV_t1.json'));

for iivol = 1:length(vol)
    for iiModel = 1:length(ModelName)
        jsonData.voltage = vol(iivol);
        jsonData.model_path = ['./',ModelName{iiModel},'.mat'];
        jsonData.savepath_root = ['./t2/',ModelName{iiModel},'_',num2str(vol(iivol)),'kV_t1/'];
        Simu4dstem_main(jsonData);
        if exist([jsonData.savepath_root,'dps_epdp_0.mat'],'file') == 2
            delete([jsonData.savepath_root,'dps_epdp_0.mat']);
        end
    end
end