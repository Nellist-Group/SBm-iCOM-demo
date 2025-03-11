function Simu4dstem_COM_contrast_scan(dps_4d, dk, sample_num)
% calculate COM_contrast of a series of rotation angle, and
% plot the contrast image
% input: sample_num: the number of rotation angle calculated
AngleList = 360./sample_num .* (0:sample_num-1);
res_list = Simu4dstem_cal_COM(dps_4d, dk, AngleList);
dCOM_range_list = zeros(length(AngleList),1);
iCOM_range_list = zeros(length(AngleList),1);
for ii = 1:length(AngleList)
    dCOM_range_list(ii) = max(res_list(ii).dCOM(:)) - min(res_list(ii).dCOM(:));
    iCOM_range_list(ii) = max(res_list(ii).iCOM(:)) - min(res_list(ii).iCOM(:));
end
figure();
subplot(1,2,1);
plot(AngleList, dCOM_range_list);
title('dCOM max-min'); xlabel('angle');
subplot(1,2,2);
plot(AngleList, iCOM_range_list);
title('iCOM max-min'); xlabel('angle');
end

