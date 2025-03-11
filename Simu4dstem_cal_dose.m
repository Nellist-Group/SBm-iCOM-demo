function doseSeries = Simu4dstem_cal_dose(jsonData)
%SIMU4DSTEM_CAL_DOSE Summary of this function goes here
%   Detailed explanation goes here

% scanSpotNum = jsonData.array_size_x .* jsonData.array_size_y;
epdp = jsonData.electron_num_per_dp; % electron per diffraction pattern
% scanArea = jsonData.arra_size_x .* jsonData.arra_size_x .* jsonData.step_size.^2; % in A2
% doseSeries = scanSpotNum .* epdp ./ scanArea;

doseSeries = epdp ./ jsonData.step_size .^ 2; % e/A2



end

