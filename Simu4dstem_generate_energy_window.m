function eels_win = Simu4dstem_generate_energy_window(eels,dE,energy_window)
%SIMU4DSTEM_GENERATE_ENERGY_WINDOW generate energy windows by eels, eD and
% window settings
%   input eels starts from 0 eV
%   eels_win{1,:} --> windowed eels spectrums
%   eels_win{2,:} --> energy loss of 1st channel in windowed spectrum in eV
%   eels_win{3,:} --> counts ratio of windowed spectrums and all counts

deltaE = (0:length(eels)-1) .* dE;
eels_win = cell(3,size(energy_window,1));

for ii = 1:size(energy_window,1)
    indts = find(deltaE >= energy_window(ii,1),1);
    indte = find(deltaE > energy_window(ii,2),1)-1;
    if isempty(indte)
        indte = length(deltaE);
    end
    eels_win{1,ii} = eels(indts:indte);
    eels_win{2,ii} = deltaE(indts);
    eels_win{3,ii} = sum(eels_win{1,ii}) ./ sum(eels);
end


end

