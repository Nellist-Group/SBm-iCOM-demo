function [EW,IW] = Simu4dstem_mul_slice_prop_forward(probe, obj_cut, sliceloc, lambda, dxy)
%EPIE_MUL_SLICE_PROP_FORWARD multi-slice propagation forward
%   sliceloc in A
%   lambda in m
%   dxy in A
CurrentLoc = 0; % current probe plane, in A; as the defocus value of the probe;
IW = cell(1,length(obj_cut)); % IW: incident wavefunction before each slice
EW = cell(1,length(obj_cut)); % EW: exit wave function after each slice
EWcl = probe; % exit wave at current location;
for ii = 1:length(sliceloc)
    PropDistance = sliceloc(ii) - CurrentLoc; % in A
    EWcl = Simu4dstem_fresnel_prop_TF(EWcl, lambda, PropDistance .* 1e-10, [dxy,dxy].*1e-10);
    IW{ii} = EWcl;
    EWcl = EWcl .* obj_cut{ii};
    EW{ii} = EWcl;
    CurrentLoc = sliceloc(ii);
end

end


