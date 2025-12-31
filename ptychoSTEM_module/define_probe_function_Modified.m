function [ptycho] = define_probe_function_Modified(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% this function needs comments

% ptycho.varfunctions.define_probe_function: Flag to indicate the function
% has been executed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  probe function
% Set up the lens aberrations:
% Lens Aberrations (sorted by radial order):
%     aberr_input(1) = ptycho.aberr.C1;
%     aberr_input(2) = ptycho.aberr.C12a;
%     aberr_input(3) = ptycho.aberr.C12b;
%     aberr_input(4) = ptycho.aberr.C23a; 
%     aberr_input(5) = ptycho.aberr.C23b; 
%     aberr_input(6) = ptycho.aberr.C21a; 
%     aberr_input(7) = ptycho.aberr.C21b; 
%     aberr_input(8) = ptycho.aberr.C30a;
%     aberr_input(9) = ptycho.aberr.C34a;
%     aberr_input(10) = ptycho.aberr.C34b;
%     aberr_input(11) = ptycho.aberr.C32a;
%     aberr_input(12) = ptycho.aberr.C32b;
    
    aberr_input(1) = ptycho.aberr.C1;
    aberr_input(2) = real(ptycho.aberr.C12a.*exp(1i.*ptycho.aberr.C12b));
    aberr_input(3) = imag(ptycho.aberr.C12a.*exp(1i.*ptycho.aberr.C12b));
    aberr_input(4) = real(ptycho.aberr.C23a.*exp(1i.*ptycho.aberr.C23b));
    aberr_input(5) = imag(ptycho.aberr.C23a.*exp(1i.*ptycho.aberr.C23b));
    aberr_input(6) = real(ptycho.aberr.C21a.*exp(1i.*ptycho.aberr.C21b));
    aberr_input(7) = imag(ptycho.aberr.C21a.*exp(1i.*ptycho.aberr.C21b));
    aberr_input(8) = ptycho.aberr.C30a;
    aberr_input(9) = real(ptycho.aberr.C34a.*exp(1i.*ptycho.aberr.C34b));
    aberr_input(10) = imag(ptycho.aberr.C34a.*exp(1i.*ptycho.aberr.C34b));
    aberr_input(11) = real(ptycho.aberr.C32a.*exp(1i.*ptycho.aberr.C32b));
    aberr_input(12) = imag(ptycho.aberr.C32a.*exp(1i.*ptycho.aberr.C32b));
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% probe function
% chi function
func_aber = FuncAberrUV(ptycho.Kx_wp,ptycho.Ky_wp,aberr_input);

% transfer function
func_transfer=exp(-1i*2*pi/ (ptycho.wavelength*1e-10) .* func_aber);

% aperture function
func_ObjApt = ones(size(ptycho.Kwp));
func_ObjApt( ptycho.Kwp > ptycho.ObjApt_angle) = 0;

% dose equals to the summed intensity of the average ronchigram.
dose = sum(ptycho.pacbed(:)) *1.0;
% for resampled ronchigram
% pacbed_rs = interp2(tx_wp,ty_wp,mean_m_wp,Kx,Ky,'cubic',0);
% dose = sum(pacbed_rs(:)) *1.0;

% normalize the Objective Aperture to the desired number of electrons
scaling = sqrt(dose/sum(func_ObjApt(:)));
func_ObjApt = func_ObjApt.* scaling;

% convergence aperture function; to filter out the modified  probe func
% during ptychography processing. 
% ObjAptMask = func_ObjApt./sum(func_ObjApt(:));

% probe function - reciprocal space
A = func_ObjApt.*func_transfer;
% probe function - real space
func_probe=fftshift(ifft2(ifftshift(A)));

if ptycho.plot_figures
    figure;
    subplot(1,2,1);
    imagesc(angle(A));axis square;colorbar; axis off
    title('Aperture Phase Surface');
    subplot(1,2,2);
    imagesc(abs(func_probe));axis image; colorbar; axis off
    title('Probe Function');
    
end

ptycho.func_aber = func_aber;
ptycho.func_transfer = func_transfer;
ptycho.func_ObjApt = func_ObjApt;
ptycho.A = A;
ptycho.func_probe = func_probe;
ptycho.aberr_input = aberr_input;
ptycho.scaling = scaling;

ptycho.varfunctions.define_probe_function = 1;

end
