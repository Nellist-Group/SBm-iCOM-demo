function [ptycho] =  Aberration_corrected_single_side_band_reconstruction(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To perform the single side band reconstruction. 
% It looks at the Kf space (detector plane) at every Q value, and extract
% the phase and amplitude information.
% The algorithm focuses on the overlapping discs between scattered beam and
% central bright discs, therefore frequencies that do not have overlapping
% regions with bright field discs are neglected. Using the aberrations
% determined either automatically using SVD or manually, the wave aberration function is adjusted to
% correct for aberrations.

% Function requires structured array ptycho. Figures are displayed if
% ptycho.plot_figures is set to 1.

% Input variables required are:
% ptycho.scanwindow: vector with the pixel values of the scanned area 
% ptycho.ObjSize: vector with the pixel of the scanned area + zero-padding
% ptycho.ObjApt_angle: Probe convergence angle
% ptycho.rot_angle: Rotation angle 
% ptycho.X0: vector of probe position coordinates in x 
% ptycho.Y0: vector of probe position coordinates in y 
% ptycho.kx_wp: reciprocal space vector in x for truncated ronchigram
% ptycho.ky_wp: reciprocal space vector in y for truncated ronchigram
% ptycho.Q_x: vector of reciprocal space coordinates with respect to probe
% position in x
% ptycho.Q_y: vector of reciprocal space coordinates with respect to probe
% position in y 
% ptycho.Q : matrix of reciprocal space coodinates with respect to probe
% position (Q-space)
% ptycho.G_wp: 4D array that contains the ronchigram for each spatial frecuency of probe position  
% ptycho.IBFimg: matrix of Incoherent Bright Field image  
% ptycho.ABFimg: matrix of Annular Bright Field image
% ptycho.DFimg: matrix of Annular Dark Field image

% Output variables added to ptycho are:
% ptycho.trotterimgL = trotterimgL;-left SSB phase before aberration correction
% ptycho.trotterimgR = trotterimgR;-right SSB phase before aberration correction
% ptycho.trotterimgLC = trotterimgLC;-left SSB phase after aberration correction
% ptycho.trotterimgRC = trotterimgRC;-right SSB phaseafter aberration correction
% ptycho.varfunctions.single_side_band_reconstruction: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc
disp('running single side band (SSB) reconstruction')

trotterphasesumL  = zeros(ptycho.scanwindow);
trotterphasesumR  = zeros(ptycho.scanwindow);
trotterphasesumLC  = zeros(ptycho.scanwindow);
trotterphasesumRC  = zeros(ptycho.scanwindow);



for yy=1:ptycho.ObjSize(1)
%         disp(['SSB: ',num2str(yy/ptycho.ObjSize(1)*100),'%'])
    for xx=1:ptycho.ObjSize(2) 

%         Define coordinates with respect to rotation angle, for applying
%         masks to G(K,Q)
        
        Q_x_rot = ptycho.Q_x(xx).*cos(ptycho.rot_angle) - ptycho.Q_y(yy).*sin(ptycho.rot_angle);
        Q_y_rot = ptycho.Q_x(xx).*sin(ptycho.rot_angle) + ptycho.Q_y(yy).*cos(ptycho.rot_angle);
        [dx,dy] = meshgrid(ptycho.kx_wp + Q_x_rot , ptycho.ky_wp + Q_y_rot);
        d1 = sqrt(dx.*dx+dy.*dy);
        [dx,dy] = meshgrid(ptycho.kx_wp - Q_x_rot , ptycho.ky_wp - Q_y_rot);
        d2 = sqrt(dx.*dx+dy.*dy);
        [dx,dy] = meshgrid(ptycho.kx_wp, ptycho.ky_wp);
        d3 = sqrt(dx.*dx+dy.*dy);
        
        g =  squeeze(ptycho.G_wp(:,:,yy,xx));
        
        Kx_wp_plusQ  = ptycho.Kx_wp + Q_x_rot;
        Ky_wp_plusQ  = ptycho.Ky_wp + Q_y_rot;

          
        Kx_wp_minusQ = ptycho.Kx_wp - Q_x_rot;
        Ky_wp_minusQ = ptycho.Ky_wp - Q_y_rot;

        
        
%         initially diagnose aberrations from initial value of G(K,Q) before
%         applying focus steps
        
        func_aber = FuncAberrUV(ptycho.Kx_wp,ptycho.Ky_wp,ptycho.aberr_input);
        func_aber_plusQ = FuncAberrUV(Kx_wp_plusQ,Ky_wp_plusQ,ptycho.aberr_input);
        func_aber_minusQ = FuncAberrUV(Kx_wp_minusQ,Ky_wp_minusQ,ptycho.aberr_input);
        
%         calculate wave aberration function from determined aberrations
        
        func_transfer = exp(-sqrt(-1)*2*pi/ (ptycho.wavelength*1e-10) .* func_aber);
        func_transfer_plusQ = exp(-sqrt(-1)*2*pi/ (ptycho.wavelength*1e-10) .* func_aber_plusQ);
        func_transfer_minusQ = exp(-sqrt(-1)*2*pi/ (ptycho.wavelength*1e-10) .* func_aber_minusQ);
        
%         define mask for trotter amplitude:
        
        g_ampL = abs(g);
        g_ampL(d1>ptycho.ObjApt_angle)=0;
        g_ampL(d2<ptycho.ObjApt_angle)=0;
        g_ampL(d3>ptycho.ObjApt_angle)=0;
        g_ampR = abs(g);
        g_ampR(d1<ptycho.ObjApt_angle)=0;
        g_ampR(d2>ptycho.ObjApt_angle)=0;
        g_ampR(d3>ptycho.ObjApt_angle)=0;
        if ptycho.Q(yy,xx) == 0
            g_ampL = abs(g);
            g_ampR = abs(g);
        end
        
%         define mask for trotter phase

        g_phaseL = angle(g);
        g_phaseL(d1>ptycho.ObjApt_angle)=0;
        g_phaseL(d2<ptycho.ObjApt_angle)=0;
        g_phaseL(d3>ptycho.ObjApt_angle)=0;
        g_phaseR = angle(g);
        g_phaseR(d1<ptycho.ObjApt_angle)=0;
        g_phaseR(d2>ptycho.ObjApt_angle)=0;
        g_phaseR(d3>ptycho.ObjApt_angle)=0;

        
%        Left trotter before aberration correction
        trotterL = g_ampL.*exp(sqrt(-1).*g_phaseL);
%        Left trotter after aberration correction
        trotterLC = g_ampL.*exp(sqrt(-1).*g_phaseL).*(func_transfer.*conj(func_transfer_plusQ));
       
%       Sum the phase of each left trotter
        trotterphasesumL(yy,xx) = sum(trotterL(:)); % uncorrected
        trotterphasesumLC(yy,xx) = sum(trotterLC(:)); % corrected

%        Right trotter before aberration correction
        trotterR = g_ampR.*exp(sqrt(-1).*g_phaseR);
%        Right trotter after aberration correction
        trotterRC = g_ampR.*exp(sqrt(-1).*g_phaseR).*(func_transfer_minusQ.*conj(func_transfer));

%       Sum the phase of each right trotter
        trotterphasesumR(yy,xx) = sum(trotterR(:));
        trotterphasesumRC(yy,xx) = sum(trotterRC(:));

%         if xx == 31 && yy == 25
%             a = 1;
%         end
    end
%     clc
end

% Phase images before aberration correction, for left and right trotters
trotterimgL = fftshift(ifft2(ifftshift(trotterphasesumL)));
trotterimgR = fftshift(ifft2(ifftshift(trotterphasesumR)));

% Phase images after aberration correction, for left and right trotters
trotterimgLC = fftshift(ifft2(ifftshift(trotterphasesumLC)));
trotterimgRC = fftshift(ifft2(ifftshift(trotterphasesumRC)));

ptycho.Q_x_rot = Q_x_rot;
ptycho.Q_y_rot = Q_y_rot;


ptycho.trotterimgL = trotterimgL;
ptycho.trotterimgR = trotterimgR;
ptycho.trotterimgLC = trotterimgLC;
ptycho.trotterimgRC = trotterimgRC;
ptycho.varfunctions.single_side_band_reconstruction = 1;
ptycho.func_transfer_minusQ = func_transfer_minusQ;
ptycho.func_transfer = func_transfer;

if ptycho.plot_figures
    show_SSB_results(ptycho);
end

end