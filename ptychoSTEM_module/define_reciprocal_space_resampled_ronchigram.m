function [ptycho] = define_reciprocal_space_resampled_ronchigram(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% this function needs comments

% ptycho.varfunctions.define_reciprocal_space_resampled_ronchigram: Flag to indicate the function has been run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define angle space of the resampled Ronchigram:  kx, ky

% pixel size in Ronchigram is defined by real space FOV (scan area + padding)
kx_delta = ptycho.wavelength/ptycho.pix_size(2)/ptycho.ObjSize(2);
ky_delta = ptycho.wavelength/ptycho.pix_size(1)/ptycho.ObjSize(1);

% Ronchigram space meshgrid: theta_x_rs, theta_y_rs
if mod(ptycho.ObjSize(2),2) == 0
    kx = linspace(-fix(ptycho.ObjSize(2)/2),fix(ptycho.ObjSize(2)/2)-1,ptycho.ObjSize(2)).*kx_delta;
else
    kx = linspace(-fix(ptycho.ObjSize(2)/2),fix(ptycho.ObjSize(2)/2),ptycho.ObjSize(2)).*kx_delta;
end
if mod(ptycho.ObjSize(1),2) == 0
    ky = linspace(-fix(ptycho.ObjSize(1)/2),fix(ptycho.ObjSize(1)/2)-1,ptycho.ObjSize(1)).*ky_delta;
else
    ky = linspace(-fix(ptycho.ObjSize(1)/2),fix(ptycho.ObjSize(1)/2),ptycho.ObjSize(1)).*ky_delta;
end
[Kx, Ky] = meshgrid(kx,ky);
% K = sqrt(Kx.^2 + Ky.^2);
% 
% resmaple the Ronchigram intensities.
Iju_rs = interp2(ptycho.Kx_wp,ptycho.Ky_wp,ptycho.mean_m_wp,Kx,Ky,'cubic',0);
% rotate the Ronchigrams by the known amount(diagnosis from trotters)
Iju_rs = imrotate(Iju_rs, ptycho.rot_angle /pi * 180 ,'bicubic','crop');
% threshold the intensity of Iju_rs
pacbed2 = Iju_rs; 
pacbed2(pacbed2< ptycho.pacbedthreshold *max(pacbed2(:)))=0;
% binary thresholded pacbed
%BW = im2bw(pacbed2);
BW = pacbed2>0;

% fit ellipsis to the re-sampled roated CBED Iju_rs to calibrate the reciprocal space.
%Fit Ellipsis using the matlab built-in function regionprops
% s is a structural file that contains the properties of the ellipsis:
% - Centroid, MajorAxisLength ,MinorAxisLength, Orientation
s = fitEllipsis(BW,'regionprops');
if length(s)>1
    display('warning: more than one ellipsis in the pacbed')
    [~,index] = max([s.MajorAxisLength]);
    s = s(index);
end

kx = ((1:ptycho.ObjSize(2))-s.Centroid(1)) ./ ((s.MajorAxisLength+s.MinorAxisLength)/2) .* ptycho.ObjApt_angle*2;
ky = ((1:ptycho.ObjSize(1))-s.Centroid(2)) ./ ((s.MajorAxisLength+s.MinorAxisLength)/2) .* ptycho.ObjApt_angle*2;
[Kx, Ky] = meshgrid(kx,ky);
K = sqrt(Kx.^2 + Ky.^2);

ptycho.kx = kx;
ptycho.ky = ky;
ptycho.K = K;
ptycho.s_resampled = s;
ptycho.varfunctions.define_reciprocal_space_resampled_ronchigram = 1;
end
