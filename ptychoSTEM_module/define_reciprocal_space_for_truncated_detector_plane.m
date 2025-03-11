function [ptycho] = define_reciprocal_space_for_truncated_detector_plane(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To define the reciprocal space at the detector plane for the truncated
% ronchigrams that contain the BF disk only.
% Function requires structured array ptycho. 
% Input variables required are:
% ptycho.theta_x: vector of scattering angle in x
% ptycho.theta_y: vector of scattering angle in y
% ptycho.theta: matrix of scattering angles
% ptycho.det_range_x: vector with the detector range in x for truncated
% ronchigrams
% ptycho.det_range_y: vector with the detector range in y for truncated
% ronchigrams
% Output variables added to ptycho are:
% ptycho.theta_wp: matrix of scattering angles for the truncated
% ronchigrams
% ptycho.phi_wp: matrix of?
% ptycho.Kwp: matrix of reciprocal space vectors for the truncated
% ronchigrams
% ptycho.kx_wp: reciprocal space vector in x for truncated ronchigram
% ptycho.ky_wp: reciprocal space vector in y for truncated ronchigram
% ptycho.varfunctions.define_reciprocal_space_for_truncated_detector_plane: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the reciprocal space for the truncated detector plane.
kx_wp = ptycho.theta_x(ptycho.det_range_x);
ky_wp = ptycho.theta_y(ptycho.det_range_y);
theta_wp   = ptycho.theta(ptycho.det_range_y,ptycho.det_range_x);
[Kx_wp,Ky_wp] = meshgrid(kx_wp,ky_wp);
Kwp = sqrt(Kx_wp.^2 + Ky_wp.^2);
phi_wp = acos(Kx_wp./Kwp);

% singularity issue when theta=0;
phi_wp(Kwp==0)=0;
% when the angle is over 180 degrees,
% phi(theta_y>0) = 2*pi-phi(theta_y>0);
for i=1:size(phi_wp,1)
    for j=1:size(phi_wp,2)
        if Ky_wp(i,j)>0
            phi_wp(i,j) = 2*pi - phi_wp(i,j);
        end
    end
end

ptycho.theta_wp = theta_wp;
ptycho.phi_wp = phi_wp;
ptycho.Kwp = Kwp;
ptycho.kx_wp = kx_wp;
ptycho.ky_wp = ky_wp;
ptycho.Kx_wp = Kx_wp;
ptycho.Ky_wp = Ky_wp;
ptycho.varfunctions.define_reciprocal_space_for_truncated_detector_plane = 1;

end
