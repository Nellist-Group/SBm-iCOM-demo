function [ptycho] = define_reciprocal_space(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To define reciprocal space angle for each pixel in the detector plane
% Function requires structured array ptycho.
% Input variables required are:
% ptycho.pacbedthreshold: Treshold value to select the pacbed bright field
% disk area
% ptycho.m: 4D array that contains the ronchigram for each probe position
% Output variables added to ptycho are:
% ptycho.pacbed: matrix containing the Position Averaged Convergent Beam
% Electron Diffraction (PACBED) pattern.
% ptycho.s: structure object that contains the parameters of fitting an
% ellipsis to the pacbed.
% ptycho.theta_x: vector of scattering angle in x
% ptycho.theta_y: vector of scattering angle in y
% ptycho.tx: matrix of scattering angle in x from meshgrid
% ptycho.ty: matrix of scattering angle in y from meshgrid
% ptycho.theta: matrix of scattering angles
% ptycho.varfunctions.define_reciprocal_space: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the reciprocal space angle for each pixel in the detector plane

display('calibrating the detector plane')
% find the centre of the BF disc and its radius
% it works for both circles using unbinned detector, and ellipsis using binned detector.

% calculate pacbed
pacbed = mean(mean(ptycho.m,4),3);
% thresholding the pacbed
pacbed2 = pacbed;  pacbed2(pacbed2<ptycho.pacbedthreshold *max(pacbed2(:)))=0;
% binary thresholded pacbed
%BW = im2bw(pacbed2);
BW = pacbed2>0; %added by Gerardo
%BW = edge(pacbed_thres,'Sobel');

%Fit Ellipsis using the matlab built-in function regionprops
% s is a structural file that contains the properties of the ellipsis:
% - Centroid, MajorAxisLength ,MinorAxisLength, Orientation
%%% note: sometimes the position of BF disc is determined from a separate
%%% experiment (recording the BF of e-beam passing through vacuum), 
%%% then the BF position information(variable s) is saved into a file
%%% called "s_vac.mat".  check if there's such a file before calculating
%%% the s variable.
if ~exist('s_vac.mat','file')
    s = fitEllipsis(BW,'regionprops');
    if length(s)>1
        display('warning: more than one ellipsis in the pacbed')
        [~,index] = max([s.MajorAxisLength]);
        s = s(index);
    end
else
    load s_vac.mat
end

% calculate the elipsis edge shape in coordinate [x,y].
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
xbar = s.Centroid(1);
ybar = s.Centroid(2);
a = s.MajorAxisLength/2;
b = s.MinorAxisLength/2;
thetaElipsis = pi*s.Orientation/180;
R = [ cos(thetaElipsis)   sin(thetaElipsis)
    -sin(thetaElipsis)   cos(thetaElipsis)];
xy = [a*cosphi; b*sinphi];
xy = R*xy;
x = xy(1,:) + xbar;
y = xy(2,:) + ybar;

if ptycho.plot_figures
figure; 
subplot(3,1,1); 
imagesc(pacbed); title('pacbed');
hold on; plot(x,y,'.r'); hold off
axis image
subplot(3,1,2); imagesc(pacbed2);title('thresholded pacbed')
hold on; plot(x,y,'.r'); hold off
axis image
subplot(3,1,3); imagesc(BW); title('binary image of pacbed')
hold on; plot(x,y,'.r'); hold off
axis image
end
% Based on the fitted ellipsis, the detector space angle is defined below using a given probe convergence angle.
% Note: the convergence angle here may not be calibrated !!!! Calibration of the convergence angle is a separate task.
if abs(s.Orientation) < 90-abs(s.Orientation)
    %x-rediction is the major axis
    theta_x = ((1:size(ptycho.m,2))-s.Centroid(1))./s.MajorAxisLength*ptycho.ObjApt_angle*2;
    theta_y = ((1:size(ptycho.m,1))-s.Centroid(2))./s.MinorAxisLength*ptycho.ObjApt_angle*2;
else
    %y-rediction is the major axis
    theta_x = ((1:size(ptycho.m,2))-s.Centroid(1))./s.MinorAxisLength*ptycho.ObjApt_angle*2;
    theta_y = ((1:size(ptycho.m,1))-s.Centroid(2))./s.MajorAxisLength*ptycho.ObjApt_angle*2;
end
[tx,ty] = meshgrid(theta_x, theta_y);
theta = sqrt(tx.^2 + ty.^2);

%%%%%%%%%% Zhiyuan add %%%%%%%%
if isfield(ptycho,'set_dk_rad') && ~isempty(ptycho.set_dk_rad)
    s = [];
    theta_x = ((0:size(ptycho.m,2)-1) - size(ptycho.m,2)/2) .* ptycho.set_dk_rad;
    theta_y = ((0:size(ptycho.m,1)-1) - size(ptycho.m,1)/2) .* ptycho.set_dk_rad;
    [tx,ty] = meshgrid(theta_x, theta_y);
    theta = sqrt(tx.^2 + ty.^2);
end

ptycho.s = s;
ptycho.theta= theta;
ptycho.theta_x= theta_x;
ptycho.theta_y= theta_y;
ptycho.tx = tx;
ptycho.ty = ty;
ptycho.pacbed = pacbed;
ptycho.varfunctions.define_reciprocal_space = 1;


end
