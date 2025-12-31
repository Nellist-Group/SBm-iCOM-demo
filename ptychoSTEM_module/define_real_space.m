function [ptycho] = define_real_space(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To define real space (in Angstrom) for each pixel in the image
% Function requires structured array ptycho, with fields obtained from parameter
% file using load_parameters and load_data functions
% Input variables required are:
% ptycho.probe_range: vector containing the coordinates of the probe scan
% ptycho.use_padding: Flag to define is zero-padding is used
% ptycho.pix_size: Image pixel size in Angstrom
% ptycho.wavelength: Wavelength in Angstrom
% Output variables added to ptycho are:
% ptycho.scanwindow: vector with the pixel values of the scanned area 
% ptycho.ObjSize: vector with the pixel of the scanned area + zero-padding (if
% stated in parameter file). Padding is useful to avoid artifacts when
% doing FFTs.
% ptycho.X0: vector of probe position coordinates in x 
% ptycho.Y0: vector of probe position coordinates in y 
% ptycho.Rx: matrix of probe position coordinates in x from meshgrid 
% ptycho.Ry: matrix of probe position coordinates in y from meshgrid 
% ptycho.Q_x: vector of reciprocal space coordinates with respect to probe
% position in x
% ptycho.Q_y: vector of reciprocal space coordinates with respect to probe
% position in y 
% ptycho.Q : matrix of reciprocal space coodinates with respect to probe
% position (Q-space)
% ptycho.varfunctions.define_real_space: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%scanwindow = [ptycho.probe_range(3),ptycho.probe_range(4)] %[yy, xx]
scanwindow = [ptycho.probe_range(4),ptycho.probe_range(3)]; %[yy, xx] needs
%revision

if ptycho.use_padding
ObjSize = scanwindow + fix( 2 ./ptycho.pix_size ); %[yy,xx]
% ObjSize = [ 1, 1 ] .* max(scanwindow + fix( 2 ./ptycho.pix_size )) %[y,x]
else
ObjSize = scanwindow;
end

% if isempty(dir('params.mat'))
%     save('params.mat','ptycho');
% else
%     save('params.mat','ptycho','-append');
% end

% generate a meshgrid for probe positions
if mod(ObjSize(2),2) == 0
    X0 = linspace(-fix(ObjSize(2)/2),fix(ObjSize(2)/2)-1,ObjSize(2)).*ptycho.pix_size(2);
else
    X0 = linspace(-fix(ObjSize(2)/2),fix(ObjSize(2)/2),ObjSize(2)).*ptycho.pix_size(2);
end
if mod(ObjSize(1),2) == 0
    Y0 = linspace(-fix(ObjSize(1)/2),fix(ObjSize(1)/2)-1,ObjSize(1)).*ptycho.pix_size(1);
else
    Y0 = linspace(-fix(ObjSize(1)/2),fix(ObjSize(1)/2),ObjSize(1)).*ptycho.pix_size(1);
end
[Rx,Ry] = meshgrid(X0,Y0);

% define the reciprocal space vector Q with respect to probe position:
Q_x_delta = ptycho.wavelength/ptycho.pix_size(2)/(fix(ObjSize(2)/2)*2);
Q_y_delta = ptycho.wavelength/ptycho.pix_size(1)/(fix(ObjSize(1)/2)*2);


if mod(ObjSize(2),2) == 0
    Q_x = linspace(-fix(ObjSize(2)/2),fix(ObjSize(2)/2)-1,ObjSize(2)).*Q_x_delta;
else
    Q_x = linspace(-fix(ObjSize(2)/2),fix(ObjSize(2)/2),ObjSize(2)).*Q_x_delta;
end
if mod(ObjSize(1),2) == 0
    Q_y = linspace(-fix(ObjSize(1)/2),fix(ObjSize(1)/2)-1,ObjSize(1)).*Q_y_delta;
else
    Q_y = linspace(-fix(ObjSize(1)/2),fix(ObjSize(1)/2),ObjSize(1)).*Q_y_delta;
end
[q1,q2] = meshgrid(Q_x,Q_y);
Q = sqrt(q1.*q1+q2.*q2);


ptycho.scanwindow = scanwindow;
ptycho.ObjSize = ObjSize;
ptycho.X0 = X0;
ptycho.Y0 = Y0;
ptycho.Rx = Rx;
ptycho.Ry = Ry;
ptycho.Q_x = Q_x;
ptycho.Q_y = Q_y;
ptycho.Q = Q;
ptycho.varfunctions.define_real_space = 1;

end
