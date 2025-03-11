function [ptycho] = ronchigram_resampling(ptycho)
% this function needs comments


% Ronchigram resampling
if exist('RonchiResampling','var') && RonchiResampling==1
   
    SizeRonchiRS = round(sqrt( length(det_range_x)*length(det_range_y)  ));
    size4D(1) = SizeRonchiRS;
    size4D(2) = SizeRonchiRS;
    input_params.size4D = size4D;
    if isempty(dir('params.mat'))
        save('params.mat','input_params');
    else
        save('params.mat','input_params','-append');
    end
    
    if abs(s.Orientation) < abs(s.Orientation-90)
        %x-rediction is the major axis
        theta_x = ((1:size(m,2))-s.Centroid(1))./s.MajorAxisLength*input_params.ObjApt_angle*2;
        theta_y = ((1:size(m,1))-s.Centroid(2))./s.MinorAxisLength*input_params.ObjApt_angle*2;
    else
        %y-rediction is the major axis
        theta_x = ((1:size(m,2))-s.Centroid(1))./s.MinorAxisLength*input_params.ObjApt_angle*2;
        theta_y = ((1:size(m,1))-s.Centroid(2))./s.MajorAxisLength*input_params.ObjApt_angle*2;
    end
    [tx,ty] = meshgrid(theta_x, theta_y);
    theta = sqrt(tx.^2 + ty.^2);

        
    % redefine the real and reciprocal space angles
    kx_wp = theta_x(det_range_x);
    ky_wp = theta_y(det_range_y);
    theta_wp   = theta(det_range_y,det_range_x);
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

    % resample the 4D dataset m_wp, and calculate G again.
    display('FT{M} w.r.t Probe Positions...')
    % allocate memory for G matrix
    G_wp = zeros(size4D);
    for i = 1:size4D(1)
        for j=1:size4D(2)
            disp(['FT{M}:',num2str(i), '/', num2str(j)]);
            G_wp(i,j,:,:) = fftshift(fft2(ifftshift( ImagePadMean(squeeze(m_wp(i,j,:,:)),ObjSize) )));
        end
    end
    
end


end