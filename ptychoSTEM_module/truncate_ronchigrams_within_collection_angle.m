function [ptycho] = truncate_ronchigrams_within_collection_angle(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To select only the bright field disk from the ronchigram that is going to
% be used in the ptychographyc reconstruction. This function will define a
% detector range to truncate the ronchigram to include the bright field 
% disk only.
% Function requires structured array ptycho. 
% Input variables required are:
% ptycho.m: 4D array that contains the ronchigram for each probe position
% ptycho.maskWP: matrix of mask for Weak Phase Ptychography
% If ptycho.maskWP is not available, it will calculate it. This will
% require the following variables:
%   ptycho.ObjApt_angle: Probe convergence angle
%   ptycho.theta: matrix of scattering angles
% Output variables added to ptycho are:
% ptycho.m_wp: m-matrix with truncated ronchigrams
% ptycho.det_range_x: vector with the detector range in x for truncated
% ronchigrams
% ptycho.det_range_y: vector with the detector range in y for truncated
% ronchigrams
% ptycho.varfunctions.truncate_ronchigrams_within_collection_angle: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MaxCollectionAngle: maximum scattring angle used in ptychography
try 
    % if ptycho.MaxCollectionAngle is specified in the paramfile
    MaxCollectionAngle = ptycho.MaxCollectionAngle;
catch
    % if ptycho.MaxCollectionAngle is not specified in the paramfile, set
    % default values as 1.1 times of the convergence semi-angle.
    MaxCollectionAngle = ptycho.ObjApt_angle*1.1;
end

if isfield(ptycho,'maskWP')
    maskWP = ptycho.maskWP;
else
    maskWP = ones(size(ptycho.theta));
    maskWP(ptycho.theta>MaxCollectionAngle)=0;
end

% Truncate Ronchigrams to be within the specified collection angle 

% ---------------------------------
% for single side band, use BF only
% for PIE and ePIE,     use BF+DF 
% ---------------------------------

% seems to be a potential bug here when BFdisk very close to edge 
% of detector already
% for weak phase
% % This just looks for where on the mask the the row or column sum
% % go to zero to find left, right top bottom of the BF disk
% % So if it is already cropped, can just set 
% % ptycho.det_range_x=1:size(ptycho.m,1) or is it ,2)
% % ptycho.det_range_y=1:size(ptycho.m,2) or is it ,1)
% % ptycho.m_wp=ptycho.m
% % and skip this function entirely
% % 

maskWP_x = sum(maskWP,1);  maskWP_y = sum(maskWP,2);
%%%%%%%% changed by Zhiyuan %%%%%%%%
if isfield(ptycho,'mask_mode') && strcmp(ptycho.mask_mode,'none')
    ptycho.det_range_wp_x = [1,length(maskWP_x)];
    ptycho.det_range_wp_y = [1,length(maskWP_y)];
else % default: no field ptycho.mask_mode, or mask_mode = 'semiangle'
    for i=1:length(maskWP_x)-1
        if maskWP_x(i)==0 && maskWP_x(i+1)>0
            ptycho.det_range_wp_x(1)=i+1;
        elseif maskWP_x(i)>0 && maskWP_x(i+1)==0
            ptycho.det_range_wp_x(2)=i;
        end
    end
    for i=1:length(maskWP_y)-1
        if maskWP_y(i)==0 && maskWP_y(i+1)>0
            ptycho.det_range_wp_y(1)=i+1;
        elseif maskWP_y(i)>0 && maskWP_y(i+1)==0
            ptycho.det_range_wp_y(2)=i;
        end
    end
end 
%%%%%%%% changed by Zhiyuan END %%%%%%%%%%%

if ptycho.det_range_wp_x(1)== 0
    ptycho.det_range_wp_x(1) = 1;
end
if ptycho.det_range_wp_y(1)== 0
    ptycho.det_range_wp_y(1) = 1;
end

% mask and crop the M matrix
det_range_x = ptycho.det_range_wp_x(1):ptycho.det_range_wp_x(2);
det_range_y = ptycho.det_range_wp_y(1):ptycho.det_range_wp_y(2);

%%%%%%%%%%% Changed by Zhiyuan %%%%%%%%%%%%
if isfield(ptycho,'mask_mode') && strcmp(ptycho.mask_mode,'none')
    m_wp = ptycho.m;
else % default: no field ptycho.mask_mode, or mask_mode = 'semiangle'
    m_wp = zeros(length(det_range_y),length(det_range_x),size(ptycho.m,3),size(ptycho.m,4));
    for yy=1:size(ptycho.m,3)
        for xx=1:size(ptycho.m,4)
            temp = ptycho.m(:,:,yy,xx).*maskWP;
            m_wp(:,:,yy,xx) = temp(det_range_y,det_range_x);
        end
    end
end
%%%%%%%%%% Changed by Zhiyuan END %%%%%%%%%%
temp = mean(mean(m_wp,4),3);

if ptycho.plot_figures
    figure; imagesc(log(temp-min(temp(:))+1));axis image
end

% if ptycho.remove_m_to_save_space
%     ptycho.m = 'removed to save space';
%     ptycho.varfunctions.remove_m_to_save_space = 1;
% end;

ptycho.m_wp = m_wp;
ptycho.det_range_x = det_range_x;
ptycho.det_range_y = det_range_y;
ptycho.varfunctions.truncate_ronchigrams_within_collection_angle = 1;
end
