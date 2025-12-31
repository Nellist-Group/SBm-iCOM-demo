%function to create the m matrix needed for ptycho
% input parameters:
% filename: name of the files created by MULTEM in the following format:
%       name_counterX_counterY
% where counterX and counterY corresponds to the probe positions
%
% variable size4D: vector with the settings of the m-matrix
% the values should be arranged as:
% [detectorsize_x detectorsize_y probeposition_x probeposition_y]
%
% crop is how many pixels should be retained in the center - assumes square
%
function [m] = make_4D_matrix_m_qstem(filename,size4D,crop)

if nargin < 3
    crop = 0;
end

m=zeros(crop,crop,size4D(3),size4D(4));

tic

for ii=0:1:size4D(3)-1
    for jj =0:1:size4D(4)-1

        loadname = strcat(filename,'_',num2str(ii),'_',num2str(jj),'.img')

        data=binread2D(loadname);
        data=data(1+size4D(1)/2 - crop/2:size4D(1)/2 + crop/2,1+ size4D(2)/2 - crop/2:size4D(2)/2 + crop/2);
        m(:,:,ii+1,jj+1)=data;

        [ii+1 jj+1]

    end

    clear data
end

display(strcat(filename,' 4D matrix has been created'));

toc

end

