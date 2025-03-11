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
function [m] = make_4D_matrix_m(filename,size4D)


m=zeros(size4D);

tic

for ii=1:1:size4D(3)
    for jj =1:1:size4D(4)

        loadname = strcat(filename,'_',num2str(ii),'_',num2str(jj),'.mat');

        data=load(loadname,'crys_thick');

        m(:,:,ii,jj)=data.crys_thick(1).cbed(:,:);

        [ii jj]

    end

    clear data
end

display(strcat(filename,' 4D matrix has been created'));

toc

end

