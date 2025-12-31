%function prepare arrays of doubles into a format that matlab outputs nicely as 16
%bit tif
%
% input parameters: 
% matrix of doubles
%
function [outimage] = prep_matrix_for_imwrite(matrix)

outimage=matrix-min(min(matrix));
outimage=outimage/max(max(outimage));
outimage=uint16(outimage*65535);

end
