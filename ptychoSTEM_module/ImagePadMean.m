function padded = ImagePadMean(input_file, padsize)
% this function adds the padding around the scan if requested
% (ptycho.use_padding=1)
%input_file is the array of intensity values across the whole unpadded scan
%padsize is the size of the scan plus the padding

original_size = size(input_file);

% check if the padsize is larger than original size
if original_size(1)>padsize(1) || original_size(2)>padsize(2)
    error('pad size smaller than original image size')
end
%if not, you are probably loading your data in incorrectly somehow
%for instance 
%   ptycho.probe_range = [0 0 size(ptycho.m,3) size(ptycho.m,4)]; 
%rather than
%   ptycho.probe_range = [0 0 size(ptycho.m,4) size(ptycho.m,3)]; 

% make the padded area have the value of the averge intensity of this
% detector pixel over the scan
padded = ones(padsize).* mean(input_file(:));

%find the difference in the size of the padded scan and original
diff_size = padsize - original_size;

%set things inside orignal scan range to the original scan values, leaving
%the padding at the average value around the edges of the padded data.
start = round(diff_size/2);
padded(start(1)+1:original_size(1)+start(1),start(2)+1:original_size(2)+start(2)) = input_file;

end