function m=read_mib(data_basefilename,numprobeposx,numprobeposy,...
    num_extra_images_per_row,imagesperfile,ptycho)

if nargin<5 %check we have a minium number of variables
    error('not enough input variables\n');
end

%%%%%% Check if a subsection of probe positions is specified
%%%%%% if not, read in the whole scan
if isfield(ptycho,'xmin_probepos')
    xmin_probepos=ptycho.xmin_probepos
else
    xmin_probepos=1;
end
if isfield(ptycho,'xmax_probepos')
    xmax_probepos=ptycho.xmax_probepos
else
    xmax_probepos=numprobeposx;
end
if isfield(ptycho,'ymin_probepos')
    ymin_probepos=ptycho.ymin_probepos
else
    ymin_probepos=1;
end
if isfield(ptycho,'ymax_probepos')
    ymax_probepos=ptycho.ymax_probepos
else
    ymax_probepos=numprobeposy;
end
%%%%%%%%%%%%

%%read the header in the mib (assumes there is one...)
%open binary ('b') file to read ('r')
file=fopen(strcat(data_basefilename,'1','.mib'),'r','b');
%get header
header=char(fread(file,50));
%fclose(file);
%find the information and turn into values
celled_header=textscan(header, '%s %d %d %d %d %d  %s %s','delimiter', ',');
%get relavent values
data_offset=cell2mat(celled_header(3)); %length of header (offset to data)
data_offset=data_offset(1);
num_chips=cell2mat(celled_header(4)); %number of medipix chips
xsize=cell2mat(celled_header(5)); %image size, num columns
ysize=cell2mat(celled_header(6)); %image size, num rows
data_type=celled_header{7};

if strcmp(data_type,'U32')
    mlab_data_type='uint32'
    bytes=xsize*ysize*4;
elseif strcmp(data_type,'U16')
    mlab_data_type='uint16'
    imagesize_bytes=xsize*ysize*2;
elseif strcmp(data_type,'U08')
    mlab_data_type='uint8'
    imagesize_bytes=xsize*ysize*1;
else
    error('Error: data type not recognized - raw (R64) not yet added');
end

%number of "actual" images per row, including any extra throw away images
% for syncing
numimagesperrowraw= numprobeposx+num_extra_images_per_row;
%total number of images
numberRonchigrams=numprobeposy*numimagesperrowraw;

%what about the headers when loading the data?

% output datasize
yscanlengthout=ymax_probepos-ymin_probepos;
xscanlengthout=xmax_probepos-xmin_probepos;

% allocate space before loading for decent performance
m = zeros(ysize,xsize,yscanlengthout,xscanlengthout);

%initialize imgfilenum to zero so auto incrementation works
imagefilenum=0;

%%%%%%% now read the data %%%%%%%%
%file=fopen(strcat(data_basefilename,num2str(imagefilenum),'.mib'),'r','b');
for i=1:numberRonchigrams
    i
    if i==fix(i/imagesperfile)*imagesperfile +1
        disp("inif:");
        i
        fclose(file);
        imagefilenum=imagefilenum+1
        file=fopen(strcat(data_basefilename,num2str(imagefilenum),'.mib'),'r','b');
    end
    yy = fix((i-1)/numimagesperrowraw)+1;
    xx = i - (yy-1)*numimagesperrowraw;
     %header=char(fread(file,data_offset));
     fseek(file,data_offset,'cof');
    %if probe position (yy,xx) is within region of interest, load data
    if ymin_probepos<=yy && yy<=ymax_probepos && xmin_probepos<=xx && xx<=xmax_probepos
        image=fread(file,[xsize ysize],mlab_data_type);
        m(:,:,yy-ymin_probepos+1,xx-xmin_probepos+1) = image;
    else %if not just move pointer along
        fseek(file,imagesize_bytes,'cof');
    end
end


end
