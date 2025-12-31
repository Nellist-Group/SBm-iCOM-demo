function [ ptycho ] = load_data(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To load the ptychography data depending on the format.
% Function requires structured array ptycho.
% Input variables required are:
% ptycho.dataformat: indicates the type of data that will be loaded.
%   Options are:
%           'PN' : for data acquired using PNdetector
%           'Sim': for data that was simulated using software
%           'tif': for data in the tif format 
%           'hdf5_hyperspy': for data using HDF5 format as in hyperspy
%           'hdf5_glasgow': for data using HDF5 format as used in glasgow
%           'EMD': for data in the Electron Microscopy Datasets style hdf5 
%           'mib': for data using the mib format (QD's Merlin)
% Output variables added to ptycho are:
% ptycho.probe_range: vector containing the coordinates of the probe scan
% ptycho.m: 4D array that contains the ronchigram for each probe position
% ptycho.varfunctions.load_data: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ptycho.varfunctions.load_data
    
    display(['matrix m exists, loaded from file type ',ptycho.dataformat])
    display(['file identifier ',ptycho.identifier]);
else
    
    switch ptycho.dataformat
        case 'PN'
            %function to load PN data
            
            % generate/load the Ronchigram Roadmap
            RonchiRoadmap = RonchiRoadmapFunc_v1();
            
            % Select which image to analyze
            % total number of serial acquisitions in the datasets;
            NumImg = size(RonchiRoadmap.SurveyImage,3);
            figure;
            for i=1:size(RonchiRoadmap.SurveyImage,3)
                subplot(1,size(RonchiRoadmap.SurveyImage,3),i)
                imagesc(RonchiRoadmap.SurveyImage(:,:,i));axis image;
            end
            colormap jet
            if NumImg==1
                SelImFrame = 1;
            else
                SelImFrame = input([num2str(NumImg),'images exist, Which one to work on?']);
                if ~sum(ismember(SelImFrame,1:max(RonchiRoadmap.Flag(:))))
                    error('Error: Invalid Input.')
                end
            end
            close;
            
            % Select the Probe Range from the IBF Survey Image
            %  calibrate the real space based on given sampling rate
            %  Load from raw-data the Ronchi frames within the specified probe range
            %if isempty(dir('params.mat'))
            %    probe_range=RangeFunc(RonchiRoadmap.SurveyImage(:,:,SelImFrame));
            %    ptycho.probe_range = probe_range;
            %  save('params.mat','input_params');
            
            %else
            %   load params.mat;
            % if no scanning area (probe_range) saved in previous run
            if  ~isfield(ptycho,'probe_range')
                ptycho.probe_range=RangeFunc(RonchiRoadmap.SurveyImage(:,:,SelImFrame));
            else
                %comment all this to do automated upload
                % if probe_range exists, show the scanning area and ask if keep using it.
                %                 figure;
                %                 imagesc(RonchiRoadmap.SurveyImage(:,:,SelImFrame));
                %                 axis image; hold on;
                %                 rectangle('Position',ptycho.probe_range,'LineWidth',2,'LineStyle','--');
                %                 hold off;
                %                 IfChangeProbeRange = input('change probe area? [y]/[n]: ', 's');
                %                 if IfChangeProbeRange == 'y'
                %                     ptycho.probe_range=RangeFunc(RonchiRoadmap.SurveyImage(:,:,SelImFrame));
                %                 end
                %                 close;
            end
            %end
            
            ptycho.m = RoadmapScanner(RonchiRoadmap,SelImFrame,ptycho.probe_range);
            % the four dimensions of matrix m are: [ronchi-y,ronchi-x, scan-y, scan-x]
            
        case 'Sim'
            
            %function to load simulated data
            if isfield(ptycho,'data_filename')
            	data_as_struct=load(ptycho.data_filename);
                fns=fieldnames(data_as_struct);
                m=data_as_struct.(fns{1});
                clear data_as_struct;
	elseif isfield(ptycho,'simulation_name')
		load(ptycho.simulation_name)
		warning('Warning: ptycho.simulation_name is deprecated, use ptycho.data_filename in the future');
	else
		error('error: No filename');
	%	return
            end
    display('disp');
            %   creating m matrix from simulated data
            for ii=1:1:size(m,3)
                for jj=1:1:size(m,4)
                    m(:,:,ii,jj)=squeeze(m(:,:,ii,jj))';
                    %m(:,:,ii,jj)=squeeze(m(:,:,ii,jj));
                end
            end
            
            %to repeat unit cell of simulations
            if ptycho.repeat_unit_cell_simulation
                maux = m;
                for countery = 1:1:ptycho.repeat_number_unit_cell_in_x-1
                    m = cat(3,m,maux);
                end
                maux = m;
                for counterx = 1:1:ptycho.repeat_number_unit_cell_in_y-1
                    m = cat(4,m,maux);
                end
            end
            
            ptycho.probe_range = [0 0 size(m,4) size(m,3)]; %whole simulation
            ptycho.m = m;
            
        case 'mib'
            %note that ptycho.data_filename should be path/filename only up 
            %to before the 1.mib not the full filename
            ptycho.m=read_mib(ptycho.data_filename,ptycho.numprobeposx,...
            ptycho.numprobeposy,ptycho.num_extra_images_per_row,...
            ptycho.imagesperfile,ptycho);
            ptycho.probe_range = [0 0 size(ptycho.m,4) size(ptycho.m,3)]; 
              
        case 'tif'
            
            %function to load tif data. Should be made more generic.
            
%            [FileTif, pathname] = uigetfile({'*.tif','Select the Ronchigram Series File'}) ;
            %FileTif='~/data/Andy/tuning_testing/75/andy_75_128pix.tif';
            FileTif=ptycho.data_filename;

            InfoImage=imfinfo(FileTif);
            %InfoImage=imfinfo(FileTif);
            ronchigramSizex=InfoImage(1).Width;
            ronchigramSizey=InfoImage(1).Height;
            numberRonchigrams=length(InfoImage);
            
            %number of probe positions x and y including probe fly-back
            numpixelx0 = 34%134 %66%258;
            numpixely0 = 32%128 %64%256;
            %number of "actual" probe positions x and y
            numpixelx = 32% 64%256;
            numpixely = 32 %64%256;
            % scanwindow = [probe_range(4),probe_range(3)];  %[yy, xx]
            
            % allocate memory for m
            %     m = zeros(ronchigramSizey,ronchigramSizex,numpixelx,numpixely);
            for i=1:numberRonchigrams
                i
                yy = fix((i-1)/numpixelx0)+1;
                xx = i - (yy-1)*numpixelx0;
                % if probe position (yy,xx) is within the "actual" range
                if yy<=numpixely && xx <= numpixelx
                    m(:,:,yy,xx) = imread(FileTif,'Index',i,'Info',InfoImage);
                end
            end
            
            m = im2double(m);
            
            ptycho.probe_range = [0 0 size(m,3) size(m,4)]; %whole simulation
            ptycho.m = m;
            
                        
        case 'hdf5_glasgow'
            %function to load Glasgow style HDF5-file
            [fileHDF5, pathname] = uigetfile({'*.hdf5','Select HDF5 file'}) ;
            filedata = struct('data',h5read(fileHDF5, '/fpd_expt/fpd_data/data'));
            
            m = double(squeeze(filedata.data));
            
            ptycho.probe_range = [0 0 size(m,4) size(m,3)]; %whole simulation
            ptycho.m = m;

        case 'hdf5_hyperspy'
            %function to load HyperSpy style HDF5-file
            [fileHDF5, pathname] = uigetfile({'*.hdf5','Select HDF5 file'});
            %filedata = struct('data',h5read(fileHDF5,'/Experiments/__unnamed__/data'));
            filedata = struct('data',h5read(fileHDF5,'/Experiments/-fpd_expt-fpd_data/data'));
            m = double(squeeze(filedata.data));
            %ptycho.probe_range = [0 0 size(m,3) size(m,4)];
            ptycho.probe_range = [0 0 size(m,4) size(m,3)];
            ptycho.m = m;

       case 'EMD'
            [fileHDF5, pathname] = uigetfile({'*.emd','Select EMD file'});
            s = load_emd_file(fileHDF5);
            m = double(squeeze(s.signal1.data));
            m = permute(m,[3 4 1 2]);
            ptycho.probe_range = [0 0 size(m,4) size(m,3)];
            ptycho.m = m;
    end
    
end

ptycho.varfunctions.load_data = 1;

end


