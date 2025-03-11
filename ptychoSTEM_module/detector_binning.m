function [ptycho] = detector_binning(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to adapt matrix m in case the detector was binned.
% Function requires structured array ptycho.
% Input variables required are:
% ptycho.m: 4D array that contains the ronchigram for each probe position
% Output variables modifed to ptycho are:
% ptycho.m: m-matrix, in which the ronchigrams have been binned
% ptycho.varfunctions.detector_binning: Flag to indicate the function has
% been executed
% Detector Binning to reduce the number of pixels in ronchigrams;

try 
    if ptycho.varfunctions.detector_binning >0
        
        if size(ptycho.m,2)>=size(ptycho.m,1)
            binningy = ptycho.varfunctions.detector_binning;
            binningx = round(size(ptycho.m,2)./size(ptycho.m,1))*binningy;
        else
            binningx = ptycho.varfunctions.detector_binning;
            binningy = round(size(ptycho.m,1)./size(ptycho.m,2))*binningx;
        end
        
        by = 1:binningy:size(ptycho.m,1);
        bx = 1:binningx:size(ptycho.m,2);
        
        display(['detector binning: [',num2str(size(ptycho.m,1)),'x',num2str(size(ptycho.m,2)), ...
            '] -> [', num2str(length(by)),'x',num2str(length(bx)), ']'])
        
        mbin1 = zeros(length(by),size(ptycho.m,2),size(ptycho.m,3),size(ptycho.m,4));
        mbin2 = zeros(length(by),length(bx),size(ptycho.m,3),size(ptycho.m,4));
        
        for i = 1:(length(by)-1)
            mbin1(i,:,:,:) = sum(ptycho.m(by(i):(by(i+1)-1),:,:,:),1);
        end
        mbin1(length(by),:,:,:) = sum(ptycho.m(by(end):size(ptycho.m,1),:,:,:),1);
        
        for i = 1:(length(bx)-1)
            mbin2(:,i,:,:) = sum(mbin1(:,bx(i):(bx(i+1)-1),:,:),2);
        end
        mbin2(:,length(bx),:,:) = sum(mbin1(:,bx(end):size(ptycho.m,2),:,:),2);
        
        ptycho.m = mbin2;
        
        % to make sure this function is only executed once, change this flag
        % back to zerol
        ptycho.varfunctions.detector_binning = 0;
        
    end
    
catch
    % in case ptycho.varfunctions.detector_binning is not specified
    display('specified variable ptycho.varfunctions.detector_binning')
end

end




