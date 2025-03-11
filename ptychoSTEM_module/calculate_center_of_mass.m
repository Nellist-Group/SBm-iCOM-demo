function [ptycho] = calculate_center_of_mass(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% To calculate and show Center of Mass of ronchigrams. 
% Function requires structured array ptycho. Figures are displayed if
% ptycho.plot_figures is set to 1.  
% Input variables required are:
% ptycho.m: 4D array that contains the ronchigram for each probe position
% ptycho.s: structured array that contains the centre and radius of the 
% BF disc. In case of binned detectors, it outputs the centre, major and 
% minor axis length of the ellipsis.
% ptycho.IBFimg: matrix of Incoherent Bright Field image
% ptycho.maskPIE: matrix of mask for PIE
% If ptycho.maskPIE is not available, it will calculate it. This will
% require the following variables:
%   ptycho.ObjApt_angle: Probe convergence angle
%   ptycho.theta: matrix of scattering angles
% Output variables added to ptycho are:
% ptycho.BFDiscComX: matrix of center of mass coordinates in x
% ptycho.BFDiscComY: matrix of center of mass coordinates in y
% ptycho.BFDiscCom: matrix of center of mass coordinates 
% ptycho.varfunctions.calculate_center_of_mass: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(ptycho,'maskPIE')
    maskPIE = ptycho.maskPIE;
else
    maskPIE  = ones(size(ptycho.theta)); 
    maskPIE(ptycho.theta>ptycho.ObjApt_angle*1.3)=0;
end


% Quantifying the Centre of Mass of each Ronchigram.

   Com = zeros(2,size(ptycho.m,3),size(ptycho.m,4));
   %solve the problem of M containing negative values for COM calcs.
   Moffset = min(ptycho.m(:)); 
   
   for yy=1:size(ptycho.m,3)
       for xx = 1:size(ptycho.m,4)
           Com(:,yy,xx) = centerOfMass( (ptycho.m(:,:,yy,xx)-Moffset).*maskPIE );
           %edge(M(:,:,xx,1),'prewitt');
       end
   end
   
   BFDiscComX = squeeze(Com(2,:,:))-ptycho.s.Centroid(1);
   BFDiscComY = squeeze(Com(1,:,:))-ptycho.s.Centroid(2);
   BFDiscCom = sqrt(BFDiscComX.^2+BFDiscComY.^2);

   if ptycho.plot_figures
       figure;
       subplot(1,3,1)
       imagesc(BFDiscCom);
       axis image; axis off
       colormap hot;
       colorbar;
       subplot(1,3,2)
       imagesc(BFDiscComX);
       axis image; axis off
       colormap hot;
       colorbar;
       title('X')
       subplot(1,3,3)
       imagesc(BFDiscComY);
       axis image; axis off
       colormap hot;
       colorbar;
       title('Y');
       
       figure;
       imagesc(ptycho.IBFimg); axis image; colormap gray;  axis off;
       hold on
       h=quiver(BFDiscComX,BFDiscComY);
       h.Color = 'blue';
       %         axis image
       axis tight;
       set(gca,'YDir','reverse')
       xlim([1 size(ptycho.m,4)]); ylim([1 size(ptycho.m,3)]);
       hold off
       
       figure;
       h=quiver(BFDiscComX,BFDiscComY);
       h.Color = 'blue';
       axis tight; axis off
       set(gca,'YDir','reverse')
       %         xlim([1 size(m,4)]); ylim([1 size(m,3)]);
       hold off
       axis image
   end

ptycho.BFDiscComX = BFDiscComX;
ptycho.BFDiscComY = BFDiscComY;
ptycho.BFDiscCom = BFDiscCom;
ptycho.varfunctions.calculate_center_of_mass = 1;

end
