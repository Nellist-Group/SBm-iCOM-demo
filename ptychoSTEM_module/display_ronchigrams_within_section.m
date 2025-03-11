function [ptycho] = display_ronchigrams_within_section(ptycho)

% this function needs comments

m = ptycho.m;
IBFimg = ptycho.IBFimg;
DFimg = ptycho.DFimg;
maskPIE = ptycho.maskPIE;

%display all ronchigrams within a specified region

[displayRonchi]=RangeFunc(IBFimg);
figure;
hold all
nx = displayRonchi(3);
ny = displayRonchi(4);

count = 1;
for ix=1:nx
    for iy=1:ny   
        RonchiframeROI(:,:,count) = m(:,:,iy+displayRonchi(2)-1,ix+displayRonchi(1)-1);
        count = count + 1;
        subplot(ny,nx,(iy-1)*nx+ix);
        imagesc(m(:,:,iy+displayRonchi(2)-1,ix+displayRonchi(1)-1).*maskPIE);axis off;
    end
end
hold off
colormap jet

figure;
imagesc(DFimg(displayRonchi(2):displayRonchi(2)+displayRonchi(4),displayRonchi(1):displayRonchi(1)+displayRonchi(3)));
axis image; axis off


end