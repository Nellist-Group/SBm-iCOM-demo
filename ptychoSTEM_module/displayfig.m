function displayfig(input_image)

    if isreal(input_image)
%         figure; 
        imagesc(input_image); 
        axis image; 
        axis off
%         colorbar; 
    else
%         figure;
        subplot(1,2,1); imagesc(abs(input_image)); 
        axis image; 
        axis off
%         colorbar;
        
        subplot(1,2,2); imagesc(angle(input_image)); 
        axis image;  
        axis off
%         colorbar;
    end
    
%     colormap(violetFire())
    colormap(gray(256))
    set(gcf, 'Color', 'w');
    
end