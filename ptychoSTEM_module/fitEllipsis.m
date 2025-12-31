function s = fitEllipsis(BW,method)

% this function needs comments

if strcmp(method,'regionprops')
    
    s = regionprops(BW, 'Orientation', 'MajorAxisLength', 'MinorAxisLength','Centroid');  
    % check the quality of the ellipsis fitting.
%     figure;
%         imshow(BW)
%         hold on
%         phi = linspace(0,2*pi,50);
%         cosphi = cos(phi);
%         sinphi = sin(phi);
% 
%         for k = 1:length(s)
%             xbar = s(k).Centroid(1);
%             ybar = s(k).Centroid(2);
% 
%             a = s(k).MajorAxisLength/2;
%             b = s(k).MinorAxisLength/2;
% 
%             theta = pi*s(k).Orientation/180;
% %             theta = pi*0/180;
%             R = [ cos(theta)   sin(theta)
%                  -sin(theta)   cos(theta)];
% 
%             xy = [a*cosphi; b*sinphi];
%             xy = R*xy;
% 
%             x = xy(1,:) + xbar;
%             y = xy(2,:) + ybar;
% 
%             plot(x,y,'r','LineWidth',2);
%         end
%         hold off
end

end