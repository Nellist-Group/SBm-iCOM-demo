function [ptycho] = FT_from_G_to_H(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% this function needs comments

% ptycho.varfunctions.FT_from_G_to_H: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
display('Calculating matrix H= IFT{G}, w.r.t. Kf')

H = zeros(size(ptycho.G_wp),'like',ptycho.G_wp);
% tmp = zeros(1,1,ptycho.ObjSize(1),ptycho.ObjSize(2),'like',G_wp);
for yy=1:size(ptycho.G_wp,3)
    %display(num2str(yy));
    for xx = 1:size(ptycho.G_wp,4)
        H(:,:,yy,xx) = fftshift(ifft2(ifftshift( squeeze(ptycho.G_wp(:,:,yy,xx)) )));
    end
end
display('H calculated..')

ptycho.H = H;
ptycho.varfunctions.FT_from_G_to_H = 1;

end