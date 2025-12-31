function [ptycho] = parallel_FT_from_G_to_H(ptycho)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% this function needs comments

% ptycho.varfunctions.FT_from_G_to_H: Flag to indicate the function has been run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
display('Calculating matrix H= IFT{G}, w.r.t. Kf')

size4D = size(ptycho.G_wp);
size3D = [size(ptycho.G_wp,1),size(ptycho.G_wp,2), ptycho.ObjSize(1)*ptycho.ObjSize(2)];

Gres = reshape(ptycho.G_wp,size3D);

H = zeros(size3D);

parfor yy=1:size3D(3)
    %display(num2str(yy));
       H(:,:,yy) = fftshift(ifft2(ifftshift( squeeze(Gres(:,:,yy)) )));
end

H = reshape(H,size4D);

display('H calculated..')

ptycho.H = H;
ptycho.varfunctions.FT_from_G_to_H = 1;

end