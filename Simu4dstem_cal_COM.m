function Result = Simu4dstem_cal_COM(dps_4d, dk, rotate_angle)
% calculate COMx, COMy, dCOM and iCOM
% input: dk in mrad, rotate_angle in degree
dpPxlNum = [size(dps_4d,1),size(dps_4d,2)];
dpNum = [size(dps_4d,3),size(dps_4d,4)];

COMx = zeros(dpNum(1),dpNum(2));
COMy = zeros(dpNum(1),dpNum(2));

rotate_angle_rad = rotate_angle/180*pi;
dk_rad = dk * 1e-3; % dk in rad

% iCOM filter
fm = Simu4dstem_generate_mask('center-dark-dot',dpNum,min(dpNum).*0.05);

for ii = 1:dpNum(1)
    for jj = 1:dpNum(2)
        tmpx = dk_rad*((-dpPxlNum(1)/2+0.5)):dk_rad:dk_rad*((dpPxlNum(1)/2-0.5));
        tmpy = dk_rad*((-dpPxlNum(2)/2+0.5)):dk_rad:dk_rad*((dpPxlNum(2)/2-0.5));
        [k1,k2] = meshgrid(tmpy,tmpx);
        COMx(ii,jj) = sum(sum(k1 .* dps_4d(:,:,ii,jj)));
        COMy(ii,jj) = sum(sum(k2 .* dps_4d(:,:,ii,jj)));
    end
end

for kk = 1:length(rotate_angle)
    IcomxT = COMx * cos(rotate_angle_rad(kk)) + COMy * sin(rotate_angle_rad(kk));
    IcomyT = - COMx * sin(rotate_angle_rad(kk)) + COMy * cos(rotate_angle_rad(kk));
    
    % normalize
%     IcomxT = IcomxT - mean(IcomxT(:));
%     IcomyT = IcomyT - mean(IcomyT(:));
    
    Result(kk).COMx = IcomxT;
    Result(kk).COMy = IcomyT;
    
    [gxx,gxy] = gradient(IcomxT);
    [gyx,gyy] = gradient(IcomyT);
    Result(kk).dCOM = gxx + gyy;
    
    fCX = fftshift(fft2(fftshift(IcomxT)));
    fCY = fftshift(fft2(fftshift(IcomyT)));
    KX = size(fCX,2);
    KY = size(fCY,1);
    kxran = linspace(-1, 1, KX);
    kyran = linspace(-1, 1, KY);
    [kx, ky] = meshgrid(kxran, kyran);
    fCKX = fCX .* kx;
    fCKY = fCY .* ky;
    fnum = (fCKX + fCKY);
    fdenom = pi * 2 * (0 + 1j) .* (kx .^ 2 + ky .^ 2);
    fnum(fnum == 0) = 1e-6;
    fdenom(fdenom == 0) = 1e-6;
    fK = fnum ./ fdenom;
    Result(kk).iCOM = real(fftshift(ifft2(ifftshift(fK))));
    
    iCOMf = Simu4dstem_fft2(Result(kk).iCOM);
    iCOMf = abs(iCOMf) .* fm .* exp(1j .* angle(iCOMf));
    Result(kk).iCOMf = real(Simu4dstem_ifft2(iCOMf));
end




end

