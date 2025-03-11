function Result = Simu4dstem_cal_DPC(dps_4d, dk, rotate_angle, innerAngle, outerAngle)
% calculate COMx, COMy, dCOM and iCOM
% input: dk in mrad, rotate_angle in degree,
%   innerAngle: segmentied detector innerAngle in mrad
%   outerAngle: segmentied detector outerAngle in mrad
dpPxlNum = [size(dps_4d,1),size(dps_4d,2)];
dpNum = [size(dps_4d,3),size(dps_4d,4)];

if dpPxlNum(1) ~= dpPxlNum(2)
    error('Simu4dstem DPC error: diffraction patterns should be square');
end

COMx = zeros(dpNum(1),dpNum(2));
COMy = zeros(dpNum(1),dpNum(2));

rotate_angle_rad = rotate_angle/180*pi;
dk_rad = dk * 1e-3; % dk in rad

% generate masks for 4 detectors
maskRing = Simu4dstem_generate_mask('center-ring',dpPxlNum,innerAngle/dk,outerAngle/dk);
[my,mx] = meshgrid(1:dpPxlNum(1),1:dpPxlNum(2));
mask1 = mx > my;
mask2 = my > dpPxlNum(1) - mx;
maskYp = maskRing & mask1 & mask2;
maskYn = maskRing & (~mask1) & (~mask2);
maskXn = maskRing & mask1 & (~mask2);
maskXp = maskRing & (~mask1) & mask2;

% iCOM filter
fm = Simu4dstem_generate_mask('center-dark-dot',dpNum,min(dpNum).*0.05);

for ii = 1:dpNum(1)
    for jj = 1:dpNum(2)
%         tmpx = dk_rad*((-dpPxlNum(1)/2+0.5)):dk_rad:dk_rad*((dpPxlNum(1)/2-0.5));
%         tmpy = dk_rad*((-dpPxlNum(2)/2+0.5)):dk_rad:dk_rad*((dpPxlNum(2)/2-0.5));
%         [k1,k2] = meshgrid(tmpy,tmpx);
%         COMx(ii,jj) = sum(sum(k1 .* dps_4d(:,:,ii,jj)));
%         COMy(ii,jj) = sum(sum(k2 .* dps_4d(:,:,ii,jj)));
        COMx(ii,jj) = sum(sum(maskXp.*dps_4d(:,:,ii,jj))) - sum(sum(maskXn.*dps_4d(:,:,ii,jj)));
        COMy(ii,jj) = sum(sum(maskYp.*dps_4d(:,:,ii,jj))) - sum(sum(maskYn.*dps_4d(:,:,ii,jj)));
    end
end

for kk = 1:length(rotate_angle)
    IcomxT = COMx * cos(rotate_angle_rad(kk)) + COMy * sin(rotate_angle_rad(kk));
    IcomyT = - COMx * sin(rotate_angle_rad(kk)) + COMy * cos(rotate_angle_rad(kk));
    
    % normalize
%     IcomxT = IcomxT - mean(IcomxT(:));
%     IcomyT = IcomyT - mean(IcomyT(:));
    
    Result(kk).DPCx = IcomxT;
    Result(kk).DPCy = IcomyT;
    
    [gxx,gxy] = gradient(IcomxT);
    [gyx,gyy] = gradient(IcomyT);
    Result(kk).dDPC = gxx + gyy;
    
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
    fdenom = pi * 2 * (0 + 1j) .* ( (kx .^ 2 + ky .^ 2) +  (kx .^ 2 + ky .^ 2) .^ 2);
    fnum(fnum == 0) = 1e-6;
    fdenom(fdenom == 0) = 1e-6;
    fK = fnum ./ fdenom;
    Result(kk).iDPC = real(fftshift(ifft2(ifftshift(fK))));
    
    iDPCf = Simu4dstem_fft2(Result(kk).iDPC);
    iDPCf = abs(iDPCf) .* fm .* exp(1j .* angle(iDPCf));
    Result(kk).iDPCf = real(Simu4dstem_ifft2(iDPCf));
end
end

