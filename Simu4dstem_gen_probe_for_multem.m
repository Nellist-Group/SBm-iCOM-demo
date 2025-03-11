function probe = Simu4dstem_gen_probe_for_multem(jsonData)
%MULTEM_4DSTEM_GEN_EPIE_PROBE_FOR_MULTEM Summary of this function goes here
%   Detailed explanation goes here

if jsonData.nx ~= jsonData.ny
    error('multem_4dstem_gen_ePIE_probe_for_multem: jsonData.nx ~= jsonData.ny');
end

p.beam_voltage = jsonData.voltage;
p.semiangle = jsonData.alpha;
if isfield(jsonData,'inner_angle')
    p.innerangle = jsonData.inner_angle;
else
    p.innerangle = 0;
end
p.dp.dk = jsonData.dk;
p.recon.guess_init_probe = 1;
p.coeffi = zeros(12,2);
p.coeffi(1,:) = [jsonData.defocus/10, 0];%nm      C10     df
p.coeffi(5,:) = [jsonData.cond_lens_c_30_mm*1e3, 0];%    um      C3      C3
if jsonData.cond_lens_c_12_amp_A == 0
    p.coeffi(2,:) = [0,0];
else
    p.coeffi(2,:) = [jsonData.cond_lens_c_12_amp_A/10, jsonData.cond_lens_c_12_angle];%    nm      C12     A1
end
if jsonData.cond_lens_c_23_amp_A == 0
    p.coeffi(4,:) = [0,0];
else
    p.coeffi(4,:) = [jsonData.cond_lens_c_23_amp_A/10, jsonData.cond_lens_c_23_angle];%    nm      C23     A2
end

p.dp.size = jsonData.nx;
p.lambda = 1.23984244/sqrt(p.beam_voltage*(2*510.99906+p.beam_voltage))*1e-9; % m
p.dp.dxy = p.lambda .* 1e10 ./ (p.dp.dk .* 1e-3) / p.dp.size;

p.lambda = 1.23984244/sqrt(p.beam_voltage*(2*510.99906+p.beam_voltage))*1e-9;% m
temp1= (-1/2*p.dp.size):(1/2*p.dp.size-1);
temp1=temp1*1e10*p.lambda/p.dp.dxy/p.dp.size;
temp2= (-1/2*p.dp.size):(1/2*p.dp.size-1);
temp2=temp2*1e10*p.lambda/p.dp.dxy/p.dp.size;
[k1,k2] = meshgrid(temp2,temp1);
w=k1+1i*k2;
wi=k1-1i*k2;
fai=atan2(k2,k1);
%------ calculate aberration ------%
%------ JEOL form ------%
C1   = p.coeffi(1,1)*exp(1i*p.coeffi(1,2)*pi/180)*1e-9;
A1   = p.coeffi(2,1)*exp(2i*p.coeffi(2,2)*pi/180)*1e-9;
B2   = p.coeffi(3,1)*exp(1i*p.coeffi(3,2)*pi/180)*1e-9;
A2   = p.coeffi(4,1)*exp(3*1i*p.coeffi(4,2)*pi/180)*1e-9;
C3   = p.coeffi(5,1)*exp(1i*p.coeffi(5,2)*pi/180)*1e-6;
S3   = p.coeffi(6,1)*exp(2*1i*p.coeffi(6,2)*pi/180)*1e-6;
A3   = p.coeffi(7,1)*exp(4*1i*p.coeffi(7,2)*pi/180)*1e-6;
B4   = p.coeffi(8,1)*exp(1i*p.coeffi(8,2)*pi/180)*1e-6;
D4   = p.coeffi(9,1)*exp(3*1i*p.coeffi(9,2)*pi/180)*1e-6;
A4   = p.coeffi(10,1)*exp(5*1i*p.coeffi(10,2)*pi/180)*1e-6;
C5   = p.coeffi(11,1)*exp(1i*p.coeffi(11,2)*pi/180)*1e-3;
A5   = p.coeffi(12,1)*exp(6*1i*p.coeffi(12,2)*pi/180)*1e-3;
% CEOS notation
kai=real(1/2*w.*wi*C1+...
    1/2*wi.^2*A1+...
    (w.^2).*wi*B2+1/3*wi.^3*A2+...
    1/4*(w.*wi).^2*C3+(w.^3).*wi*S3+1/4*wi.^4*A3+...
    (w.^3).*(wi.^2)*B4+(w.^4).*wi*D4+1/5*wi.^5*A4+...
    1/6*(w.*wi).^3*C5+1/6*wi.^6*A5);
%------ calculate phase shift ------%
% ixmid=floor(p.dp.size(1)/2);
% iymid=floor(p.dp.size(2)/2);
% x=p.dp.dxy(1)*p.dp.size(1)/2;
% y=p.dp.dxy(2)*p.dp.size(2)/2;
% ixoff=floor(x/p.dp.dxy(1))-ixmid;
% iyoff=floor(y/p.dp.dxy(2))-iymid;
% xoff=x-ixoff*p.dp.dxy(1);
% yoff=y-iyoff*p.dp.dxy(1);
% phaseshift=2*pi*(xoff*1e-10*k1+yoff*1e-10*k2)/lamda;
%------ END calculate phase shift ------%
% aberr = -2*pi/lamda*kai+phaseshift;
aberr = -2*pi/p.lambda*kai;
%------ calculate apeture ------%
apeture=double((k1.^2+k2.^2)<=(p.semiangle*1e-3).^2);
apeture = apeture .* double((k1.^2+k2.^2)>=(p.innerangle*1e-3).^2);
% temp=sqrt((k1.^2+k2.^2));
% Nedge=2;
% dEdge = Nedge*p.lambda*1e10/(p.semiangle*1e-3)/p.dp.size(1)/p.dp.dxy;
% ind = find((temp/(p.semiangle*1e-3) > 1-dEdge) & (temp/(p.semiangle*1e-3) < 1+dEdge));
% apeture(ind) = 0.5*(1-sin(pi/(2*dEdge)*(temp(ind)/(p.semiangle*1e-3)-1)));
probe = exp(1j.*aberr).*apeture;
probe = Simu4dstem_fft2(probe);
%     probe=probe/sqrt(sum(sum(abs(probe.*conj(probe)))))*sqrt(sum(sum(abs(dp_ref.*conj(dp_ref)))));
dp_ref = ones(size(probe,1),size(probe,2));
probe=probe/sqrt(sum(sum(abs(probe.*conj(probe)))))*sqrt(sum(sum(abs(sqrt(dp_ref).*conj(sqrt(dp_ref))))));

end

