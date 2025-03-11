function func_aberr = FuncAberrUV2(u,v,aberrcoeff)
% this function outputs the value of the aberration function, Chi, at the
% input position on the detector, u, v, aka Kx, Ky, for a given set of
% aberration coefficients
% use JEOL and CEOS notation
 
% input:
% u: Kx
% v: Ky

% u2 = u.*u;
% u3 = u2.*u;
% u4 = u3.*u;
% 
% v2 = v.*v;
% v3 = v2.*v;
% v4 = v3.*v;
% 
% % aberr are in unit of meter.
aberr.C1   = aberrcoeff(1);
aberr.C12a = aberrcoeff(2);
aberr.C12b = aberrcoeff(3);
aberr.C23a = aberrcoeff(4);
aberr.C23b = aberrcoeff(5);
aberr.C21a = aberrcoeff(6);
aberr.C21b = aberrcoeff(7);
aberr.C3   = aberrcoeff(8);
aberr.C34a = aberrcoeff(9);
aberr.C34b = aberrcoeff(10);
aberr.C32a = aberrcoeff(11);
aberr.C32b = aberrcoeff(12);
% 
% % output:  chi function. in unit of meter*radian.  multiply by 2pi/lambda to get dimensionless
% func_aberr =  1/2.*aberr.C1.*(u2+v2)...
%         + 1/2.*(aberr.C12a.*(u2-v2) + 2.*aberr.C12b.*u.*v)...
%         + 1/3.*(aberr.C23a.*(u3-3*u.*v2) + aberr.C23b.*(3.*u2.*v - v3))...
%         + 1/3.*(aberr.C21a.*(u3+u.*v2) + aberr.C21b.*(v3+u2.*v))...
%         + 1/4.* aberr.C3.*(u4+v4+2*u2.*v2)...
%         + 1/4.* aberr.C34a.*(u4-6.*u2.*v2+v4)...
%         + 1/4.* aberr.C34b.*(4.*u3.*v-4.*u.*v3)...
%         + 1/4.* aberr.C32a.*(u4-v4)...
%         + 1/4.* aberr.C32b.*(2.*u3.*v + 2.*u.*v3);

p.coeffi(1,:) = [aberr.C1*1e9, 0];%nm      C10     df
p.coeffi(2,:) = [aberr.C12a*1e9, aberr.C12b];%    nm      C12     A1
p.coeffi(3,:) = [aberr.C21a*1e9, aberr.C21b];%    nm      C21     B2
p.coeffi(4,:) = [aberr.C23a*1e9, aberr.C23b];%    nm      C23     A2
p.coeffi(5,:) = [aberr.C3*1e6, 0];%    um      C3      C3
p.coeffi(6,:) = [aberr.C32a*1e6, aberr.C32b];%    um      C32     S3
p.coeffi(7,:) = [aberr.C34a*1e6, aberr.C34b];%    um      C34     A3
p.coeffi(8,:) = [0 0];%    um      C41     A4
p.coeffi(9,:) = [0 0];%    um      C43     B4
p.coeffi(10,:) = [0 0];%   um      C45     D4
p.coeffi(11,:) = [0 0];%   mm      C5      C5
p.coeffi(12,:) = [0 0];%   mm      C56     A5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% probe function
% chi function
% func_aber = FuncAberrUV(ptycho.Kx_wp,ptycho.Ky_wp,aberr_input);
% temp1= (-1/2*p.dp.size):(1/2*p.dp.size-1);
% temp1=temp1*1e10*p.lambda/p.dp.dxy/p.dp.size;
% temp2= (-1/2*p.dp.size):(1/2*p.dp.size-1);
% temp2=temp2*1e10*p.lambda/p.dp.dxy/p.dp.size;
% [k1,k2] = meshgrid(temp2,temp1);
% w=k1+1i*k2;
% wi=k1-1i*k2;
% fai=atan2(k2,k1);
w = u + 1j .* v;
wi = u - 1j .* v;
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
func_aberr=real(1/2*w.*wi*C1+...
    1/2*wi.^2*A1+...
    (w.^2).*wi*B2+1/3*wi.^3*A2+...
    1/4*(w.*wi).^2*C3+(w.^3).*wi*S3+1/4*wi.^4*A3+...
    (w.^3).*(wi.^2)*B4+(w.^4).*wi*D4+1/5*wi.^5*A4+...
    1/6*(w.*wi).^3*C5+1/6*wi.^6*A5);

end

