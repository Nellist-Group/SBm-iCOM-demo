function func_aberr = FuncAberrUV(u,v,aberrcoeff)
% this function outputs the value of the aberration function, Chi, at the
% input position on the detector, u, v, aka Kx, Ky, for a given set of
% aberration coefficients
% this function needs comments
 
% input:
% u: Kx
% v: Ky

u2 = u.*u;
u3 = u2.*u;
u4 = u3.*u;

v2 = v.*v;
v3 = v2.*v;
v4 = v3.*v;

% aberr are in unit of meter.
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

% output:  chi function. in unit of meter*radian.  multiply by 2pi/lambda to get dimensionless
func_aberr =  1/2.*aberr.C1.*(u2+v2)...
        + 1/2.*(aberr.C12a.*(u2-v2) + 2.*aberr.C12b.*u.*v)...
        + 1/3.*(aberr.C23a.*(u3-3*u.*v2) + aberr.C23b.*(3.*u2.*v - v3))...
        + 1/3.*(aberr.C21a.*(u3+u.*v2) + aberr.C21b.*(v3+u2.*v))...
        + 1/4.* aberr.C3.*(u4+v4+2*u2.*v2)...
        + 1/4.* aberr.C34a.*(u4-6.*u2.*v2+v4)...
        + 1/4.* aberr.C34b.*(4.*u3.*v-4.*u.*v3)...
        + 1/4.* aberr.C32a.*(u4-v4)...
        + 1/4.* aberr.C32b.*(2.*u3.*v + 2.*u.*v3);

end

