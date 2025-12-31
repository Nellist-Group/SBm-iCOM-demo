function chival = ChiuvFunc( u, v, aberrfit )

% this function needs comments

% u - x-direction (kx)
% v - y-direction (ky)

u2 = u.*u;
v2 = v.*v;
u3 = u.*u.*u;
v3 = v.*v.*v;
u4 = u.*u.*u.*u;
v4 = v.*v.*v.*v;

% aberrfit = [C0, C01a, C01b, C1, C12a, C12b, C23a, C23b, C21a, C21b, C3, C34a, C34b, C32a, C32b]';
    C01a    = 0;
    C01b    = 0;
    C1      = 0;
    C12a    = 0;
    C12b    = 0;
    C23a    = 0;
    C23b    = 0;
    C21a    = 0;
    C21b    = 0;
    C3      = 0;
    C34a    = 0;
    C34b    = 0;
    C32a    = 0;
    C32b    = 0;
    
if length(aberrfit) >= 6
%     C01a    = aberrfit(2);
%     C01b    = aberrfit(3);
    C1      = aberrfit(4);
    C12a    = aberrfit(5);
    C12b    = aberrfit(6);
elseif length(aberrfit) >= 10
    C23a    = aberrfit(7);
    C23b    = aberrfit(8);
    C21a    = aberrfit(9);
    C21b    = aberrfit(10);
elseif length(aberrfit) >= 15
    C3      = aberrfit(11);
    C34a    = aberrfit(12);
    C34b    = aberrfit(13);
    C32a    = aberrfit(14);
    C32b    = aberrfit(15);
end


chival = C01a.*u + C01b.*v + 1/2.*C1.*(u2+v2)...
        + 1/2.*(C12a.*(u2-v2) + 2.*C12b.*u.*v)...
        + 1/3.*(C23a.*(u3-3*u.*v2) + C23b.*(3.*u2.*v - v3))...
        + 1/3.*(C21a.*(u3+u.*v2) + C21b.*(v3+u2.*v))...
        + 1/4.* C3.*(u4+v4+2*u2.*v2)...
        + 1/4.* C34a.*(u4-6.*u2.*v2+v4)...
        + 1/4.* C34b.*(4.*u3.*v-4.*u.*v3)...
        + 1/4.* C32a.*(u4-v4)...
        + 1/4.* C32b.*(2.*u3.*v + 2.*u.*v3);


end

