function omniline = OmniLinePrep(kx,ky,qx,qy,nth)

% this function needs comments

% AX = b;
% A is the OmniMatrix, being prepared line by line using this function.
% X = [C1, C12a, C12b, C23a, C23b, C21a, C21b, C3, C34a, C34b, C32a, C32b]';
% b = phase{G(Kf,Qp)};

% G(Kf,Qp) = |A(Kf)|^2.*delta(Qp) + A(Kf-Qp)A.*(Kf)Psi_s(Qp) +  A(Kf)A.*(Kf+Qp)Psi_s.*(-Qp);

% the phase of A(Kf)A.*(Kf+Qp)Psi_s.*(-Qp)equals to:
% chi(Kf) - chi(Kf+Qp) - phase(Psi_s(-Qp))

kx2 = kx.*kx;
ky2 = ky.*ky;
kx3 = kx2.*kx;
ky3 = ky2.*ky;
kx4 = kx3.*kx;
ky4 = ky3.*ky;

sx = kx + qx;
sy = ky + qy;
sx2 = sx.*sx;
sy2 = sy.*sy;
sx3 = sx2.*sx;
sy3 = sy2.*sy;
sx4 = sx3.*sx;
sy4 = sy3.*sy;


if nth == 1
    omniline = -[1/2.*((sx2+sy2)-(kx2+ky2)),...  %C1
        1/2.*((sx2-sy2)-(kx2-ky2)), sx.*sy-kx.*ky];            %A1
elseif nth == 2
    omniline = -[1/2.*((sx2+sy2)-(kx2+ky2)),...  %C1
        1/2.*((sx2-sy2)-(kx2-ky2)), sx.*sy-kx.*ky,...            %A1
        1/3.*((sx3-3.*sx.*sy2)-(kx3-3.*kx.*ky2)), 1/3.*((3.*sx2.*sy-sy3)-(3.*kx2.*ky-ky3)),...    % A2
        1/3.*((sx3+sx.*sy2)-(kx3+kx.*ky2)), 1/3.*((sy3+sx2.*sy)-(ky3+kx2.*ky)) ];      %Coma B2
elseif nth == 3
    omniline = -[1/2.*((sx2+sy2)-(kx2+ky2)),...  %C1
        1/2.*((sx2-sy2)-(kx2-ky2)), sx.*sy-kx.*ky,...           %A1
        1/3.*((sx3-3.*sx.*sy2)-(kx3-3.*kx.*ky2)), 1/3.*((3.*sx2.*sy-sy3)-(3.*kx2.*ky-ky3)),...    % A2
        1/3.*((sx3+sx.*sy2)-(kx3+kx.*ky2)), 1/3.*((sy3+sx2.*sy)-(ky3+kx2.*ky)),...       %Coma B2
        1/4.*((sx4+sy4+2.*sx2.*sy2)-(kx4+ky4+2.*kx2.*ky2)), ...  % C3
        1/4.*((sx4-6.*sx2.*sy2+sy4)-(kx4-6.*kx2.*ky2+ky4)), 1/4.*((4.*sx3.*sy-4.*sx.*sy3)-(4.*kx3.*ky-4.*kx.*ky3)), ...  %A3
        1/4.*((sx4-sy4)-(kx4-ky4)), 1/4.*((2.*sx3.*sy+2.*sx.*sy3)-(2.*kx3.*ky+2.*kx.*ky3))];  % S3
elseif nth>=4
    error('4th order and above is under development...')
end


end
