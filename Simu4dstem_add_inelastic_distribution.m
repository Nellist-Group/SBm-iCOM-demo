function dpse = Simu4dstem_add_inelastic_distribution(dps,dk,energyDistribution,Estart,dE,E0)
%SIMU4DSTEM_ADD_INELASTIC_DISTRIBUTION generate the inelastic distribution
% on diffraction patters, by the EELS spectrum
%   input:
%       dps --> diffraction patterns, 3D matrix. 
%       dk --> pixel size in mrad for dps
%       energyDistribution --> eels spectrum, simulated or experimental, start from 0 eV.
%           the counts on different energy chanels
%       Estart --> the energy loss of first channel in eV
%       dE --> the chanel size in eV for energy Distribution
%       E0 --> original Energy in keV. (Accelerate voltage for EM)
%   output:
%       dpse --> the distribution of electrons on dp considering inelastic
%           distribution. It is a dp normalized to intensity sum as 1.

deltaE = (0:length(energyDistribution)-1) .* dE + Estart; % in eV, energy loss for each channel
thetaE = (deltaE ./ E0 ./ 2); % mrad, charactoristic scattering angle for each channel
eRatio = energyDistribution ./ sum(energyDistribution(:));

% init dpse
dpse = zeros(size(dps),'like',dps);

% init wait bar
wb = Simu4dstem_WaitBarStr;
wb.Init('calculating inelastic distribution...');
for iiE = length(thetaE):-1:1
    if iiE >= 2
        gama = thetaE(iiE);
        patchR = min(round(min(size(dps,1),size(dps,2))./2.5),2.*ceil(gama./dk)); % in pixel
        [gy,gx] = meshgrid((-patchR:patchR).*dk, (-patchR:patchR).*dk);
        theta = sqrt(gx.^2 + gy.^2);
        lorentzDstb = gama ./ (theta.^2 + gama.^2) ./ pi;
        lorentzDstb = lorentzDstb ./ sum(lorentzDstb(:));
    else
        patchR = 0;
        lorentzDstb = 1;
    end
    for iidp = 1:size(dps,3)
        dpt = zeros(size(dps,1),size(dps,2),'like',dps);
        dpt = padarray(dpt,[patchR,patchR]);
        dpn = dps(:,:,iidp) ./ sum(sum(dps(:,:,iidp))) .* eRatio(iiE);
        for iikx = 1:size(dps,1)
            for iiky = 1:size(dps,2)
                dstb = lorentzDstb .* dpn(iikx,iiky);
                cutRangeDp = [iikx,iikx+patchR.*2;
                              iiky,iiky+patchR.*2];
                dpt(cutRangeDp(1,1):cutRangeDp(1,2),cutRangeDp(2,1):cutRangeDp(2,2)) =...
                    dstb + dpt(cutRangeDp(1,1):cutRangeDp(1,2),cutRangeDp(2,1):cutRangeDp(2,2));
            end
        end
%         dpt = poissrnd(dpt .* epdp);
        dpse(:,:,iidp) = dpse(:,:,iidp) + dpt(patchR+1:size(dpt,1)-patchR, patchR+1:size(dpt,2)-patchR);
    end
    wb.Update(1-(iiE./length(thetaE)));
end
wb.Update(1);
end

