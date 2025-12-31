function STEMImg = Simu4dstem_cal_stem(dps_4d,dk,stem_angle_mrad)
% dk:  diffraction sampleing in mrad
% input: stem_angle_mrad --> [inner angle, outer angle] in mrad
Mask = Simu4dstem_generate_mask('center-ring',size(dps_4d,[1,2]),stem_angle_mrad(1)./dk,stem_angle_mrad(2)./dk);
STEMImg = zeros(size(dps_4d,3),size(dps_4d,4));
for jj = 1:size(STEMImg,2)
    for ii = 1:size(STEMImg,1)
        STEMImg(ii,jj) = sum(sum(single(dps_4d(:,:,ii,jj)).*Mask));
    end
end

end
