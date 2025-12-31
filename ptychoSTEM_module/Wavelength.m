function w = Wavelength( kev )
% this function needs comments

emass = 510.99906;
hc = 12.3984244;

w = hc/sqrt(kev * (2*emass + kev));
end

