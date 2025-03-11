function trans = generate_trans(scan_num,scan_step,start_point,theta)
%GENERATE_TRANS Summary of this function goes here
%   Detailed explanation goes here
% generate trans 
count = 0;
trans = zeros(scan_num(1)*scan_num(2),2);
for jj=1:scan_num(2)
    for ii=1:scan_num(1)
        count=count+1;
        trans(count,1)=ii*scan_step(1)*cos(theta(1)) - jj*scan_step(2)*sin(theta(2)) + start_point(1);
        trans(count,2)=ii*scan_step(1)*sin(theta(1)) + jj*scan_step(2)*cos(theta(2)) + start_point(2);
    end
end
end

