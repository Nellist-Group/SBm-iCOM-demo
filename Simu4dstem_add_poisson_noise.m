function dps = Simu4dstem_add_poisson_noise(dps, AverageCounts)
%ADDPOISNOISEINDP  add poission noise in dps by the average electron count
%on dps
%   dps: 4D array, [:,:,x,y] is a dp indexed as [x,y]
%   AverageCount: a number, the average number of electron counts on each dp

ScanNum = size(dps,[3,4]);

for ii = 1:ScanNum(1)
    for jj = 1:ScanNum(2)
        % normalized dps
        dps(:,:,ii,jj) = dps(:,:,ii,jj) ./ sum(sum(dps(:,:,ii,jj))) .* AverageCounts;
        % add poission noise
        dps(:,:,ii,jj) = poissrnd(dps(:,:,ii,jj));
    end
end



end

