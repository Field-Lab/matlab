function dist = compare_rasters(spks1,spks2,dt,convolve)
% same condition 
    %convolve=100;
    len =size(spks1,2);
    binSz=dt;
    
    
useRaster1 = spks1;useRaster2 = spks2;
eps= 1;
meanll=[];samell=[];diffll=[];
figure;

for irepeat = 1:size(spks2,1)
    % mean firing rate
    
spkMat =useRaster1([1:irepeat-1,irepeat+1:size(useRaster1,1)],:);
[PSTH_rec,time]=calculate_psth_fcn_gauss(convolve,binSz,len,spkMat);
meanll(irepeat) = ll(useRaster1(irepeat,:),mean(PSTH_rec)*ones(1,size(useRaster1,2)));
plot(mean(PSTH_rec)*ones(1,size(useRaster1,2)));hold on;

% same condition
spkMat =useRaster1([1:irepeat-1,irepeat+1:size(useRaster1,1)],:);
[PSTH_rec,time]=calculate_psth_fcn_gauss(convolve,binSz,len,spkMat);
samell(irepeat) = ll(useRaster1(irepeat,:),PSTH_rec);
plot(PSTH_rec);hold on;

% different condition
spkMat =useRaster2([1:irepeat-1,irepeat+1:size(useRaster2,1)],:);
[PSTH_rec,time]=calculate_psth_fcn_gauss(convolve,binSz,len,spkMat);
diffll(irepeat) = ll(useRaster1(irepeat,:),PSTH_rec);
plot(PSTH_rec);hold on;

dists(irepeat)  =  (diffll(irepeat) - meanll(irepeat))/(samell(irepeat) - meanll(irepeat));
end

dist = mean(dists);

end

function llval = ll(repeat,ll)
llval = -sum(ll) + log(ll)*repeat';
end