function [dist,dist_sd] = compare_rasters(spks1,spks2,dt,convolve)
% mean firing rate correction was done for raster from other condition!

% same condition 
    %convolve=100;
    len =size(spks1,2);
    binSz=dt;
    
    
useRaster1 = spks1;useRaster2 = spks2;
eps= 1;
meanll=[];samell=[];diffll=[];

rasLen = size(useRaster1,2);
iidxUse = ceil(0.1*rasLen):floor(0.9*rasLen);

figure; 
for irepeat = 1:size(spks2,1)
    % mean firing rate
    
spkMat =useRaster1([1:irepeat-1,irepeat+1:size(useRaster1,1)],:);
[PSTH_rec,time]=calculate_psth_fcn_gauss(convolve,binSz,len,spkMat);PSTH_rec=PSTH_rec(iidxUse);
meanll(irepeat) = ll(useRaster1(irepeat,iidxUse),mean(PSTH_rec)*ones(1,size(useRaster1(:,iidxUse),2)));
plot(mean(PSTH_rec)*ones(1,size(useRaster1(:,iidxUse),2)),'k');hold on;
meanfr = mean(PSTH_rec);

% same condition
spkMat =useRaster1([1:irepeat-1,irepeat+1:size(useRaster1,1)],:);
[PSTH_rec,time]=calculate_psth_fcn_gauss(convolve,binSz,len,spkMat);PSTH_rec=PSTH_rec(iidxUse);
samell(irepeat) = ll(useRaster1(irepeat,iidxUse),PSTH_rec);
plot(PSTH_rec,'r');hold on;

% different condition
spkMat =useRaster2([1:irepeat-1,irepeat+1:size(useRaster2,1)],:);
[PSTH_rec,time]=calculate_psth_fcn_gauss(convolve,binSz,len,spkMat);PSTH_rec=PSTH_rec(iidxUse);
PSTH_rec = PSTH_rec - mean(PSTH_rec) + meanfr;

diffll(irepeat) = ll(useRaster1(irepeat,iidxUse),PSTH_rec);
plot(PSTH_rec,'b');hold on;

dists(irepeat)  =  (diffll(irepeat) - meanll(irepeat))/(samell(irepeat) - meanll(irepeat));
end

dist = mean(dists);
dist_sd = sqrt(var(dists))/sqrt(length(dists));

end

function llval = ll(repeat,ll)
llval = -sum(ll) + log(ll)*repeat';
end