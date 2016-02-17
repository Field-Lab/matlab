function [metrics,metrics_description] = compute_raster_metrics (raster1,raster2)

%% correlation between PSTHs
convolve=7;
binSz = 1/120;len=1200;
[PSTH1,time]=calculate_psth_fcn2(convolve,binSz,len,raster1);
[PSTH2,time]=calculate_psth_fcn2(convolve,binSz,len,raster2);
R2_PSTHs = R_2_value(PSTH1(0.2*end:end)',PSTH2(0.2*end:end)');

metrics(1) = R2_PSTHs;
metrics_description(1).name = 'correlation between PSTHs';
%% fraction of explained variance 

nr = var(PSTH1-PSTH2);
dr = var(PSTH1);
fev = 1-(nr/dr);

metrics(2) =fev;
metrics_description(2).name = 'fraction of explained variance';

%% correlation between PSTHs from odd/even trials of same stimuli

[PSTH1a,time]=calculate_psth_fcn2(convolve,binSz,len,raster1([1:2:end],:));
[PSTH1b,time]=calculate_psth_fcn2(convolve,binSz,len,raster1([2:2:end],:));
R2_variation1 = R_2_value(PSTH1a(0.2*end:end)',PSTH1b(0.2*end:end)');

metrics(3) = R2_variation1;
metrics_description(3).name = 'correlation between PSTHs from odd/even trials';

%% fraction of explained variance (between odd/even trials of same stimuli)

nr = var(PSTH1a-PSTH1b);
dr = var(PSTH1a);
fev_odd_even = 1-(nr/dr);

metrics(4) = fev_odd_even;
metrics_description(4).name = 'fraction of explained variance between odd and even trials of same stimuli';

end