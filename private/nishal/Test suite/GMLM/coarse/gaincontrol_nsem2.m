% Fit with fixed SU from WN data run and a bank of bkg luminence signals .. 
useSU=3;

movie_full_sliced = mov_slice_flattened(movie_full,WN_data.totalMaskSliced);

maskedMovdd_sliced= 1000*filterMov(movie_full,WN_data.totalMaskSliced,squeeze(WN_data.ttf));
movinp = prefilter(maskedMovdd_sliced,fitWN.fitGMLM_log{useSU}.Linear.filter);


maskedMovdd_sliced2= 1000*filterMov(movie_full,WN_data.totalMaskSliced,downsample(squeeze(WN_data.ttf),2));
movinp2 = prefilter(maskedMovdd_sliced2,fitWN.fitGMLM_log{useSU}.Linear.filter);


maskedMovdd_sliced3= 1000*filterMov(movie_full,WN_data.totalMaskSliced,downsample(squeeze(WN_data.ttf),4));
movinp3 = prefilter(maskedMovdd_sliced3,fitWN.fitGMLM_log{useSU}.Linear.filter);

movinp = [movinp;movinp2;movinp3];


[filtered_lum,filter_log,times,luminence_movie]=luminence_signal_bank(movie_full_sliced,fitWN.fitGMLM_log{useSU},WN_data);

%  
% [filtered_lum_inp,~]=cleanup_segments(filtered_lum,binnedSpikeResponses,size(filter_log,1),nmovies);
% [movinp,binnedSpikeResponses_tr_use]=cleanup_segments(movinp,binnedSpikeResponses,size(filter_log,1),nmovies);


[filtered_lum_inp,~]=cleanup_segments_btwn_saccades(filtered_lum,binnedSpikeResponses,size(filter_log,1),nmovies,30);
[movinp,binnedSpikeResponses_tr_use]=cleanup_segments_btwn_saccades(movinp,binnedSpikeResponses,size(filter_log,1),nmovies,30);


%movinp = [movinp;ones(1,size(movinp,2))];

nSU =2*useSU;    
interval=1;
% 
f=@(x) max(x,0).^2 +0.001;
fd = @(x) 2*max(x,0).^(1);
fdd = @(x)2;
% f=@(x) exp(x);
% fd= @(x) exp(x);
fitASM = fit_divGainASM2_coordinatedescent(movinp,filtered_lum_inp,binnedSpikeResponses_tr_use,nSU,interval,f,fd,[]);
fitASM.f = f;
fitASM.fd = fd;


%% 

movie_full_test = double(test.testmovie.matrix)/255 - 0.5;
movie_full_test_sliced = mov_slice_flattened(movie_full_test,WN_data.totalMaskSliced);

sta = movinp*binnedSpikeResponses_tr_use'/sum(binnedSpikeResponses_tr_use);
[filtered_lum_test,filter_log,times,luminence_movie_test]=luminence_signal_bank(movie_full_test_sliced,fitWN.fitGMLM_log{useSU},WN_data)


%movinp_test = maskedMovdd_sliced_test;

maskedMovdd_sliced_test= 1000*filterMov(movie_full_test,WN_data.totalMaskSliced,squeeze(WN_data.ttf));
movinp_test = prefilter(maskedMovdd_sliced_test,fitWN.fitGMLM_log{useSU}.Linear.filter);


maskedMovdd_sliced_test2= 1000*filterMov(movie_full_test,WN_data.totalMaskSliced,downsample(squeeze(WN_data.ttf),2));
movinp_test2 = prefilter(maskedMovdd_sliced_test2,fitWN.fitGMLM_log{useSU}.Linear.filter);

maskedMovdd_sliced_test3= 1000*filterMov(movie_full_test,WN_data.totalMaskSliced,downsample(squeeze(WN_data.ttf),4));
movinp_test3 = prefilter(maskedMovdd_sliced_test3,fitWN.fitGMLM_log{useSU}.Linear.filter);
movinp_test = [movinp_test;movinp_test2;movinp_test3];


% Predict test response
nTrials=60;
interval=1;

f=@(x) max(x,0).^2+0.001 ;
fd = @(x) 2*max(x,0).^(1);
[predictedResponse,lam,nr,dr] = predict_divGainASM2(fitASM,movinp_test,filtered_lum_test,interval,f,nTrials);
mov_test = maskedMovdd_sliced_test;
[predictedResponseexp,lamexp] = predictGMLM_bias_lr(fitGMLMexp,mov_test,nTrials,interval);
predictedResponseexp = predictedResponseexp';
mov_test = [mov_test;ones(1,size(mov_test,2))];
[predictedResponsepower,lampower] = predictGMLM_gamma2_lr(fitGMLMpower,mov_test,nTrials,2,interval);
predictedResponsepower = predictedResponsepower';



% Test response
spikes_test = organizedspikes.block.t_sp_withinblock(1:2:end);

binnedSpikeResponses_test =[];
for imov=1:length(spikes_test)
    
    sss{1}=spikes_test{imov}*20000; % Remove spikes which are too soon after change of movie.
    sss{1}(sss{1}<20000*30/120)=0;
    spkMat = makeSpikeMat(sss,1/120,3600);    

    binnedSpikeResponses_test = [binnedSpikeResponses_test;spkMat];
end


[x1,y1] = plotSpikeRaster(predictedResponse~=0);

[x3,y3] = plotSpikeRaster(predictedResponseexp~=0);

[x4,y4] = plotSpikeRaster(predictedResponsepower~=0);
[x2,y2] = plotSpikeRaster(binnedSpikeResponses_test~=0);

saccades = [1:120:120*30]/120;
colors=distinguishable_colors(10);

figure;
plot(x1/120,y1,'Color',colors(1,:));
hold on;
plot(x2/120,y2+60,'Color',colors(2,:));
hold on;
plot(x3/120,y3-60,'Color',colors(5,:));
hold on;
plot(x4/120,y4+120,'Color',colors(4,:))

plot([1:length(lam)]/120,90+200*(lam-mean(lam))/max(abs(lam-mean(lam))),'Color',colors(1,:));
hold on;

plot([1:length(lam)]/120,250+60*(nr-mean(nr(:)))/max(abs(nr(:)-mean(nr(:)))),'Color',colors(4,:))
hold on;

plot([1:length(lam)]/120,(250+40*(0-mean(nr(:)))/max(abs(0-mean(nr(:)))))*ones(length(lam),1),'Color',colors(4,:));
hold on;

plot([1:length(lam)]/120,350+40*((dr-mean(dr(:))))/max(abs(dr(:)-mean(dr(:)))),'Color',colors(5,:));
hold on;
plot([1:length(lam)]/120,(350+40*(0-mean(dr(:)))/max(abs(0-mean(dr(:)))))*ones(length(lam),1),'Color',colors(5,:));

% hold on;
% plot([1:length(lam)]/120,-100+50*(filtered_lum_test - mean(filtered_lum_test))/max(abs(filtered_lum_test - mean(filtered_lum_test))),'Color',colors(6,:));

hold on; 
xx=mean(movie_full_test_sliced,1);
plot([1:length(lam)]/120,-200+30*(xx - mean(xx))/max(xx-mean(xx)),'Color',colors(7,:));

for isacc= 1:length(saccades)
plot(saccades(isacc)*[1,1],[-250,450],'--m');
hold on;
end

xlim([0,30])


%% Between-Saccades STA




