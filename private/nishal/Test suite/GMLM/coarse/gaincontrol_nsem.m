%% Load data 
     
cellID=3843;
load(sprintf('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ/organizedspikes_OFFPar_%d.mat',cellID));
train = load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat');
test = load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat');

movie_frames_train = train.NSEMmovie.fitmovie.movie_byblock;
movie_frames_test = test.testmovie.matrix;


nmovies=length(movie_frames_train);
movie_full=(zeros(80,40,(7200)*nmovies));
icnt=1;
for imov=1:nmovies
    imov
movie_full(:,:,icnt:icnt+7200-1) = double(movie_frames_train{imov}.matrix)/255 - 0.5;
icnt = icnt+7200;
% movie_full(:,:,icnt:icnt+3600-1) = movie_frames_test;
% icnt = icnt+3600;
end

spikes_by_block = organizedspikes.block.t_sp_withinblock(2:2:end);

binnedSpikeResponses =[];
for imov=1:length(spikes_by_block)
    
    ss{1}=spikes_by_block{imov}*20000; % Remove spikes which are too soon after change of movie.
    
    ss{1}(ss{1}<20000*30/120)=0;
%     if(rem(imov,2)==1)
%     spkMat = makeSpikeMat(ss,1/120,3600);
%     else
    spkMat = makeSpikeMat(ss,1/120,7200);    
%     end
     
    binnedSpikeResponses = [binnedSpikeResponses,spkMat];
end

%% STA computation 

times = 1:length(binnedSpikeResponses);
t_sp=times(binnedSpikeResponses~=0);

STA=zeros(80,40,30);
for it_sp=t_sp(t_sp>40)
    if(rem(it_sp,120)<60) % second half of each second .. to avoid saccades!
STA = STA+movie_full(:,:,it_sp:-1:it_sp-30+1);
    end
end

STA =STA/numel(t_sp);

average_mov = repmat(mean(movie_full,3),[1,1,30]);


STA2 = STA - average_mov;
STA3 = STA./average_mov;

figure;
subplot(4,1,1);
imagesc(STA(:,:,6)');
colormap gray
axis image

subplot(4,1,2);
imagesc(STA2(:,:,6)');
colormap gray;
axis image

subplot(4,1,3);
imagesc(STA3(:,:,6)');
colormap gray;
axis image


subplot(4,1,4);
imagesc(-STA3(:,:,6)');
colormap gray;
axis image

%% Load mask, etc from WN data
WN_data = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/longBW_offpar/cell_%d.mat',cellID));

maskedMovdd_sliced= 1000*filterMov(movie_full,WN_data.totalMaskSliced,squeeze(WN_data.ttf));

%% play movie for fun
figure
for itime = 1:10000
   imagesc(movie_full(:,:,itime)'.*WN_data.totalMaskSliced');colormap gray;axis image;caxis([-0.5,0.5])
    pause(1/120)
end
%% load WN fits
fitWN = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/longBW_offpar/fitcell_%d.mat',cellID));
%% fit with luminence
movie_full_sliced = mov_slice_flattened(movie_full,WN_data.totalMaskSliced);

[filtered_lum,filter_log,times,luminence_movie]=luminence_signal(-movie_full_sliced);

movinp = maskedMovdd_sliced;
%movinp = prefilter(maskedMovdd_sliced,fitWN.fitGMLM_log{3}.Linear.filter);
%movinp = [];

mov_fit = movinp;
%mov_fit = [movinp;filtered_lum];

[mov_fit,binnedSpikeResponses_tr_use]=cleanup_segments(mov_fit,binnedSpikeResponses,size(filter_log,1),nmovies);


nSU =1;    
interval=1;
%filteredStimDim2 = size(mov_fit,1);
%fitGMLM = fitGMLM_EM_bias(binnedSpikeResponses_tr_use',mov_fit,filteredStimDim2,nSU,interval);


mov_fit = [mov_fit;ones(1,size(mov_fit,2))];
filteredStimDim2 = size(mov_fit,1);
fitGMLM = fitGMLM_EM_accelerate_power2(binnedSpikeResponses_tr_use',mov_fit,filteredStimDim2,nSU,interval,2);

%
idxlum = size(movinp,1)+1;
figure;
for isu=1:nSU
   
filter_total = zeros(size(filter_log,1),1);

for icnt =1:size(filter_log,2)
filter_total  = filter_total+ filter_log(:,icnt) * fitGMLM.Linear.filter{isu}(idxlum+icnt-1);
end

plot(times,filter_total);
hold on;
end

%% Remove bias in GMLM


%% test 
% Get test movie
movie_full_test = double(test.testmovie.matrix)/255 - 0.5;
movie_full_test_sliced = mov_slice_flattened(movie_full_test,WN_data.totalMaskSliced);

maskedMovdd_sliced_test= 1000*filterMov(movie_full_test,WN_data.totalMaskSliced,squeeze(WN_data.ttf));
[filtered_lum_test,filter_log,times,luminence_movie_test]=luminence_signal(-movie_full_test_sliced);

movinp = maskedMovdd_sliced_test;
%movinp = prefilter(maskedMovdd_sliced_test,fitWN.fitGMLM_log{3}.Linear.filter);
%movinp = [];

%mov_test = [maskedMovdd_sliced_test;filtered_lum_test];
mov_test = movinp;
%mov_test = [movinp;filtered_lum_test];

% Predict test response
nTrials=60;
interval=1;
%[predictedResponse,lam] = predictGMLM_bias_lr(fitWN.fitGMLM_log{4},mov_test*2000,nTrials,interval);
%[predictedResponse,lam] = predictGMLM_bias_lr(fitGMLM,mov_test,nTrials,interval);
mov_test = [mov_test;ones(1,size(mov_test,2))];
[predictedResponse,lam] = predictGMLM_gamma2_lr(fitGMLM,mov_test,nTrials,2,interval);
predictedResponse = predictedResponse';

% Test response
spikes_test = organizedspikes.block.t_sp_withinblock(1:2:end);

binnedSpikeResponses_test =[];
for imov=1:length(spikes_test)
    
    ss{1}=spikes_test{imov}*20000; % Remove spikes which are too soon after change of movie.
    ss{1}(ss{1}<20000*30/120)=0;
    spkMat = makeSpikeMat(ss,1/120,3600);    

    binnedSpikeResponses_test = [binnedSpikeResponses_test;spkMat];
end


[x1,y1] = plotSpikeRaster(predictedResponse~=0);
[x2,y2] = plotSpikeRaster(binnedSpikeResponses_test~=0);

saccades = [1:120:120*30];

figure;
plot(x1,y1,'r');
hold on;
plot(x2,y2+60,'k');
hold on;
plot(90+400*(lam-mean(lam))/max(lam-mean(lam)))
hold on;
plot(30+50*(luminence_movie_test - mean(luminence_movie_test))/max(luminence_movie_test - mean(luminence_movie_test)),'b')

hold on; xx=mean(movie_full_test_sliced,1);
plot(0+30*(xx - mean(xx))/max(xx-mean(xx)),'m');

for isacc= 1:length(saccades)
plot(saccades(isacc)*[1,1],[0,120],'--m');
hold on;
end

% % play test movie
% figure
% for itime = 120*15-20:120*15+30
%    imagesc(double(movie_full_test(:,:,itime)').*WN_data.totalMaskSliced');colormap gray;axis image;caxis([-0.5,0.5])
%     pause
% end

% plot movie frames every 1 sec next to each other
[r,c] = find(WN_data.totalMaskSliced==1);

rsz = max(r)-min(r)+1;
csz = (max(c)-min(c)+1);
x = zeros(rsz*2+1,csz*30);
blksz = 120;
for iframe = 1:30
    for isample=1:2
x((isample-1)*rsz+1 + isample-1 : isample-1+ isample*rsz,(iframe-1)*csz+1 : iframe*csz) = movie_full_test(min(r):max(r),min(c):max(c),(iframe-1)*blksz+(isample-1)*118+1);
    end
end

figure;
imagesc(x);
colormap gray
caxis([-0.5,0.5])
axis image

%% Use a/b, ab,a features ..

% train

[features,filter_log,~] = make_features(maskedMovdd_sliced,movie_full_sliced);
[mov_fit,binnedSpikeResponses_tr_use]=cleanup_segments(features,binnedSpikeResponses,size(filter_log,1),nmovies);


nSU =3;    
interval=1;
% filteredStimDim2 = size(mov_fit,1);
%fitGMLM = fitGMLM_EM_bias(binnedSpikeResponses_tr_use',mov_fit,filteredStimDim2,nSU,interval);

mov_fit = [mov_fit;ones(1,size(mov_fit,2))];
filteredStimDim2 = size(mov_fit,1);
fitGMLM = fitGMLM_EM_accelerate_power2(binnedSpikeResponses_tr_use',mov_fit,filteredStimDim2,nSU,interval,2);
%



% test


[features_test,~,~] = make_features(maskedMovdd_sliced_test,movie_full_test_sliced);

features_test = [features_test;ones(1,size(features_test,2))];
[predictedResponse,lam] = predictGMLM_gamma2_lr(fitGMLM,features_test,nTrials,2,interval);
predictedResponse = predictedResponse';

% Test response
spikes_test = organizedspikes.block.t_sp_withinblock(1:2:end);

binnedSpikeResponses_test =[];
for imov=1:length(spikes_test)
    
    ss{1}=spikes_test{imov}*20000; % Remove spikes which are too soon after change of movie.
    ss{1}(ss{1}<20000*30/120)=0;
    spkMat = makeSpikeMat(ss,1/120,3600);    

    binnedSpikeResponses_test = [binnedSpikeResponses_test;spkMat];
end


[x1,y1] = plotSpikeRaster(predictedResponse~=0);
[x2,y2] = plotSpikeRaster(binnedSpikeResponses_test~=0);

saccades = [1:120:120*30];

figure;
plot(x1,y1,'r');
hold on;
plot(x2,y2+60,'k');
hold on;
plot(90+400*(lam-mean(lam))/max(lam-mean(lam)))
hold on;
plot(30+50*(luminence_movie_test - mean(luminence_movie_test))/max(luminence_movie_test - mean(luminence_movie_test)),'b')

hold on; xx=mean(movie_full_test_sliced,1);
plot(0+30*(xx - mean(xx))/max(xx-mean(xx)),'m');

for isacc= 1:length(saccades)
plot(saccades(isacc)*[1,1],[0,120],'--m');
hold on;
end



%% fit div gain ASM

movie_full_sliced = mov_slice_flattened(movie_full,WN_data.totalMaskSliced);

movinp = maskedMovdd_sliced;
%movinp = prefilter(maskedMovdd_sliced,fitWN.fitGMLM_log{3}.Linear.filter);
%movinp = [];


sta = movinp*binnedSpikeResponses'/sum(binnedSpikeResponses);
[filtered_lum,filter_log,times,luminence_movie]=luminence_signal(repmat(sta,[1,size(movie_full_sliced,2)]).*(movie_full_sliced+0.5));
[filtered_lum_inp,~]=cleanup_segments(filtered_lum,binnedSpikeResponses,size(filter_log,1),nmovies);

[movinp,binnedSpikeResponses_tr_use]=cleanup_segments(movinp,binnedSpikeResponses,size(filter_log,1),nmovies);

%movinp = [movinp;ones(1,size(movinp,2))];

nSU =2;    
interval=1;
% 
f=@(x) max(x,0).^2 +0.001;
fd = @(x) 2*max(x,0).^(1);
fdd = @(x)2;
% f=@(x) exp(x);
% fd= @(x) exp(x);
fitASM = fit_divGainASM2_coordinatedescent(movinp,-filtered_lum_inp,binnedSpikeResponses_tr_use,nSU,interval,f,fd,fitWN.fitGMLM_log{nSU}.Linear.filter);
fitASM.f = f;
fitASM.fd = fd;

%
h=figure;
for isu=1:nSU
   
filter_total = zeros(size(filter_log,1),1);

for icnt =1 : size(filter_log,2)
filter_total  = filter_total+ filter_log(:,icnt) * fitASM.params{isu}.a.value(icnt+size(filter_log,2));
end

plot(times,filter_total);
hold on;
end
xlabel('past (in s)')
title('history filter');

%hgexport(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/cellID_3843_tf.eps');

% show kx filter shape
[r,c]=find(WN_data.totalMaskSliced==1);
iidx = reshape(1:numel(WN_data.totalMaskSliced),size(WN_data.totalMaskSliced));
maskiidx = iidx(logical(WN_data.totalMaskSliced));
u=zeros(size(WN_data.totalMaskSliced));
u(maskiidx)=fitASM.params{2}.k.value;
h=figure;
subplot(1,3,1);
imagesc(repelem(u(min(r):max(r),min(c):max(c)),100,100));
colormap gray
axis image
title('Gain model spatial filter')
set(gca,'xTick',[]);
set(gca,'yTick',[]);

% sta b
%sta = movinp*binnedSpikeResponses_tr_use'/sum(binnedSpikeResponses_tr_use);
sta = STA3(:,:,6);
u = sta;
subplot(1,3,2);
imagesc(repelem(u(min(r):max(r),min(c):max(c)),100,100));
colormap gray
axis image
title('STA')
set(gca,'yTick',[]);
set(gca,'xTick',[]);

% a filter
ss = fitASM.params{1}.a.value;
u(maskiidx)=ss;

subplot(1,3,3);
imagesc(repelem(u(min(r):max(r),min(c):max(c)),100,100));
colormap gray
axis image
title('Denominator spatial filter')
set(gca,'xTick',[]);
set(gca,'yTick',[]);

%hgexport(h,'/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/cellID_3843_STA2.eps');


%% Fit GLM, ASM
nSU =2;    
interval=1;
mov_fit = movinp;
filteredStimDim2 = size(mov_fit,1);
fitGMLMexp = fitGMLM_EM_bias(binnedSpikeResponses_tr_use',mov_fit,filteredStimDim2,nSU,interval);

mov_fit = [mov_fit;ones(1,size(mov_fit,2))];
filteredStimDim2 = size(mov_fit,1);
fitGMLMpower = fitGMLM_EM_accelerate_power2(binnedSpikeResponses_tr_use',mov_fit,filteredStimDim2,nSU,interval,2);
%

%% test 

movie_full_test = double(test.testmovie.matrix)/255 - 0.5;
movie_full_test_sliced = mov_slice_flattened(movie_full_test,WN_data.totalMaskSliced);

sta = movinp*binnedSpikeResponses_tr_use'/sum(binnedSpikeResponses_tr_use);
maskedMovdd_sliced_test= 1000*filterMov(movie_full_test,WN_data.totalMaskSliced,squeeze(WN_data.ttf));
[filtered_lum_test,filter_log,times,luminence_movie_test]=luminence_signal(repmat(sta,[1,size(movie_full_test_sliced,2)]).*(movie_full_test_sliced+0.5));

movinp_test = maskedMovdd_sliced_test;


% Predict test response
nTrials=60;
interval=1;

f=@(x) max(x,0).^2+0.001 ;
fd = @(x) 2*max(x,0).^(1);
[predictedResponse,lam,nr,dr] = predict_divGainASM2(fitASM,movinp_test,-filtered_lum_test,interval,f,nTrials);
mov_test = movinp_test;
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

plot([1:length(lam)]/120,250+60*(nr-mean(nr))/max(abs(nr-mean(nr))),'Color',colors(4,:))
hold on;

plot([1:length(lam)]/120,(250+40*(0-mean(nr))/max(abs(0-mean(nr))))*ones(length(lam),1),'Color',colors(4,:));
hold on;

plot([1:length(lam)]/120,350+40*((dr-mean(dr)))/max(abs(dr-mean(dr))),'Color',colors(5,:));
hold on;
plot([1:length(lam)]/120,(350+40*(0-mean(dr))/max(abs(0-mean(dr))))*ones(length(lam),1),'Color',colors(5,:));

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




% plot movie frames every 1 sec next to each other
[r,c] = find(WN_data.totalMaskSliced==1);

rsz = max(r)-min(r)+1;
csz = (max(c)-min(c)+1);
x = zeros(rsz*2+1,csz*30);
blksz = 120;
for iframe = 1:30
    for isample=1:2
x((isample-1)*rsz+1 + isample-1 : isample-1+ isample*rsz,(iframe-1)*csz+1 : iframe*csz) = movie_full_test(min(r):max(r),min(c):max(c),(iframe-1)*blksz+(isample-1)*118+1);
    end
end

figure;
imagesc(x);
colormap gray
caxis([-0.5,0.5])
axis image

%% make a movie of input-response

imov=1;

[tim,pst] = psth_calc(binnedSpikeResponses_test,10,'overlap');

writerObj = VideoWriter(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/NSEM_cellID_%d.avi',cellID),'Uncompressed AVI');
writerObj.FrameRate = 120;

open(writerObj);


for itime = (imov-1)*120+1:imov+30*120-1
  close all  
h=figure;
   subplot(3,4,[1,8]);
    xx= test.testmovie.matrix(:,:,itime);
   xx=(repelem(xx,1,1,3));
   xx(min(r),min(c):max(c),1) =200;
   xx(min(r):max(r),min(c),1) =200;
   xx(min(r):max(r),max(c),1) =200;
   xx(max(r),min(c):max(c),1) =200;
   
    imagesc((permute(xx,[2,1,3])));axis image;caxis([0,255]);colormap gray
  % pause(1/120)



    subplot(3,4,[9,12]);
    
plot(x2/120,y2,'Color',colors(1,:));
hold on; 

xx=mean(movie_full_test_sliced,1);
plot([1:length(lam)]/120,-50+30*(xx - mean(xx))/max(xx-mean(xx)),'Color',colors(7,:));
hold on;
for isacc= 1:length(saccades)
plot(saccades(isacc)*[1,1],[-100,125],'--m');
hold on;
end

plot((itime/120)*[1,1],[-100,125],'--k')

xlim([0,30])
set(gca,'yTick',[]);

hold on;
plot(tim/120,30+90*(pst-mean(pst))/max(abs(pst-mean(pst))),'r')
ylim([-100,125])

F = getframe(h);
writeVideo(writerObj,F);

end

close(writerObj);