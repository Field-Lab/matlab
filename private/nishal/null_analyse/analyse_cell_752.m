
startup_null_analyse_tenessee
b=load('/Volumes/Analysis/nishal/cell752_2.mat');
% Have data for trial start indices.
c=load('/Volumes/Analysis/nishal/cell752.mat');
%%
allElecData=c.allElecData;
waveforms=cell(9,1);
spike_times=cell(9,1);
spike_amps=cell(9,1);
recon_snippets=cell(9,1);

for icell=3:9
waveforms{icell}=b.waveforms{icell};
spike_times{icell}=b.spike_times{icell};
spike_amps{icell}=b.spike_amps{icell};
recon_snippets{icell}=b.recon_snippets{icell};
end

%%
figure;
icnt=0;
for imov=3:9
for icell=1:2
icnt=icnt+1;
subplot(8,2,icnt);
plot(waveforms{imov}{icell});
title(sprintf('Mov: %d Cell : %d ',imov, icell));
end

end
subplot(8,2,15);
plot(c.init_waveform_all_channels');
title('Initial waveforms')

%% Cell 1 in all movies is the target cell, most probably .. 
%% Raster, PSTH, Variance stuff for cell 

dataSpikeTrialsSorted=cell(9,1);
spikeTrialBinnedb=cell(9,1);


time_len = 330*20000;
idx=[1:time_len];
correctCell=1;

for imov=3:9
    imov
    dummyTrain=zeros(1,time_len);
 
    
    
    dummyTrain(round(spike_times{imov}{correctCell})+1)=1;
    trialStartIndices=c.trialStartIndices{imov};
    trialIntervalSamples = 11*20000;
    numTrials=30;
    
    % NOTE: use only 29 samples
    dataTrials=zeros(numTrials-1,trialIntervalSamples);
    for iTrial=1:numTrials-1 % 30th sample less data?
    dataTrials(iTrial,:)=dummyTrain(1,trialStartIndices(iTrial):trialStartIndices(iTrial)+trialIntervalSamples-1);
    end

    dataSpikeTrialsSorted{imov} = dataTrials;
end

binSize=20;
for imov=3:9
    imov
[~,spikeTrialBinnedb{imov}]=getTrialSpikes( dataSpikeTrialsSorted{imov},0.00,'max',binSize);
end

figure; 
for imov=3:9
subplot(7,1,imov-2);
plotSpikeRaster(logical(spikeTrialBinnedb{imov}));
title(sprintf('Movie %d',imov));
end
%%
% PSTH
psthBinSize=10;
psthData=cell(9,1);
timeLog=cell(9,1);
for imov=3:9
[timeLog{imov},psthData{imov}]=  psth_calc((spikeTrialBinnedb{imov}),psthBinSize,'nonoverlap');%mean(dataTrialsSpike{imov},1);
end


figure;
for imov=3:9
subplot(7,1,imov-2);
plot(timeLog{imov}*binSize/(20000),psthData{imov});
ylim([0,25]);
title(sprintf('Movie %d .. %0.4f',imov,sqrt(var(psthData{imov}))));
end

pltIdx=[1:length(psthData{3})];
figure;
plot(timeLog{3}*binSize/(20000),psthData{3}(pltIdx),'b') % PSTH and logPSTH - logPSTH is k.x according to LNP model
hold on
plot(timeLog{4}*binSize/(20000),psthData{4}(pltIdx),'r')
hold on
plot(timeLog{6}*binSize/(20000),psthData{6}(pltIdx),'g')
legend('3','4','6')
title('Movies 3,4,6, PSTH');

figure;
plot(timeLog{7}*binSize/(20000),psthData{7}(pltIdx),'b')
hold on
plot(timeLog{8}*binSize/(20000),psthData{8}(pltIdx),'r')
hold on
plot(timeLog{9}*binSize/(20000),psthData{9}(pltIdx),'g')
legend('7','8','9')
title('Movies 7,8,9, PSTH');

pltIdx=[1:length(psthData{3})];
figure;
semilogy(timeLog{3}*binSize/(20000),psthData{3}(pltIdx),'b') % PSTH and logPSTH - logPSTH is k.x according to LNP model
hold on
semilogy(timeLog{4}*binSize/(20000),psthData{4}(pltIdx),'r')
hold on
semilogy(timeLog{6}*binSize/(20000),psthData{6}(pltIdx),'g')
legend('3','4','6')
title('Movies 3,4,6, log PSTH');

figure;
semilogy(timeLog{7}*binSize/(20000),psthData{7}(pltIdx),'b')
hold on
semilogy(timeLog{8}*binSize/(20000),psthData{8}(pltIdx),'r')
hold on
semilogy(timeLog{9}*binSize/(20000),psthData{9}(pltIdx),'g')
legend('7','8','9')
title('Movies 7,8,9, log PSTH');
%% New STA
% load movies
movies=cell(9,1);
subtract_movies=cell(9,1);
% if I had loaded the movie for rig 4, it would give 1200x64x32 movie .. do transpose,etc accordingly 

[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on/orig_s10.rawMovie',1200,1);
%subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movies{3}=stim-repmat(subtract_movies{3},[1200,1,1]);

[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on/modi_s10.rawMovie',1200,1);
%subtract_movies{4}=mean(stim,1);
subtract_movies{4}=mean(stim,1)*0+127.5;
movies{4}=stim-repmat(subtract_movies{4},[1200,1,1]);


[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on/orig_d_s2.rawMovie',1200,1);
%subtract_movies{5}=mean(stim,1);
subtract_movies{5}=mean(stim,1)*0+127.5;
movies{5}=stim-repmat(subtract_movies{5},[1200,1,1]);

[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on/modi_on_con_s10.rawMovie',1200,1);
%subtract_movies{6}=mean(stim,1);
subtract_movies{6}=mean(stim,1)*0+127.5;
movies{6}=stim-repmat(subtract_movies{6},[1200,1,1]);

[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on_parasol/ori_pon_s10.rawMovie',1200,1);
%subtract_movies{7}=mean(stim,1);
subtract_movies{7}=mean(stim,1)*0+127.5;
movies{7}=stim-repmat(subtract_movies{7},[1200,1,1]);

[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on_parasol/modi_pon_s10.rawMovie',1200,1);
%subtract_movies{8}=mean(stim,1);
subtract_movies{8}=mean(stim,1)*0+127.5;
movies{8}=stim-repmat(subtract_movies{8},[1200,1,1]);

[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-08-20-0/Visual/null_on_parasol/modi_pon_con_s10.rawMovie',1200,1);
%subtract_movies{9}=mean(stim,1);
subtract_movies{9}=mean(stim,1)*0+127.5;
movies{9}=stim-repmat(subtract_movies{9},[1200,1,1]);








%%
staReEst=cell(9,1);

% calculate new STA

staLen=30;
for imov=[3,4,6,7,8,9]
width=size(movies{imov},2);
height=size(movies{imov},3);
staReEst{imov} = zeros(staLen,width,height);
spkCount=0;
for itrial=1:numTrials-1
    itrial
    startSample = round((staLen+1)*20000/120);
    for itime=startSample:10*20000
        if(dataSpikeTrialsSorted{imov}(itrial,itime)==1) % spike at that time
            spkCount=spkCount+1;
            startFrameNo = round(itime*120/20000);
            for i=1:staLen
            staReEst{imov}(i,:,:) = staReEst{imov}(i,:,:)+movies{imov}(startFrameNo-i+1,:,:);
            end
        end
    end
end
staReEst{imov}=staReEst{imov}/spkCount;
end

%% 
stas=cell(1,1);
stas{1}= c.datarun.stas.stas{c.matlab_id};

%staCell=sum(staCell,3);
stas_new=cell(length(stas),1);
for icell=1:length(stas)
    st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new{icell}=st_temp;
end

stas=stas_new;
stas=stas{1};

%plot STA
figure;
for itime=1:30
imagesc(stas(:,:,1,itime));
colormap gray
pause(0.1);
end

xx=[];
for imov=[3,4,6,7,8,9]
innerProd=0;
for iLen=1:staLen
innerProd= innerProd + sum(sum(reshape(stas(:,:,1,iLen),[32,32]).*(reshape(staReEst{imov}(iLen,:,:),[32,32])-127.5)));
end
xx=[xx;[imov,acosd(innerProd/(norm(stas(:))*norm(staReEst{imov}(:))))]]%/(norm(staCell(:))*norm(staReEst{imov}(:)))
end

plot(xx(:,1),xx(:,2),'*')
xlabel('Movie')
ylabel('Inner Product')
%%
cl='rgbcmk';
figure;
icnt=0;
for imov=[4,6,8,9,3,7]
    icnt=icnt+1;
staR=zeros(32,32,30);
for iitime=1:staLen
staR(:,:,iitime)=(staReEst{imov}(iitime,:,:));
end
imov
staNewCross=convn(staR,stas(end:-1:1,end:-1:1,:),'valid')
% Do for whole movie and verify .. 


movi=zeros(width,height,1200+29);
for itime=1:1200
movi(:,:,itime+29)=movies{imov}(itime,:,:);
end
movCross=convn(movi,stas(end:-1:1,end:-1:1,:),'valid');
normMovCross = norm(movCross(:))
plot([1:1200]/120,movCross(:),cl(icnt));
hold on;
end
legend('4','6','8','9','3','7'); % Verify that list above and list in legend same .. 
s%%
figure;
imov=4;
for i=1:staLen
    subplot(2,1,1);
imagesc(reshape(staReEst{imov}(i,:,:),[32,32]));
colormap gray
colorbar
caxis([min(staReEst{imov}(:)),max(staReEst{imov}(:))]);
axis image
%title(sprintf('%d',i));
%title('STA reestimated')
subplot(2,1,2)
imagesc(stas(:,:,1,i));
colormap gray
%title(sprintf('%d',i));
%title('STA')
colorbar
axis image
caxis([min(stas(:)),max(stas(:))])
pause(0.1);
end