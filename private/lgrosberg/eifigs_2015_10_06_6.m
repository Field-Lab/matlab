
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 25);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 27);

[~, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    308,'movieNo',194,'plotElecWaveforms',[73 69 59 55 52 282 278 269],'showElecNums',false);
[~, amps1] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/', 308,'movieNo',202);

[~, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    308,'movieNo',202,'plotElecWaveforms',[73 69 59 55 52 282 278 269],'showElecNums',false);

%%
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 316, 'movieIndex', 25);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 27);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 86, 'movieNo', 299);


playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 197);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',197);
alltrials = 1:25; 
spiketrials = [2 3 4 6 10 11 14 15 16 17 18 19 21 23];
nospiketrials = setxor(alltrials,spiketrials);
figure; 
toplot = squeeze(rawdata(:,93,:))';
subM = repmat(mean(toplot,2),1,size(toplot,2)); 
data = toplot-subM;
figure; 
plot(data(:,spiketrials),'r'); 
hold on; plot(data(:,nospiketrials),'k');

spikedata = squeeze(mean(rawdata(spiketrials,:,:),1)); 
subdata = squeeze(mean(rawdata(nospiketrials,:,:),1));
data = spikedata - subdata; 
amps = max(data,[],2) - min(data,[],2);

