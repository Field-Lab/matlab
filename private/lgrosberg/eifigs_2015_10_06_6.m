
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 25);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 27);

[~, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    308,'movieNo',194,'plotElecWaveforms',[73 69 59 55 52 282 278 269],'showElecNums',false);
[~, amps1] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/', 308,'movieNo',202);

[~, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    308,'movieNo',202,'plotElecWaveforms',[73 69 59 55 52 282 278 269],'showElecNums',false);
% There is a slight axonal activation at the lower current amplitude. 

%%
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 316, 'movieIndex', 25);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 27);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 86, 'movieNo', 299);


playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 197);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',197);
eiContour_wLinFit(amps0);
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
eiContour_wLinFit(amps)

playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 221);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',221);


[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',197,'plotElecWaveforms',[238 230 231 232 259 264 276 105 93 86 83 79 76 72 64]);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',221,'plotElecWaveforms',[238 230 231 232 259 264 276 105 93 86 83 79 76 72 64]);

% Contour plot of n1338
[eiM,neuronIdList] = convertEiFileToMatrix('/Volumes/Analysis/2015-10-06-6/data002/data002.ei');
ei = eiM(:,:,find(neuronIdList==1338)); 
eiamps = max(ei,[],2) - min(ei,[],2); 

eiContour_wLinFit(eiamps)
patternNos = find(positions(:,1)< max(positions(:,1))/2 & ...
    positions(:,1)> -max(positions(:,1))/2 &...
    positions(:,2)< max(positions(:,2))/2 & ...
    positions(:,2)> -max(positions(:,2))/2);
figure; 
scatter(positions(:,1),positions(:,2)); 
hold on; scatter(positions(patternNos,1),positions(patternNos,2),'k');
figure;
for p = 1:length(patternNos)
    bundleMeans = getBundleVoltagesAStar('/Volumes/Analysis/2015-10-06-6/data001/', ...
        patternNos(p), 0, 0);
    hold on;
    plot(abs(bundleMeans(1:38,2)),abs(bundleMeans(1:38,1)));
    disp(['finished with pattern ' num2str(p) ' of ' num2str(length(patternNos))]); 
end
xlabel('stimulation amplitude (uA)'); ylabel('mean voltage recorded along bundle'); 

figure; 
scatter(abs(bundleMeans(1:38,2)),abs(bundleMeans(1:38,1)),100,'k','filled');
xlabel('stimulation amplitude (uA)'); ylabel('mean voltage recorded along bundle'); 