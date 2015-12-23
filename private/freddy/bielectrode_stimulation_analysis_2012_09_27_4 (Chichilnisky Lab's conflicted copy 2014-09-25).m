% Script to plot the stimulation thresholds for ON Parasol cells in
% 2012-09-24-3 512 stimulation data for bi-electrode stimulation
% experiments

% LGrosberg June 2014

%% Load data
dataPath = '/Volumes/Analysis/2012-09-27-4/data007/data007';%'/Volumes/Analysis/2012-09-24-3/data007/data007';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
positions = datarun.ei.position;


elecStimPath = '/Volumes/Analysis/2012-09-27-4/data009/'; %'/Volumes/Analysis/2012-09-24-3/data008/';
dirInfo = dir(elecStimPath);

elecsByPattern = full(getElecsFromPatternFiles([elecStimPath 'pattern_files/'])); %load the matrix that gives electrode numbers indexed by pattern

% Load thresholds for bidirection axon bundle activation (determined by visual inspection of data)
%temp = load('~/Development/matlab-standard/private/lauren/MATLAB_code/analysis/dataset-specific/axonBundleThresholds_byPattern_2012_09_24_3_data008.mat');
%axonBundleThresholds = temp.axonBundleThresholds_byPattern_2012_09_24_3_data008;
%%
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:IterationLimitExceeded');
% ON Parasol cells in quadrant 3
neurons = [1248 1281 1504 1804 1880 1729 1099 1427 1463 1968 2015 ...
            2299 2393 2480 2677 2842]; %2641, 1909, and 2162 not included due to low SNR/edge of quadrant

%[2 33 168 197 229 303 332 436 601 647 919 4909 6603 6633 6769 6936 6995 7008 7173 7218 7474 7507 7549]; %7607

activationThresh_bielec1    = zeros(size(neurons,2),512);
activationThresh_bielec2    = zeros(size(neurons,2),512);
activationThresh_singleElec = zeros(size(neurons,2),512);

thresholdDiff_bielec1 = zeros(size(neurons,2),4); %4 was highest number of patterns per neuron -- empirical
thresholdDiff_singleelec = zeros(size(neurons,2),4);

somaStimThreshs_bielec1 = zeros(size(neurons,2),4);
somaStimThreshs_singleelec = zeros(size(neurons,2),4);

patternNos_bielec1 = zeros(size(neurons,2),4);
patternNos_singleelec = zeros(size(neurons,2),4);
% x = 1;
% z = 1;
for j = 1:1:size(neurons,2)         % For all ON parasol cells
    neuron  = neurons(j);           % gets the id number of the ith neuron
    a = 1; b = 1;
    for  n = 3:size(dirInfo,1)
        if ~dirInfo(n).isdir
            fname = dirInfo(n).name;
            i = find(fname=='_',2,'first');
            if strcmp(['n' num2str(neuron)],fname(i(1)+1:i(2)-1))
                temp = load([elecStimPath fname]);
                elecResp = temp.elecResp; clear temp;
                
                if size(elecResp.stimInfo.electrodes,2) == 2     % Two electrode stimulation
                    % Polarity 1: Left electrode,  2:-3:1, Right electrode, -2:3:-1
                    % Polarity 2: Left electrode, -2:3:-1, Right electrode,  2:-3:1
                    if positions(elecResp.stimInfo.electrodes(1),1) < positions(elecResp.stimInfo.electrodes(2),1)
                        leftElectrode = elecResp.stimInfo.electrodes(1);
                        if min(elecResp.stimInfo.pulseVectors{1}(1,2,:)) < min(elecResp.stimInfo.pulseVectors{1}(2,2,:))
                            stimType = 'bielectrode_P1';
                        else
                            stimType = 'bielectrode_P2';
                        end
                    elseif positions(elecResp.stimInfo.electrodes(1),1) > positions(elecResp.stimInfo.electrodes(2),1)
                        leftElectrode = elecResp.stimInfo.electrodes(2);
                        if min(elecResp.stimInfo.pulseVectors{1}(1,2,:)) > min(elecResp.stimInfo.pulseVectors{1}(2,2,:))
                            stimType = 'bielectrode_P1';
                        else
                            stimType = 'bielectrode_P2';
                        end
                    end
                elseif size(elecResp.stimInfo.electrodes,2) == 1 % Single electrode stimulation
                    stimType = 'singleElectrode';
                end
                
                
                responseProb = elecResp.analysis.successRates;
                stimAmps     = abs(elecResp.stimInfo.stimAmps);
                
                % Define function that will be used to fit data
                % (F is a vector of fitting parameters)
                f = @(F,x) (1 +exp(-F(1)*(x - F(2)))).^(-1); % sigmoid
                F_fitted = nlinfit(stimAmps,responseProb,f,[1 1]);
                
                % Plot data fit
                y = f(F_fitted,stimAmps);
                
                % Find the threshold voltage for 50% response probability
                ii = find(y>0.5, 1,'first');
                if ii<size(stimAmps,1)-1
                    xx = stimAmps(ii-1):0.001:stimAmps(ii+2);
                    yy = f(F_fitted,xx);
                    
                    threshVoltage = xx(find(yy>0.5,1,'first'));
                    
                    switch stimType
                        case 'bielectrode_P1'
                            activationThresh_bielec1(j,leftElectrode) = threshVoltage;
                            %                             bielecDiff(x) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            thresholdDiff_bielec1(j,a) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            patternNos_bielec1(j,a) = elecResp.stimInfo.patternNo;
                            somaStimThreshs_bielec1(j,a) = threshVoltage;
                            a = a+1;
                            %                             x = x + 1;
                            disp(['Bielectrode hit for neuron ' num2str(neurons(j)) ' patternNo: ' num2str(elecResp.stimInfo.patternNo)]);
                        case 'bielectrode_P2'
                            activationThresh_bielec2(j,leftElectrode) = threshVoltage;
                        case 'singleElectrode'
                            activationThresh_singleElec(j,elecResp.stimInfo.electrodes) = threshVoltage;
                            %                             singleelecDiff(z) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            thresholdDiff_singleelec(j,b) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            patternNos_singleelec(j,b) = elecResp.stimInfo.patternNo;
                            somaStimThreshs_singleelec(j, b) = threshVoltage;
                            b = b + 1 ;
                            %                             z = z + 1;
                            disp(['Single electrode hit for neuron ' num2str(neurons(j)) ' patternNo: ' num2str(elecResp.stimInfo.patternNo)]);
                    end
                    
                    %                     actThresholds(elecResp.stimInfo.electrodes) = threshVoltage;
                    %
                    %Plotting
                    %                     figure; plot(stimAmps, responseProb,'.'); ylim([0 1])
                    %                     xlabel('Current (\muA)'); ylabel('response probability');
                    %                     title(sprintf('neuron %d; stimulating electrode %d',neuron,elecResp.stimInfo.electrodes));
                    %
                    %                     hold on; plot(stimAmps,y,'g');
                    %                     hold on; plot(threshVoltage,yy(find(yy>0.5,1,'first')),'ro');
                    %                     text(threshVoltage+0.1,0.5,[num2str(threshVoltage) '\muA'])
                    %
                    
                    
                end
                
                
                
            end
        end
    end
end


%% Plotting ..
groupnames = cell(size(neurons, 2), 1);
for i = 1:size(neurons, 2)
    groupnames{i} = num2str(neurons(i));
end

%groupnames = {'2';'33';'168';'197';'229';'303';'332';'436';'601';'647';'919';'4909';'6603';'6633';'6769';'6936';'6995';'7008';'7173';'7218';'7474';'7507';'7549'}; %'7607'
bw_colormap = [1 0 0];
gridstatus = 'y';
figure;
barweb(thresholdDiff_bielec1, zeros(size(thresholdDiff_bielec1)), [], groupnames, 'Bi-electrode stimulation difference from bundle activation', 'neuron', 'soma threshold - bundle threshold (uA)', bw_colormap, gridstatus, []);

figure;
barweb(thresholdDiff_singleelec, zeros(size(thresholdDiff_singleelec)), [], groupnames, 'Single electrode stimulation difference from bundle activation', 'neuron', 'soma threshold - bundle threshold (uA)', bw_colormap, gridstatus, []); whitebg('white');

thresholdDiff_singleelec(thresholdDiff_singleelec==0) = NaN;
thresholdDiff_bielec1(thresholdDiff_bielec1==0)       = NaN;

somaStimThreshs_singleelec(somaStimThreshs_singleelec == 0) = NaN;
somaStimThreshs_bielec1(somaStimThreshs_bielec1 == 0) = NaN;

patternNos_bielec1(patternNos_bielec1==0) = NaN;
patternNos_singleelec(patternNos_singleelec==0) = NaN;

maxPatterns_bielec1 = NaN(size(neurons,2),1);
maxPatterns_singleelec = NaN(size(neurons,2),1);

maxThreshs_singleelec = NaN(size(neurons,2), 1);
maxThreshs_bielec1 = NaN(size(neurons,2), 1);

[singleelecvalues, singleelecindices] = min(thresholdDiff_singleelec,[],2);
[bielec1values, bielec1indices] = min(thresholdDiff_bielec1,[],2);

for i = 1:1:size(singleelecindices)
    maxPatterns_bielec1(i) = patternNos_bielec1(i, bielec1indices(i));
    maxThreshs_bielec1(i) = somaStimThreshs_bielec1(i, bielec1indices(i));
    maxPatterns_singleelec(i) = patternNos_singleelec(i, singleelecindices(i));
    maxThreshs_singleelec(i) = somaStimThreshs_singleelec(i, singleelecindices(i));
end



threshDiffs = cat(2,singleelecvalues,bielec1values);
figure;
bw_legend = {'single electrode stim','bi-electrode stim'};
bw_colormap = [1 0 0; 0 0 1];
barweb(threshDiffs, zeros(size(threshDiffs)), [], groupnames, 'Somatic activation of ON Parasols VS bundle activation', 'neuron', 'soma threshold - bundle threshold (uA)', bw_colormap, gridstatus, bw_legend);

threshDiffDiffs = singleelecvalues-bielec1values;

max(thresholdDiff_singleelec,[],2)
%Set max/min color limits for plotting thresholds
act_max = 4;
act_min = 0.5;

tmp = figure;
rgbVals = colormap(tmp,gray);
rgbVals = flipud(rgbVals);
rgbVals(end,:) = 0;
close(tmp);

% Overlay rf fits and IDs
figure; set(gcf,'Color',[1 1 1]); hold on;
for n = 1:1:size(neurons,2)
    currentActThresh_single = activationThresh_singleElec(n,:);
    [~,J,val] = find(currentActThresh_single);
    minThresh = min(val);
    bestStimElec = J(find(val == min(val)));
    if ~isempty(minThresh)
        rgbIndex = (minThresh - act_min)/(act_max-act_min); %scale between 0 and 1
        rgbIndex = round(rgbIndex * size(rgbVals,1));
        if rgbIndex<1; rgbIndex = 1; elseif rgbIndex>size(rgbVals,1); rgbIndex = size(rgbVals,1); end
        
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color',rgbVals(rgbIndex,:));
    else
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color','k');
    end
    
    disp(['activation threshold for neuron ' num2str(neurons(n)) ' = ' num2str(minThresh) ' uA at electrode ' num2str(bestStimElec) ' (single electrode stim)']);
end
axis image; axis off;
colorbar; colormap(rgbVals); caxis([act_min act_max]);
title('ON parasol activation thresholds, single electrode stimulation');



figure; set(gcf,'Color',[1 1 1]); hold on;
for n = 1:1:size(neurons,2)
    currentActThresh_bielec1 = activationThresh_bielec1(n,:);
    [~,J,val] = find(currentActThresh_bielec1);
    minThresh = min(val);
    bestStimElec = J(find(val == min(val)));
    if ~isempty(minThresh)
        rgbIndex = (minThresh - act_min)/(act_max-act_min); %scale between 0 and 1
        rgbIndex = round(rgbIndex * size(rgbVals,1));
        if rgbIndex<1; rgbIndex = 1; elseif rgbIndex>size(rgbVals,1); rgbIndex = size(rgbVals,1); end
        
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color',rgbVals(rgbIndex,:));
    else
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color','k');
    end
    
    disp(['activation threshold for neuron ' num2str(neurons(n)) ' = ' num2str(minThresh) ' uA at electrode ' num2str(bestStimElec) ' (bi-electrode stim, polarity 1)']);
end
axis image; axis off;
colorbar; colormap(rgbVals); caxis([act_min act_max]);
title('ON parasol activation thresholds, bi-electrode stimulation, polarity 1');



figure; set(gcf,'Color',[1 1 1]); hold on;
for n = 1:1:size(neurons,2)
    currentActThresh_bielec2 = activationThresh_bielec2(n,:);
    [~,J,val] = find(currentActThresh_bielec2);
    minThresh = min(val);
    bestStimElec = J(find(val == min(val)));
    if ~isempty(minThresh)
        rgbIndex = (minThresh - act_min)/(act_max-act_min); %scale between 0 and 1
        rgbIndex = round(rgbIndex * size(rgbVals,1));
        if rgbIndex<1; rgbIndex = 1; elseif rgbIndex>size(rgbVals,1); rgbIndex = size(rgbVals,1); end
        
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color',rgbVals(rgbIndex,:));
    else
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'k','fill_color','k');
    end
    
    disp(['activation threshold for neuron ' num2str(neurons(n)) ' = ' num2str(minThresh) ' uA at electrode ' num2str(bestStimElec) ' (bi-electrode stim, polarity 2)']);
end
axis image; axis off;
colorbar; colormap(rgbVals); caxis([act_min act_max]);
title('ON parasol activation thresholds, bi-electrode stimulation, polarity 2');


%% Plotting EIs
printElecs = 1;
neuronID = 919;

index = find(neurons==neuronID);
singleelecpattern = maxPatterns_singleelec(index);  %get the pattern number for the best single electrode difference for this neuron
bielec1pattern = maxPatterns_bielec1(index);        %do the same for the double electrode pattern
singleelec = elecsByPattern(singleelecpattern, 1);
bielecs1 = elecsByPattern(bielec1pattern, :);
disp(['Single electrode pattern: ' num2str(singleelecpattern) '    Double electrode pattern: ' num2str(bielec1pattern)]);
disp(['Single electrode thresh:  ' num2str(activationThresh_singleElec(index, singleelec))  ' Double electrode pattern: ']);
singleelec = elecsByPattern(singleelecpattern, 1); %get the electrode used for this pattern
bielecs1 = elecsByPattern(bielec1pattern, :);      %get the  vector for the electrodes used in the double electrode stimulation


EI = datarun.ei.eis{get_cell_indices(datarun, neuronID)};
eiamps = max(EI') - min(EI');
eiamps2 = eiamps*8;
[eiMax, maxIndex] = max(eiamps);
figure('position', [0 500 900 500]); scatter(positions(:,1),positions(:,2),eiamps2,[0.5,0.5,1],'filled'); axis image; title(['neuron ' num2str(neuronID) ': '  num2str(threshDiffDiffs(index))]); whitebg('black'); axis off; set(gcf, 'InvertHardCopy', 'off');
%figure('position', [0 500 900 500]); scatter(positions(:,1),positions(:,2),350, eiamps,'filled'); colorbar; caxis([0 25]); axis image; title(['neuron ' num2str(neuronID) ': '  num2str(threshDiffDiffs(index))]); whitebg('black'); axis off; set(gcf, 'InvertHardCopy', 'off');
hold on; scatter(positions(maxIndex, 1), positions(maxIndex, 2), 350, 'white', 'filled');
if printElecs
    for x = 1:512
        if x==singleelec || any(x==bielecs1)
            th=text(positions(x,1),positions(x,2),num2str(x),'HorizontalAlignment','center');
            if x==singleelec && any(x==bielecs1)
                set(th, 'color', [1 0.5 0.5]);
            elseif x == singleelec
                set(th, 'color', [1 0.5 1]);
            elseif any(x==bielecs1)
                set(th, 'color', [0.5 1 0.5]);
            end
        end
        
    end
end
%%

figure; scatter(1:23,threshDiffDiffs); line([0 23], [0 0]); whitebg('black'); set(gcf, 'InvertHardCopy', 'off');
colorsForPlot = zeros(512, 3);

threshDiffDiffFull = zeros(512, 1);


for n = 1:size(neurons,2)
    singleelecpattern = maxPatterns_singleelec(n); %get the pattern number for the best single electrode difference for this neuron
    
    if ~isnan(singleelecpattern)
        singleelec = elecsByPattern(singleelecpattern, 1); %get the electrode used for this pattern
        if ~isnan(threshDiffDiffs(n))
            threshDiffDiffFull(singleelec) = threshDiffDiffs(n);
        end
        if threshDiffDiffs(n) > 0
            colorsForPlot(singleelec, :) = [1 0 0];
        else
            colorsForPlot(singleelec, :) = [0 1 0];
        end
    end
end

threshDiffDiffFull(threshDiffDiffFull == 0) = 0.001;

figure('position', [0 500 900 500]); scatter(positions(:,1),positions(:,2),abs(threshDiffDiffFull)*500,colorsForPlot,'filled'); axis image; title('diffdiff'); axis off; set(gcf, 'InvertHardCopy', 'off'); set(gcf, 'InvertHardCopy', 'off');

%% Get electrode numbers for best case by neuron

neuronNumber = 33;
index = get_cell_indices(datarun, neuronNumber);
singleelecpattern = maxPatterns_singleelec(index);
bielec1pattern = maxPatterns_bielec1(index);
singleelec = elecsByPattern(singleelecpattern, 1);
bielecs1 = elecsByPattern(bielec1pattern, :);

%% soma activation vs axon activations

maxPatterns_bielec1(maxPatterns_bielec1 == 1) = NaN;
maxPatterns_singleelec(maxPatterns_singleelec == 1) = NaN;

bundleThreshs_copy = axonBundleThresholds;
bundleThreshs_copy(1) = NaN;

maxPatterns_bielec1(isnan(maxPatterns_bielec1)) = 1;
maxPatterns_singleelec(isnan(maxPatterns_singleelec)) = 1;

figure; plot(cat(2, maxThreshs_bielec1, bundleThreshs_copy(maxPatterns_bielec1))); whitebg('black'); title('bi-electrode axon vs soma thresh'); axis([1 23 0 4.5]); set(gcf, 'InvertHardCopy', 'off');
figure; plot(cat(2, maxThreshs_singleelec, bundleThreshs_copy(maxPatterns_singleelec))); whitebg('black'); title('single electrode axon vs soma thresh'); axis([1 23 0 4.5]); set(gcf, 'InvertHardCopy', 'off');

%% messing with movies

singleElecPatterns = 241:368;
biElecPatterns1 = 1:2:239;
biElecPatterns2 = 2:2:240;

patternNos = biElecPatterns2;
display = false;

bundleMaxes = zeros(368, 1);

for patternNo = patternNos
    bundleMean = getBundleVoltages(patternNo, display);
    bundleMaxes(patternNo) = max(abs(bundleMean(:,1)));
end
disp('finished analyzing patterns');
electrodesByPattern = getElecsFromPatternFiles([elecStimPath 'pattern_files/']);

electrodeNumberingScheme = zeros(size(electrodesByPattern, 1), 1);

electrodeNumberingScheme(1:2:15) = electrodesByPattern(1:2:15, 1);
electrodeNumberingScheme(17:2:127) = electrodesByPattern(17:2:127, 2);
electrodeNumberingScheme(129:2:239) = electrodesByPattern(129:2:239, 1);

electrodeNumberingScheme(2:2:16) = electrodesByPattern(2:2:16, 2);
electrodeNumberingScheme(18:2:128) = electrodesByPattern(18:2:128, 1);
electrodeNumberingScheme(130:2:240) = electrodesByPattern(130:2:240, 2);

electrodeNumberingScheme(241:368) = electrodesByPattern(241:368, 1);

% 
% vbye1 = zeros(1, 512);
% vbye1(electrodeNumberingScheme(biElecPatterns1)) = bundleMaxes(biElecPatterns1);
% 
% vbye2 = zeros(1, 512);
% vbye2(electrodeNumberingScheme(biElecPatterns2)) = bundleMaxes(biElecPatterns2);
% 
% test = mean(abs(cat(1, vbye1, vbye2)), 1);
% [~, testArray] = ei2matrix(test);
% cmap = colormap(flipud(gray));
% testArray = testArray / 16.05148356;
% figure; imagesc(testArray); colorbar; colormap(cmap); title('Approx. # of axons - bi-electrode avg'); axis image; set(gcf, 'InvertHardCopy', 'Off');


axonVoltageByElectrode = zeros(1, 512);

%axonVoltageByElectrode(electrodeNumberingScheme(patternNos)) = bundleMaxes(patternNos);

axonVoltageByElectrode(electrodeNumberingScheme(singleElecPatterns)) = singleElecMaxes(singleElecPatterns);
%axonVoltageByElectrode(electrodeNumberingScheme(biElecPatterns1)) = biElec1Maxes(biElecPatterns1);
%axonVoltageByElectrode(electrodeNumberingScheme(biElecPatterns2)) = biElec2Maxes(biElecPatterns2);


[~, electrodeArray] = ei2matrix(axonVoltageByElectrode);
cmap = colormap(flipud(gray));
%cmap(1,:) = [0,0,0]

electrodeArray = electrodeArray / 16.05148356;

figure; imagesc(electrodeArray); colorbar; colormap(cmap); title('Approx. # of axons - bi-electrode 2'); axis image; set(gcf, 'InvertHardCopy', 'Off');


% Load matrix containing the electrode numbers for the 512-electrode MEA
%%

patternNos = singleElecPatterns;
patternNos = 254;
display = false;

totalThreshVolts = 0;

for patternNo = patternNos
    bundleMean = getBundleVoltages(patternNo, display);
    threshindices = find(round(abs(bundleMean(:, 2))*100)/100 == abs(axonBundleThresholds(patternNo)));
    disp(threshindices)
    singleneuronthresh = bundleMean(find(bundleMean(:, 1) < -16.05148356, 1, 'first'), 2);
    if singleneuronthresh == bundleMean(threshindices(1), 2)
        disp('valid');
    end
    avgVoltsForThreshAmps = bundleMean(round(abs(bundleMean(:, 2))*100)/100 == abs(axonBundleThresholds(patternNo)));
end


%% find axon voltage for neurons

neuronID = 2842;

EI = datarun.ei.eis{get_cell_indices(datarun, neuronID)};
eiamps = max(EI') - min(EI');
eiamps2 = eiamps*8;
[eiMax, maxIndex] = max(eiamps);

xpos = positions(maxIndex, 1);
ypos = positions(maxIndex, 2);
yrow = (450-ypos)/60;

axon = zeros(16-yrow, 2);


axon(1, 1) = eiamps(maxIndex);
axon(1, 2) = maxIndex;
i = 2;
startypos = ypos;
for ypos = fliplr(-450:60:(startypos-60))
    leftID = electrodeIDForXY(positions, xpos-30, ypos);
    rightID = electrodeIDForXY(positions, xpos+30, ypos);
    if leftID == 0 || (rightID && eiamps(leftID) < eiamps(rightID))
        xpos = xpos+30;
        axon(i, 1) = eiamps(rightID);
        axon(i, 2) = rightID;
    elseif (rightID == 0 && leftID) || (leftID && eiamps(rightID) < eiamps(leftID))
        xpos = xpos-30;
        axon(i, 1) = eiamps(leftID);
        axon(i, 2) = leftID;
    end
    i = i + 1;
    
    
end

axonMean = mean(axon(3:end, 1), 1);



%figure('position', [0 500 900 500]); scatter(positions(:,1),positions(:,2),eiamps2,[0.5,0.5,1],'filled'); axis image; title(['neuron ' num2str(neuronID) ': '  num2str(threshDiffDiffs(index))]); whitebg('black'); axis off; set(gcf, 'InvertHardCopy', 'off');
figure('position', [0 500 900 500]); scatter(positions(:,1),positions(:,2),350, eiamps,'filled'); colorbar; caxis([0 25]); axis image; title(['neuron ' num2str(neuronID) ': '  num2str(axonMean) ' mV average']); whitebg('black'); axis off; set(gcf, 'InvertHardCopy', 'off');
%hold on; scatter(positions(maxIndex, 1), positions(maxIndex, 2), 350, 'white', 'filled');
hold on; scatter(positions(axon(3:end,2), 1), positions(axon(3:end, 2), 2), 350, 'white', 'filled');