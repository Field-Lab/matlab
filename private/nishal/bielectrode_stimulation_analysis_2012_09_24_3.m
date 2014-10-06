% Script to plot the stimulation thresholds for ON Parasol cells in
% 2012-09-24-3 512 stimulation data for bi-electrode stimulation
% experiments

% LGrosberg June 2014

%% Load data
currDir = pwd;
if strcmp(currDir(1:5),'/home')
    pname = '~/Development/matlab-standard/code/lab'; 
    addpath(genpath(pname));
    visionPath = '/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar';
    javaaddpath(visionPath);
end

dataPath = '/Volumes/Analysis/2012-09-24-3/data007/data007';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
positions = datarun.ei.position;


elecStimPath = '/Volumes/Analysis/2012-09-24-3/data008/';
dirInfo = dir(elecStimPath);

% Load thresholds for bidirection axon bundle activation (determined by visual inspection of data)
temp = load('~/Development/matlab-standard/private/lauren/MATLAB_code/analysis/dataset-specific/axonBundleThresholds_byPattern_2012_09_24_3_data008.mat'); 
axonBundleThresholds = temp.axonBundleThresholds_byPattern_2012_09_24_3_data008; 
%%
warning('off','stats:nlinfit:IllConditionedJacobian'); 
warning('off','stats:nlinfit:IterationLimitExceeded'); 
% ON Parasol cells in quadrant 4
neurons = [2 33 168 197 229 303 332 436 601 647 919 4909 6603 6633 6769 ...
    6936 6995 7008 7173 7218 7474 7507 7549]; %7607

activationThresh_bielec1    = zeros(size(neurons,2),512);                   
activationThresh_bielec2    = zeros(size(neurons,2),512);
activationThresh_singleElec = zeros(size(neurons,2),512);

thresholdDiff_bielec1 = zeros(size(neurons,2),4);
thresholdDiff_singleelec = zeros(size(neurons,2),4);
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
                if ii<size(stimAmps,1)
                    xx = stimAmps(ii-1):0.001:stimAmps(ii+1);
                    yy = f(F_fitted,xx);
                    
                    threshVoltage = xx(find(yy>0.5,1,'first'));
                    switch stimType
                        case 'bielectrode_P1'
                            activationThresh_bielec1(j,leftElectrode) = threshVoltage;
%                             bielecDiff(x) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            thresholdDiff_bielec1(j,a) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            a = a+1; 
%                             x = x + 1; 
                            disp(['Bielectrode hit for neuron ' num2str(neurons(j)) 'patternNo: ' num2str(elecResp.stimInfo.patternNo)]); 
                        case 'bielectrode_P2'
                            activationThresh_bielec2(j,leftElectrode) = threshVoltage;
                        case 'singleElectrode'
                            activationThresh_singleElec(j,elecResp.stimInfo.electrodes) = threshVoltage;
%                             singleelecDiff(z) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            thresholdDiff_singleelec(j,b) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            b = b + 1 ; 
%                             z = z + 1; 
                            disp(['Single electrode hit for neuron ' num2str(neurons(j)) 'patternNo: ' num2str(elecResp.stimInfo.patternNo)]); 
                    end
                    
%                     actThresholds(elecResp.stimInfo.electrodes) = threshVoltage;
%                     
                    %Plotting
                    figure; plot(stimAmps, responseProb,'.'); ylim([0 1])
                    xlabel('Current (\muA)'); ylabel('response probability');
                    title(sprintf('neuron %d; stimulating electrode %d',neuron,elecResp.stimInfo.electrodes));
                    
                    hold on; plot(stimAmps,y,'g');
                    hold on; plot(threshVoltage,yy(find(yy>0.5,1,'first')),'ro');
                    text(threshVoltage+0.1,0.5,[num2str(threshVoltage) '\muA'])
                   
   
                   
                end
                
                
                
            end
        end
    end
    
end


%% Plotting ..


groupnames = {'2';'33';'168';'197';'229';'303';'332';'436';'601';'647';'919';'4909';'6603';'6633';'6769';'6936';'6995';'7008';'7173';'7218';'7474';'7507';'7549'}; %'7607'
bw_colormap = [1 0 0];
gridstatus = 'y';
figure; 
barweb(thresholdDiff_bielec1, zeros(size(thresholdDiff_bielec1)), [], groupnames, 'Bi-electrode stimulation difference from bundle activation', 'neuron', 'soma threshold - bundle threshold (uA)', bw_colormap, gridstatus, []);

figure; 
barweb(thresholdDiff_singleelec, zeros(size(thresholdDiff_singleelec)), [], groupnames, 'Single electrode stimulation difference from bundle activation', 'neuron', 'soma threshold - bundle threshold (uA)', bw_colormap, gridstatus, []);

thresholdDiff_singleelec(thresholdDiff_singleelec==0) = NaN; 
thresholdDiff_bielec1(thresholdDiff_bielec1==0)       = NaN; 
threshDiffs = cat(2,min(thresholdDiff_singleelec,[],2),min(thresholdDiff_bielec1,[],2));
figure; 
bw_legend = {'single electrode stim','bi-electrode stim'}; 
bw_colormap = [1 0 0; 0 0 1]; 
barweb(threshDiffs, zeros(size(threshDiffs)), [], groupnames, 'Somatic activation of ON Parasols VS bundle activation', 'neuron', 'soma threshold - bundle threshold (uA)', bw_colormap, gridstatus, bw_legend);


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
neuronID = 33; 
EI = datarun.ei.eis{get_cell_indices(datarun, neuronID)}; 
eiamps = max(EI') - min(EI'); 
figure; scatter(positions(:,1),positions(:,2),eiamps); axis image; 
if printElecs
    for x = 1:512
          
        th=text(positions(x,1),positions(x,2),num2str(x),'HorizontalAlignment','center');
        set(th,'color','g');
        
 
    end
end