dataPath = '/Volumes/Analysis/2012-09-24-3/data007/data007';%'/Volumes/Analysis/2012-09-24-3/data007/data007';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
positions = datarun.ei.position;


elecStimPath = '/Volumes/Analysis/2012-09-24-3/data008/'; %'/Volumes/Analysis/2012-09-24-3/data008/';
dirInfo = dir(elecStimPath);

elecsByPattern = full(getElecsFromPatternFiles([elecStimPath 'pattern_files/'])); %load the matrix that gives electrode numbers indexed by pattern

%%

currentThresh = 4;

warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:IterationLimitExceeded');

on_parasol = [2 33 168 197 229 303 332 436 601 647 919 4909 6603 6633 6769 ...
    6936 6995 7008 7173 7218 7474 7507 7549];

off_parasol = [212 421 903 5568 6605 7052 7639 7532];

neurons = [on_parasol off_parasol];

somaStimThreshs_bielec1_on = zeros(size(on_parasol,2),6);
somaStimThreshs_singleelec_on = zeros(size(on_parasol,2),6);

activationThresh_bielec1_on    = zeros(size(on_parasol,2),512);
activationThresh_bielec2_on    = zeros(size(on_parasol,2),512);
activationThresh_singleElec_on = zeros(size(on_parasol,2),512);

responseProb_on = zeros(size(on_parasol,2), 6);

somaStimThreshs_bielec1_off = zeros(size(off_parasol,2),6);
somaStimThreshs_singleelec_off = zeros(size(off_parasol,2),6);

activationThresh_bielec1_off    = zeros(size(off_parasol,2),512);
activationThresh_bielec2_off    = zeros(size(off_parasol,2),512);
activationThresh_singleElec_off = zeros(size(off_parasol,2),512);

responseProb_off = zeros(size(off_parasol,2), 6);

responseCurves = zeros(size(neurons,2),6, 66);
activationThresh_bielec1 = zeros(size(neurons,2),512);
activationThresh_bielec2 = zeros(size(off_parasol,2),512);
activationThresh_singleElec = zeros(size(neurons,2),512);
somaStimThreshs_bielec1 = zeros(size(neurons,2),6);
somaStimThreshs_singleElec = zeros(size(neurons,2),6);
responseProbs = zeros(size(neurons,2), 6);

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
                
                currentIndex = find(stimAmps>currentThresh, 1,'first');
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
                            responseProbs(j,a) = y(currentIndex);
                            if size(y, 1) < 60
                                y = kron(y, ones([2 1]));
                            end
                            responseCurves(j,a,1:size(y,1)) = y;
                            somaStimThreshs_bielec1(j,a) = threshVoltage;
                            a = a+1;
                            
                            %                             x = x + 1;
                            disp(['Bielectrode hit for neuron ' num2str(neurons(j)) ' patternNo: ' num2str(elecResp.stimInfo.patternNo)]);
                        case 'bielectrode_P2'
                            activationThresh_bielec2(j,leftElectrode) = threshVoltage;
                        case 'singleElectrode'
                            activationThresh_singleElec(j,elecResp.stimInfo.electrodes) = threshVoltage;
                            %                             singleelecDiff(z) = threshVoltage-axonBundleThresholds(elecResp.stimInfo.patternNo);
                            somaStimThreshs_singleElec(j, b) = threshVoltage;
                            b = b + 1 ;
                            %                             z = z + 1;
                            disp(['Single electrode hit for neuron ' num2str(neurons(j)) ' patternNo: ' num2str(elecResp.stimInfo.patternNo)]);
                    end
                end
            end
        end
    end
end



%%

act_max = 1;
act_min = 0;

tmp = figure;
rgbVals = colormap(tmp,gray);
close(tmp);

currentIndex = 40;

% Overlay rf fits and IDs
fig = figure('KeyPressFcn', @figKeyPress); set(gcf,'Color',[1 1 1]); hold on;
data = guidata(fig);
data.activationThresh = activationThresh_singleElec;
data.currentIndex = currentIndex;
data.responseProbs = responseCurves(:,:,currentIndex);
data.responseCurves = responseCurves;
data.off_parasol = off_parasol;
data.on_parasol = on_parasol;
data.neurons = neurons;
data.datarun = datarun;
guidata(fig, data);
for n = 1:1:size(neurons,2)
    currentProb = responseProbs(n,:);
    [~,~,val] = find(currentProb);
    fit_color = 'k';
    if ~isempty(val)
        rgbIndex = max(val);
        rgbIndex = round(rgbIndex * size(rgbVals,1));
        if rgbIndex<1; rgbIndex = 1; elseif rgbIndex>size(rgbVals,1); rgbIndex = size(rgbVals,1); end
        overlaps = [];
        if ismember(neurons(n), off_parasol)
            overlaps = getOverlappingNeurons(datarun, neurons(n), neurons);
            overlaps = overlaps(overlaps ~= 0);
        end
        if ~isempty(overlaps)
            for overlap = overlaps
                index = find(neurons == overlap);
                if abs(max(find(responseProbs(index,:))) - max(val)) < 0.1
                  fit_color = 'r';
                end
            end
        end
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', fit_color,'fill_color',rgbVals(rgbIndex,:));
    else
        plot_rf_summaries(datarun, neurons(n), 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', fit_color,'fill_color','k');
    end
    
end
axis image; axis off;
colorbar; colormap(rgbVals); caxis([act_min act_max]);
title('parasol response prob, bi-electrode stimulation');



%%

overlaps = getOverlappingNeurons(datarun, neurons);
