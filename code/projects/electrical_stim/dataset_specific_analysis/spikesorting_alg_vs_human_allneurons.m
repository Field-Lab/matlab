% This script compares the spatial activation sensitivy for single 
% electrode stimulation calculated using Gonzalo'salgorithm and human 
% sorting results. 

%% Load data
diffRecStim = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
    'validation/2015-06-algorithm-results/sorted_diffStimRecElecs']);
sameRecStim = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
    'validation/2015-06-algorithm-results/sorted_sameStimRecElecs']); 

sorted_diffStimRecElecs = diffRecStim.sorted_diffStimRecElecs; clear diffRecStim; 
sorted_sameStimRecElecs = sameRecStim.sorted_sameStimRecElecs; clear sameRecStim; 

%% Compare activation thresholds found by the human and the algorithm

patternNosD = zeros(size(sorted_diffStimRecElecs)); 
thresh_humD = zeros(size(sorted_diffStimRecElecs)); 
thresh_algD = zeros(size(sorted_diffStimRecElecs)); 

for p=1:size(sorted_diffStimRecElecs,2)
    
    alg_output = sorted_diffStimRecElecs(p);
    [thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(alg_output);
    patternNosD(p) = alg_output.stimInfo.patternNo;
    if (thresholdHum > 4) || (thresholdHum < 0)
         thresh_humD(p) = 0; 
    else
        thresh_humD(p) = thresholdHum;
    end
    if (thresholdAlg > 4) || (thresholdAlg < 0)
        thresh_algD(p) = 0;
    else
        thresh_algD(p) = thresholdAlg;
    end
    if p<100 && p>50
        % Display results if the patterns are not well matched.
        if abs(thresh_humD(p) - thresh_algD(p)) > 0.1
            displayResults_compareOutput(alg_output,1,0)
            suptitle([alg_output.path ' n:' num2str(alg_output.neuronInfo.neuronIds) ...
                ' pattern: ' num2str(alg_output.stimInfo.patternNo)]);
        end
    end
end



%% Compare activation thresholds found by the human and the algorithm

patternNos = zeros(size(sorted_sameStimRecElecs)); 
thresh_hum = zeros(size(sorted_sameStimRecElecs)); 
thresh_alg = zeros(size(sorted_sameStimRecElecs)); 

for p=1:size(sorted_sameStimRecElecs,2)
    
    alg_output = sorted_sameStimRecElecs(p);
    [thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(alg_output);
    patternNos(p) = alg_output.stimInfo.patternNo;
    if (thresholdHum > 4) || (thresholdHum < 0)
         thresh_hum(p) = 0; 
    else
        thresh_hum(p) = thresholdHum;
    end
    if (thresholdAlg > 4) || (thresholdAlg < 0)
        thresh_alg(p) = 0;
    else
        thresh_alg(p) = thresholdAlg;
    end
    if p<50
        % Display results if the patterns are not well matched.
        if abs(thresh_hum(p) - thresh_alg(p)) > 0.1
            displayResults_compareOutput(alg_output,1,0)
            suptitle([alg_output.path ' n:' num2str(alg_output.neuronInfo.neuronIds) ...
                ' pattern: ' num2str(alg_output.stimInfo.patternNo)]);
        end
    end
end

%% Scatter plot to compare the thresholds found using the two methods

figure; subplot(1,3,1); 
line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_humD, thresh_algD, 50, [1  .271  0],'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title(sprintf('50%% activation thresholds, %0.0f stim patterns\n recording electrode != stimulating electrode',length(thresh_humD))); 
subplot(1,3,2); line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_hum, thresh_alg, 50,[0  .231  1], 'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title(sprintf('50%% activation thresholds, %0.0f stim patterns\n recording electrode = stimulating electrode',length(thresh_hum))); 


% Calculate percent accuracy
idx_fnD = find(thresh_humD - thresh_algD > 0.1);
idx_fpD = find(thresh_algD - thresh_humD > 0.1);
false_pos = length(find(thresh_humD(idx_fpD)==0));
false_neg = length(find(thresh_algD(idx_fnD)==0));
total_trials = length(thresh_humD);
correct_trials = length(find(abs(thresh_humD - thresh_algD)<0.1));
false_negD = false_neg/total_trials; 
false_posD = false_pos/total_trials; 
agreementD = correct_trials/total_trials; 

idx_fnS = find(thresh_hum - thresh_alg > 0.1);
idx_fpS = find(thresh_alg - thresh_hum > 0.1);
false_pos = length(find(thresh_hum(idx_fpS)==0));
false_neg = length(find(thresh_alg(idx_fnS)==0));
total_trials = length(thresh_hum);
correct_trials = length(find(abs(thresh_hum - thresh_alg)<0.1));
false_negS = false_neg/total_trials; 
false_posS = false_pos/total_trials; 
agreementS = correct_trials/total_trials; 

all_valsD = [agreementD; false_negD; false_posD ];
all_valsS = [agreementS; false_negS; false_posS ];
errorsD = [1.96*sqrt(1/total_trials*agreementD*(1-agreementD));...
    1.96*sqrt(1/total_trials*false_negD*(1-false_negD));...
    1.96*sqrt(1/total_trials*false_negD*(1-false_negD))];
errorsS = [1.96*sqrt(1/total_trials*agreementS*(1-agreementS));...
    1.96*sqrt(1/total_trials*false_negS*(1-false_negS));...
    1.96*sqrt(1/total_trials*false_negS*(1-false_negS))];

subplot(1,3,3); 
groupnames = {'agreement';'false negative';'false positive'};
bw_colormap = [1  .271  0;
    0  .231  1];
gridstatus = 'y';
bw_legend = {'rec elec != stim elec','rec elec = stim elec'}; 
barweb([all_valsD all_valsS], [errorsD errorsS], [], ...
    groupnames, 'Threshold analysis of algorithm v. human', ...
    'error type', 'proportional agreement', bw_colormap, gridstatus, bw_legend);