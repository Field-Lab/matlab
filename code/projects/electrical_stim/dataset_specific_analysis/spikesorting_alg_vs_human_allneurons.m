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

patternNos = zeros(size(sorted_diffStimRecElecs)); 
thresh_hum = zeros(size(sorted_diffStimRecElecs)); 
thresh_alg = zeros(size(sorted_diffStimRecElecs)); 

for p=1:size(sorted_diffStimRecElecs,2)
    
    alg_output = sorted_diffStimRecElecs(p);
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

figure; line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_hum, thresh_alg, 100,'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title('50% activation thresholds for all axonal activation trials'); 


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

figure; line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_hum, thresh_alg, 100,'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title('50% activation thresholds for all somatic activation trials'); 