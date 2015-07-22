% This script compares the spatial activation sensitivy of n1699 from
% 2015-04-14-0 for single electrode stimulation calculated using Gonzalo's
% algorithm and human sorting results. 

%% Load data
diffRecStim = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
    'validation/2015-06-algorithm-results/sorted_diffStimRecElecs']);
sameRecStim = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
    'validation/2015-06-algorithm-results/sorted_sameStimRecElecs']); 

sorted_diffStimRecElecs = diffRecStim.sorted_diffStimRecElecs; clear diffRecStim; 
sorted_sameStimRecElecs = sameRecStim.sorted_sameStimRecElecs; clear sameRecStim; 

%% Compare activation thresholds found by the human and the algorithm
idx_n1699=[29 477:505 556:561];
patternNos = zeros(size(idx_n1699)); 
thresh_hum_n1699 = zeros(size(idx_n1699)); 
thresh_alg_n1699 = zeros(size(idx_n1699)); 

for p=1:size(idx_n1699,2)
    if idx_n1699(p) == 29
        alg_output = sorted_sameStimRecElecs(idx_n1699(p));
    else
        alg_output = sorted_diffStimRecElecs(idx_n1699(p));
    end
    [thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(alg_output);
    patternNos(p) = alg_output.stimInfo.patternNo;
    if (thresholdHum > 4) || (thresholdHum < 0)
         thresh_hum_n1699(p) = 0; 
    else
        thresh_hum_n1699(p) = thresholdHum;
    end
    if (thresholdAlg > 4) || (thresholdAlg < 0)
        thresh_alg_n1699(p) = 0;
    else
        thresh_alg_n1699(p) = thresholdAlg; 
    end
   
    % Display results if the patterns are not well matched. 
    if abs(thresh_hum_n1699(p) - thresh_alg_n1699(p)) > 0.1
        displayResults_compareOutput(alg_output,1,0)
        suptitle([alg_output.path ' n:' num2str(alg_output.neuronInfo.neuronIds) ...
            ' pattern: ' num2str(alg_output.stimInfo.patternNo)]);
    end
end
%% Scatter plot to compare the thresholds found using the two methods
allValsH = zeros(512,1); 
allValsH(patternNos) = thresh_hum_n1699; 
allValsA = zeros(512,1); 
allValsA(patternNos) = thresh_alg_n1699; 
humanVals = thresh_hum_n1699;
algorithmVals = thresh_alg_n1699;
algorithmVals(algorithmVals>4) = 0; 
algorithmVals(algorithmVals<0) = 0; 
humanVals(humanVals>4) = 0; 
humanVals(humanVals<0) = 0; 

figure; line(0:4,0:4,'Color','k'); hold on; 
scatter(humanVals, algorithmVals,100,'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title('50% activation thresholds for n1699'); 
%% Plot the EI contour with a linear fit + electrodes
[xc,yc] = getElectrodeCoords512(); 
codebase_path = matlab_code_path; 
cmap = load([codebase_path 'code/projects/electrical_stim/'...
    'resources/redtealcmap.mat']); 

% Load EI info from an elecResp file
temp = load('/Volumes/Analysis/2015-04-14-0/data001/elecResp_n1699_p410.mat');
ei = temp.elecResp.cells.mainEI;
neuronId = temp.elecResp.cells.main; clear temp; 

% Show the EI as a contour plot for algorithm output display
eiAmps = max(ei,[],2) - min(ei,[],2);
eiContour_wLinFit(eiAmps,'linFitThresh',5)
title(sprintf('algorithm sorting output for neuron %d',neuronId)); 
% Show stimulating electrodes for the algorithm output
allValsA(allValsA>4) = 0; 
allValsA(allValsA<0) = 0; 
idx = find(allValsA); 
hold on; scatter(xc(idx),yc(idx),150, allValsA(idx),'filled'); 
colormap(cmap.m); caxis([0.5 3.5]); colorbar; 

% Show the EI as a contour plot for human output display
eiContour_wLinFit(eiAmps,'linFitThresh',5)
title(sprintf('human sorting output for neuron %d',neuronId)); 
% Show stimulating electrodes for the human ouput
allValsH(allValsH>4) = 0; 
allValsH(allValsH<0) = 0; 
idx = find(allValsH); 
hold on; scatter(xc(idx),yc(idx),150, allValsH(idx),'filled'); 
colormap(cmap.m); caxis([0.5 3.5]); colorbar; 