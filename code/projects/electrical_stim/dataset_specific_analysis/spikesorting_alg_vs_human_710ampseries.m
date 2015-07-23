%% script to generate 
load('/Volumes/Lab/Users/grosberg/Gonzalo_sorting_test_data/OutputAll.mat')
%%
patternNos = zeros(size(OutputAll)); 
thresh_hum = zeros(size(OutputAll)); 
thresh_alg = zeros(size(OutputAll)); 
sameRS = zeros(size(OutputAll)); 
distantRS = zeros(size(OutputAll)); 
neighborRS = zeros(size(OutputAll)); 

for p=1:size(OutputAll,2)
    
    alg_output = OutputAll(p);
    gm_path = alg_output.path;
    ii = find(gm_path==filesep,3,'last');
    my_path = fullfile('/Volumes/Analysis/',gm_path(ii(1):ii(3)));
    alg_output.path = my_path; 
    
    [thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(alg_output);
    patternNos(p) = alg_output.stimInfo.patternNo;
    recElec = cell2mat(alg_output.neuronInfo.prefElectrodes); % How to handle multiple recording electrodes?
    if isempty(intersect(recElec, alg_output.stimInfo.listStimElecs(:)))
        % The recording electrodes are different than the stimulating elecs
        clusterElectrodes = getCluster512(recElec); 
        if isempty(intersect(clusterElectrodes, alg_output.stimInfo.listStimElecs(:)))
            % distant electrode
            distantRS(p) = 1; 
        else
            % neighboring electrode
            neighborRS(p) = 1; 
        end
    else
        % At least one recording electrode is the same as the stimulating
        % elec
        sameRS(p) = 1; 
    end
    if (thresholdHum > 5) || (thresholdHum < 0)
         thresh_hum(p) = 100; 
    else
        thresh_hum(p) = thresholdHum;
    end
    if (thresholdAlg > 5) || (thresholdAlg < 0)
        thresh_alg(p) = 100;
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

figure; %subplot(1,3,1); 
line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_hum, thresh_alg, 50, [1  .271  0],'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title(sprintf('50%% activation thresholds, %0.0f stim patterns\n recording electrode != stimulating electrode',length(thresh_hum))); 


ratios = thresh_alg./thresh_hum; 
figure; hist(ratios,[0:.2:5])