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
    recElec = alg_output.tracesInfo.recElecs; % How to handle multiple recording electrodes?
    stimElec = unique(alg_output.stimInfo.listStimElecs(:)); 
    disp(['rec electrode: ' num2str(recElec) '; stim electrode(s): ' num2str(stimElec')])

    if isempty(intersect(recElec, stimElec))
        % The recording electrodes are different than the stimulating elecs
        clusterElectrodes = getCluster512(recElec); 
        if isempty(intersect(clusterElectrodes, stimElec))
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
        compareAlgHum_dataplots(alg_output, 1);
%         if abs(thresholdHum - thresholdAlg) > 0.2
%         else
%             compareAlgHum_dataplots(alg_output, 1);
%         end
    end
    if (thresholdHum > 5) || (thresholdHum < 0)
         thresh_hum(p) = 50; 
    else
        thresh_hum(p) = thresholdHum;
    end
    if (thresholdAlg > 5) || (thresholdAlg < 0)
        thresh_alg(p) = 50;
    else
        thresh_alg(p) = thresholdAlg;
    end
%     if p<100
%         Display results if the patterns are not well matched.
%         if abs(thresh_hum(p) - thresh_alg(p)) > 0.1
%             displayResults_compareOutput(alg_output,1,0)
%             suptitle([alg_output.path ' n:' num2str(alg_output.neuronInfo.neuronIds) ...
%                 ' pattern: ' num2str(alg_output.stimInfo.patternNo)]);
%         end
%     end
end

%% Scatter plot to compare the thresholds found using the two methods

figure; %subplot(1,3,1); 
line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_hum, thresh_alg, 50, [1  .271  0],'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title(sprintf('50%% activation thresholds, %0.0f stim patterns',length(thresh_hum))); 

%% translated histogram
points = [thresh_hum; thresh_alg];
dot(thresh_hum,thresh_alg)
linecoords = ones(size(points)); 
linecoords(2,:) = -1*linecoords(2,:);
proj = dot(points, linecoords)./dot(linecoords,linecoords)*sqrt(2); 
figure; plot(proj);
x =[-35 -4:0.25:4 35]; 
figure; hist(proj,x)
figure; hist(proj,x); ylim([0 28]); xlim([-3 3]);
%%

ratios = thresh_alg./thresh_hum; 
figure; hist(ratios,[0:.2:5])

% Same stimulating and recording electrode
same_alg_vals = thresh_alg(find(sameRS)); 

pos_idx = find(thresh_hum(find(sameRS)) < 5); % Indices where human found a response
num_pos = length(find(thresh_hum(find(sameRS)) < 5)); % Number of trials where human found a response
num_pos_alg = length(find(same_alg_vals(pos_idx) < 5)); % Number of true positives for the algorithm
fprintf(['for the same stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were true positives\n'],num_pos_alg,num_pos,num_pos_alg/num_pos*100);

num_neg = length(find(thresh_hum(find(sameRS)) > 5)); % Number of trials where human found NO response
neg_idx = find(thresh_hum(find(sameRS)) > 5); % Indices where human found NO response
num_neg_alg = length(find(same_alg_vals(neg_idx) > 5));
fprintf(['for the same stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were true negatives\n'],num_neg_alg,num_neg,num_neg_alg/num_neg*100);
fprintf(['for the same stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were correct\n'], num_neg_alg + num_pos_alg,...
    length(same_alg_vals),(num_neg_alg + num_pos_alg)/length(same_alg_vals)*100);


num_false_neg = length(find(same_alg_vals(pos_idx) > 5)); % Number of false positives for the algorithm
fprintf(['for the same stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were false negatives\n'],num_false_neg,num_pos,...
    num_false_neg/num_pos*100);

num_false_pos = length(find(same_alg_vals(neg_idx) < 5));
fprintf(['for the same stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were false positives\n'],num_false_pos, num_neg,...
    num_false_pos/num_neg*100);
fprintf(['for the same stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were mistakes\n'], num_false_pos + num_false_neg,...
    length(same_alg_vals),(num_false_pos + num_false_neg)/length(same_alg_vals)*100);

%%
% neighboring stimulating and recording electrode
neigh_alg_vals = thresh_alg(find(neighborRS)); 

pos_idx = find(thresh_hum(find(neighborRS)) < 5); % Indices where human found a response
num_pos = length(find(thresh_hum(find(neighborRS)) < 5)); % Number of trials where human found a response
num_pos_alg = length(find(neigh_alg_vals(pos_idx) < 5)); % Number of true positives for the algorithm
fprintf(['for neighboring stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were true positives\n'],num_pos_alg,num_pos,num_pos_alg/num_pos*100);

num_neg = length(find(thresh_hum(find(neighborRS)) > 5)); % Number of trials where human found NO response
neg_idx = find(thresh_hum(find(neighborRS)) > 5); % Indices where human found NO response
num_neg_alg = length(find(neigh_alg_vals(neg_idx) > 5));
fprintf(['for  neighboring stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were true negatives\n'],num_neg_alg,num_neg,num_neg_alg/num_neg*100);
fprintf(['for neighboring stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were correct\n'], num_neg_alg + num_pos_alg,...
    length(neigh_alg_vals),(num_neg_alg + num_pos_alg)/length(neigh_alg_vals)*100);


num_false_neg = length(find(neigh_alg_vals(pos_idx) > 5)); % Number of false positives for the algorithm
fprintf(['for neighboring stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were false negatives\n'],num_false_neg,num_pos,...
    num_false_neg/num_pos*100);

num_false_pos = length(find(neigh_alg_vals(neg_idx) < 5));
fprintf(['for neighboring stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were false positives\n'],num_false_pos, num_neg,...
    num_false_pos/num_neg*100);
fprintf(['for neighboring stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were mistakes\n'], num_false_pos + num_false_neg,...
    length(neigh_alg_vals),(num_false_pos + num_false_neg)/length(neigh_alg_vals)*100);

%%  distant stimulating and recording electrode
dist_alg_vals = thresh_alg(find(distantRS)); 

pos_idx = find(thresh_hum(find(distantRS)) < 5); % Indices where human found a response
num_pos = length(find(thresh_hum(find(distantRS)) < 5)); % Number of trials where human found a response
num_pos_alg = length(find(dist_alg_vals(pos_idx) < 5)); % Number of true positives for the algorithm
fprintf(['for distant stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were true positives\n'],num_pos_alg,num_pos,num_pos_alg/num_pos*100);

num_neg = length(find(thresh_hum(find(distantRS)) > 5)); % Number of trials where human found NO response
neg_idx = find(thresh_hum(find(distantRS)) > 5); % Indices where human found NO response
num_neg_alg = length(find(dist_alg_vals(neg_idx) > 5));
fprintf(['for distant stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were true negatives\n'],num_neg_alg,num_neg,num_neg_alg/num_neg*100);
fprintf(['for distant stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were correct\n'], num_neg_alg + num_pos_alg,...
    length(dist_alg_vals),(num_neg_alg + num_pos_alg)/length(dist_alg_vals)*100);


num_false_neg = length(find(dist_alg_vals(pos_idx) > 5)); % Number of false positives for the algorithm
fprintf(['for distant stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were false negatives\n'],num_false_neg,num_pos,...
    num_false_neg/num_pos*100);

num_false_pos = length(find(dist_alg_vals(neg_idx) < 5));
fprintf(['for distant stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were false positives\n'],num_false_pos, num_neg,...
    num_false_pos/num_neg*100);
fprintf(['for distant stimulating and recording electrode, %0.0f / '...
    '%0.0f (%0.1f%%) were mistakes\n'], num_false_pos + num_false_neg,...
    length(dist_alg_vals),(num_false_pos + num_false_neg)/length(dist_alg_vals)*100);