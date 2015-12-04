%% script to generate 
load('/Volumes/Lab/Users/grosberg/Gonzalo_sorting_test_data/OutputAll.mat')
%%
patternNos = zeros(size(OutputAll)); 
thresh_hum = zeros(size(OutputAll)); 
thresh_alg = zeros(size(OutputAll)); 
sameRS = zeros(size(OutputAll)); 
distantRS = zeros(size(OutputAll)); 
neighborRS = zeros(size(OutputAll)); 
aa = 1; 
bb = 1; 
cc = 1; 
%%
for kk =  551:560  %1:size(OutputAll,2)
    alg_output = OutputAll(kk);
    gm_path = alg_output.path;
    ii = find(gm_path==filesep,3,'last');
    my_path = fullfile('/Volumes/Analysis/',gm_path(ii(1):ii(3)));
    alg_output.path = my_path; 
%     disp([my_path 'elecResp_n' num2str(alg_output.neuronInfo.neuronIds) ...
%         '_p' num2str(alg_output.stimInfo.patternNo) ]);

    compareAlgHum_dataplots(alg_output, 1);
    suptitle(num2str(kk)); 
%     recElec = alg_output.tracesInfo.recElecs; % How to handle multiple recording electrodes?
%     stimElec = unique(alg_output.stimInfo.listStimElecs(:));
%     if isempty(intersect(recElec, stimElec))
%         % The recording electrodes are different than the stimulating elecs
%         clusterElectrodes = getCluster512(recElec);
%         if isempty(intersect(clusterElectrodes, stimElec)) % distant electrode
% %             disp('distant');
%         else % neighboring electrode
% %             disp('neighbor'); 
%         end
%     else % At least one recording electrode is the same as the stimulating elec
% %         disp('same');
%         disp(num2str(kk)); 
%         compareAlgHum_dataplots(alg_output, 1);
%     end

end
%%
for p=1:size(OutputAll,2)
    
    alg_output = OutputAll(p);
    gm_path = alg_output.path;
    ii = find(gm_path==filesep,3,'last');
    my_path = fullfile('/Volumes/Analysis/',gm_path(ii(1):ii(3)));
    alg_output.path = my_path; 
    
    [thresholdHum, thresholdAlg, curveHum, ~,paramsHum,~] = fitToErfOutputAndHuman(alg_output);
     
    patternNos(p) = alg_output.stimInfo.patternNo;
    recElec = alg_output.tracesInfo.recElecs; % How to handle multiple recording electrodes?
    stimElec = unique(alg_output.stimInfo.listStimElecs(:)); 
%     disp(['rec electrode: ' num2str(recElec) '; stim electrode(s): ' num2str(stimElec')])
%     disp(num2str(length(curveHum))); 
    if isempty(intersect(recElec, stimElec))
        % The recording electrodes are different than the stimulating elecs
        clusterElectrodes = getCluster512(recElec); 
        if isempty(intersect(clusterElectrodes, stimElec))
            % distant electrode
            distantRS(p) = 1; 
            if (thresholdHum < 4) && (thresholdHum > 0)
            distantRS_curve{aa} = curveHum; 
            distantRS_params(aa,:) = paramsHum; 
            if aa>80 && aa<84
                disp([alg_output.path 'p' ...
                    num2str(alg_output.stimInfo.patternNo) ' n' ...
                    num2str(alg_output.neuronInfo.neuronIds)]); 
            end
            aa = aa + 1; 
            end
        else
            % neighboring electrode
            neighborRS(p) = 1; 
            if (thresholdHum < 4) && (thresholdHum > 0)
            neighborRS_curve{bb} = curveHum; 
            neighborRS_params(bb,:) = paramsHum; 
            bb = bb + 1; 
            end
        end
    else
        % At least one recording electrode is the same as the stimulating
        % elec
        sameRS(p) = 1;
        if (thresholdHum < 4) && (thresholdHum > 0)
        sameRS_curve{cc} = curveHum;
        sameRS_params(cc,:) = paramsHum; 
        cc = cc + 1;
        end
%         compareAlgHum_dataplots(alg_output, 1);
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
% %         Display results if the patterns are not well matched.
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

%% Plot the activation curves of the different cases
figure; 
subplot(1,3,1); 
acolors = summer(round(1.3*size(distantRS_curve,2))); 
for aa = [1:80 84:size(distantRS_curve,2)]
    curve = distantRS_curve{aa}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        plot(curve(2,idx-500:2:idx+500),'Color',acolors(aa,:));
    catch
        disp(['missed one at aa = ' num2str(aa)])
    end
end
title('distant rec & stim'); 

bcolors = winter(round(1.3*size(neighborRS_curve,2))); 
subplot(1,3,2); 
% figure; 
for bb = [1:31 33:size(neighborRS_curve,2)]
    curve = neighborRS_curve{bb}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        plot(curve(2,idx-500:2:idx+500),'Color',bcolors(bb,:));
    catch
        disp(['missed one at bb = ' num2str(bb)])
    end
    
end
title('neighbor rec & stim');
ccolors = autumn(round(size(sameRS_curve,2)*1.3));
subplot(1,3,3); 
% figure; 
for cc = 1:size(sameRS_curve,2)
    curve = sameRS_curve{cc}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        plot(curve(2,idx-500:2:idx+500),'Color',ccolors(cc,:));
    catch
        disp(['missed one at cc = ' num2str(cc)])
    end
end
title('Same rec & stim'); 

figure; 
for aa = [1:80 84:size(distantRS_curve,2)]
    curve = distantRS_curve{aa};
    hold on; 
    plot(curve(1,:),curve(2,:),'Color',acolors(aa,:));
%     disp(num2str(aa));
%     pause; 
end
title('distant rec & stim'); 

figure; 
for bb = [1:31 33:size(neighborRS_curve,2)]
    curve = neighborRS_curve{bb};
    hold on;
    plot(curve(1,:),curve(2,:),'Color',bcolors(bb,:));
end
title('neighbor rec & stim');

figure; 
for cc = 1:size(sameRS_curve,2)
    curve = sameRS_curve{cc};
    hold on; 
    plot(curve(1,:),curve(2,:),'Color',ccolors(cc,:));
end
title('same rec & stim'); 

%% Plot the different curves as overlays

figure; 
acolors = summer(round(1.3*size(distantRS_curve,2))); 
for aa = [1:80 84:size(distantRS_curve,2)]
    curve = distantRS_curve{aa}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        plot(curve(2,idx-500:2:idx+500),'Color',acolors(aa,:));
    catch
        disp(['missed one at aa = ' num2str(aa)])
    end
end
title('distant rec & stim'); 

bcolors = winter(round(1.3*size(neighborRS_curve,2))); 

for bb = [1:31 33:size(neighborRS_curve,2)]
    curve = neighborRS_curve{bb}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        plot(curve(2,idx-500:2:idx+500),'Color',bcolors(bb,:));
    catch
        disp(['missed one at bb = ' num2str(bb)])
    end
    
end
title('neighbor rec & stim');
ccolors = autumn(round(size(sameRS_curve,2)*1.3));
 
for cc = 1:size(sameRS_curve,2)
    curve = sameRS_curve{cc}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        plot(curve(2,idx-500:2:idx+500),'Color',ccolors(cc,:));
    catch
        disp(['missed one at cc = ' num2str(cc)])
    end
end
title('Same rec & stim'); 

%% Plot the shifted curves
pointsToPlot = 300; 
figure; 
subplot(1,3,3);
acolors = summer(round(1.3*size(distantRS_curve,2))); 
thresholds_distantRS = zeros(size([1:80 84:size(distantRS_curve,2)])); 
for aa = [1:80 84:size(distantRS_curve,2)]
    curve = distantRS_curve{aa}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        xvals = curve(1,:)/curve(1,idx); 
        plot(xvals,curve(2,:),'Color',acolors(aa,:));
        thresholds_distantRS(aa) = curve(1,idx);
    catch
        disp(['missed one at aa = ' num2str(aa)])
    end
     
end
title('distant rec & stim'); 
xlim([0.4 1.6]);
xlabel('stimulation amplitude (normalized \muA)'); 

% bcolors = winter(round(1.3*size(neighborRS_curve,2))); 
subplot(1,3,2);
thresholds_neighborRS = zeros(size([1:31 33 35:size(neighborRS_curve,2)])); 
for bb = [1:31 33 35:size(neighborRS_curve,2)]
    curve = neighborRS_curve{bb}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        xvals = curve(1,:)/curve(1,idx); 
        plot(xvals,curve(2,:),'Color',bcolors(bb,:));
        xlim([0.4 1.6]);
    catch
        disp(['missed one at bb = ' num2str(bb)])
    end
     thresholds_neighborRS(bb) = curve(1,idx);
end
xlim([0.4 1.6]);
xlabel('stimulation amplitude (normalized \muA)'); 
title('neighbor rec & stim');
% ccolors = autumn(round(size(sameRS_curve,2)*1.3));
thresholds_sameRS = zeros(1,size(sameRS_curve,2)); 
subplot(1,3,1);
for cc = 1:size(sameRS_curve,2)
    curve = sameRS_curve{cc}; 
    idx = find(curve(2,:)>0.5,1,'first');
    hold on; 
    try
        xvals = curve(1,:)/curve(1,idx); 
        plot(xvals,curve(2,:),'Color',ccolors(cc,:));
    catch
        disp(['missed one at cc = ' num2str(cc)])
    end
    thresholds_sameRS(cc) = curve(1,idx);
end
title('same rec & stim'); 
xlim([0.4 1.6]);
ylabel('activation probability')
xlabel('stimulation amplitude (normalized \muA)'); 
%% Plot the sigma of the cumulative gaussians for the different activation curves
distant_vals = [1:80 84:size(distantRS_curve,2)];
neighbor_vals = [1:31 33 35:size(neighborRS_curve,2)];

distantRS_sigma = 1./(distantRS_params(distant_vals,1)*sqrt(2)); 
neighborRS_sigma = 1./(neighborRS_params(neighbor_vals,1)*sqrt(2)); 
sameRS_sigma = 1./(sameRS_params(:,1)*sqrt(2)); 

types = [cellstr(repmat('distant',length(distantRS_sigma),1)); ...
    cellstr(repmat('neighbor',length(neighborRS_sigma),1));...
    cellstr(repmat('same',length(sameRS_sigma),1))];
figure; boxplot([distantRS_sigma; neighborRS_sigma; sameRS_sigma],types) ;



thresh_sameRS = -sameRS_params(:,2)./sameRS_params(:,1);
thresh_neighborRS = -neighborRS_params(neighbor_vals,2)./neighborRS_params(neighbor_vals,1);
thresh_distantRS = -distantRS_params(distant_vals,2)./distantRS_params(distant_vals,1);

distantRS_coeffVar = 1./(distantRS_params(distant_vals,1)*sqrt(2))./thresh_distantRS; 
neighborRS_coeffVar = 1./(neighborRS_params(neighbor_vals,1)*sqrt(2))./thresh_neighborRS; 
sameRS_coeffVar = 1./(sameRS_params(:,1)*sqrt(2))./thresh_sameRS; 

types = [cellstr(repmat('same',length(sameRS_sigma),1)); ...
    cellstr(repmat('neighbor',length(neighborRS_sigma),1));...
    cellstr(repmat('distant',length(distantRS_sigma),1))];
figure; boxplot([sameRS_coeffVar; neighborRS_coeffVar; distantRS_coeffVar],types) ;
title('Coefficients of variation')