function values = algHumCmp(input,Gibbs, n, values)

% Compares the results of Gonzalo's spike sorting algorithm to 
% if the optional argument is one, then the x axis show amplitude condition
% instead of stimulus amplitude
neuronIds = input.neuronInfo.neuronIds;
neuronId = neuronIds(n);

latencies = cell2mat(Gibbs.variables.latencies) ;

spikes = cell2mat(Gibbs.variables.spikes);

% Load elecResp file
pathname = input.names.path; 
fname = ['elecResp_n' num2str(neuronId) '_p' num2str(input.stimInfo.patternNo) '.mat']; 
filename = fullfile(pathname,fname); 
temp = load(filename); 
elecResp = temp.elecResp; 
disp(length(elecResp.stimInfo.stimAmps))
if length(elecResp.stimInfo.stimAmps)>40
    keyboard;
end
for a = 1 : size(elecResp.analysis.latencies,1)
    try
        humanLat(a,:) = elecResp.analysis.latencies{a}; %#ok<AGROW>
    catch
        humanLat(a,1:(end-1)) = elecResp.analysis.latencies{a}; %#ok<AGROW>
        humanLat(a,end) = NaN; 
    end
end
if size(latencies,2) == (size(humanLat,2) - 1)
    humanLat(:,1) = [];
end
humanSpikes = humanLat>0; 

% figure(1000); cla;  
% subplot(1,2,1); imagesc(spikes); title('Algorithm spikes')
% subplot(1,2,2); imagesc(humanSpikes); title('Human spikes'); 

agreement = (humanSpikes == spikes); 
totalTrials = numel(humanSpikes);
percentAgreement = sum(agreement(:))/totalTrials;
 
totalHumanNegatives = sum(sum(humanSpikes == 0));
numFalsePositives = sum(spikes(humanSpikes == 0));
numTrueNegatives = sum(spikes(humanSpikes == 0) == 0);
percentFalsePositive = numFalsePositives/totalHumanNegatives;
percentTrueNegative = numTrueNegatives/totalHumanNegatives;

totalHumanPositives = sum(sum(humanSpikes == 1));
numTruePositives = sum(spikes(humanSpikes == 1));
numFalseNegatives = sum(spikes(humanSpikes == 1) == 0);
percentTruePositive = numTruePositives/totalHumanPositives;
percentFalseNegatives = numFalseNegatives/totalHumanPositives;

values = nansum([values  [numTruePositives; numFalseNegatives; totalHumanPositives; ...
    numFalsePositives; numTrueNegatives; totalHumanNegatives;...
    sum(agreement(:)); totalTrials]],2); 

% 
% clear latencies
% for j=1:J
%     lats = Gibbs.variables.latencies{n}(j,1:I(j));
%     lats = lats(lats>0);
%     if(isempty(lats))
%         latencies(j,1:3)=NaN;
%     else
%         
%         latencies(j,1) = prctile(lats,25);
%         latencies(j,2) = prctile(lats,50);
%         latencies(j,3) = prctile(lats,75);
%     end
% end
% clear sigma
% sigma  = Gibbs.variables.sigma(e,:);
% 
% subplot(2,2,1)
% plot(amps,spikeProbs,'.-','linewidth',3); 
% grid('on'); 
% hold on; plot(amps,spikeLogProbs,'--','linewidth',2); 
% xlabel('stimulation amplitude (\muA)'); 
% ylabel('spike probability'); 
% title(sprintf('activation curve for n%0.0f',neuronId)); 
% 
% for b = breakAxon
%     plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','green')
% end
% for b = breakRecElec
%     plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','red')
% end
% axis([amps(1) amps(end) 0 1])
% 
% 
% subplot(2,2,2);
% dots = plot(amps,latencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
% hold on
% grid('on')
% crosses = plot(amps,latencies(:,3),'x','markersize',15); 
% plot(amps,latencies(:,1),'x','markersize',15);
% xlabel('stimulation amplitude (\muA)'); 
% ylabel('latencies (sample points)'); 
% legend([dots,crosses],'medians','quartiles'); % of responding trials
% title(sprintf('latencies for n%0.0f',neuronId)); 
% 
% for b = breakAxon
%     
%     plot([amps(b) amps(b)], [nanmin(nanmin(latencies))-3 nanmax(nanmax(latencies))+3],'linewidth',1.5,'color','green')
% end
% for b = breakRecElec
%     plot([amps(b) amps(b)], [nanmin(nanmin(latencies))-3 nanmax(nanmax(latencies))+3],'linewidth',1.5,'color','red')
% end
% grid('on')
% axis([amps(1) amps(J) nanmin(nanmin(latencies))-3 nanmax(nanmax(latencies))+3])
% 
% 
% 
% subplot(2,2,3); 
% hum = plot(amps,elecResp.analysis.successRates,...
%     'o-','LineWidth',3); grid on; 
% hold on; alg = plot(amps,spikeProbs,'^-','linewidth',3); 
% legend([hum alg],'human','algorithm'); 
% xlabel('stimulation amplitude (\muA)'); 
% title('human analysis results'); 
% ylim([0 1]); 
% 
% % Get latency information from the elecResp files. 
% for ii = 1:size(elecResp.analysis.latencies,1)
%     humanLats = elecResp.analysis.latencies{ii}; 
%     humanLats = humanLats(humanLats>0); 
%     humanLatencies(ii,1) = prctile(humanLats,25);
%     humanLatencies(ii,2) = prctile(humanLats,50);
%     humanLatencies(ii,3) = prctile(humanLats,75);
% end
% subplot(2,2,4); 
% dots = plot(amps,humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
% hold on; grid on; 
% crosses = plot(amps,humanLatencies(:,3),'x','markersize',15); 
% plot(amps,humanLatencies(:,1),'x','markersize',15);
% xlabel('stimulation amplitude (\muA)'); 
% ylabel('latencies (sample points)'); 
% legend([dots,crosses],'medians','quartiles'); % of responding trials
% title(sprintf('human analyzed latencies for n%0.0f',neuronId)); 