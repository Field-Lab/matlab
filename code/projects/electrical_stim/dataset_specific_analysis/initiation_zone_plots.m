
% load('/Volumes/Analysis/2012-09-24-3/data003-autosort/elecRespAuto_n2796_p187.mat')
% compareAlgHum_dataplots(elecRespAuto, 	1)
% [rawData, amplitudes] = generateEiFromStimPattern(pathToAnalysisData, patternNo,varargin)

% Load elecResp file
load('/Volumes/Analysis/2015-04-09-2/data003/elecResp_n4984_p2.mat');
% Get latency information from the elecResp files. 
for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
end

humanLatencies= humanLatencies/elecResp.details.sample_rate*1000;
figure; 
dots = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,1),'x','markersize',15);
xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (ms)'); 
legend([dots,crosses],'medians','quartiles'); % of responding trials
title(sprintf('Human analyzed latencies for n%0.0f',elecResp.cells.main)); 

%% modified plots for specific dataset
relevantMovNos = 45:2:51; 
clear humanLatencies
for ii = 1:length(relevantMovNos)
    humanLats = [elecResp.analysis.latencies{relevantMovNos(ii)}; ...
        elecResp.analysis.latencies{relevantMovNos(ii)+1}]; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
    humanLatencies(ii,4) = mean(humanLats);
end

figure; 
dots = plot(abs(elecResp.stimInfo.stimAmps(relevantMovNos)),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps(relevantMovNos)),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps(relevantMovNos)),humanLatencies(:,1),'x','markersize',15);
dias = plot(abs(elecResp.stimInfo.stimAmps(relevantMovNos)),humanLatencies(:,4),'d','MarkerFaceColor','red','Markersize',6); 
xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (sample points)'); 
legend([dots,crosses,dias],'medians','quartiles','means'); % of responding trials
title(sprintf('Human analyzed latencies for n%0.0f',elecResp.cells.main)); 


%% Load elecResp file #2
load('/Volumes/Analysis/2012-09-24-3/data006/elecResp_n5748_p441.mat');
clear humanLatencies; 
% Get latency information from the elecResp files. 
for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
    humanLatencies(ii,4) = mean(humanLats); 
end

humanLatencies= humanLatencies/elecResp.details.sample_rate*1000;
figure; 
dots = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,1),'x','markersize',15);
dias = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,4),'d','MarkerFaceColor','red','Markersize',6); 

xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (ms)'); 
legend([dots,crosses,dias],'medians','quartiles','means'); % of responding trials
title(sprintf('Human analyzed latencies for n%0.0f',elecResp.cells.main)); 


%% Load elecResp file #3
load('/Volumes/Analysis/2014-09-10-0/data003/elecResp_n4998_p334.mat');
clear humanLatencies; 
% Get latency information from the elecResp files. 
for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
    humanLatencies(ii,4) = mean(humanLats); 
end

humanLatencies= humanLatencies/elecResp.details.sample_rate*1000;
figure; 
dots = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,1),'x','markersize',15);
dias = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,4),'d','MarkerFaceColor','red','Markersize',6); 

xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (ms)'); 
legend([dots,crosses,dias],'medians','quartiles','means'); % of responding trials
title(sprintf('Human analyzed latencies for n%0.0f',elecResp.cells.main)); 


%% Load elecResp file #4
load('/Volumes/Analysis/2015-04-14-0/data001/elecResp_n3947_p264.mat');
clear humanLatencies; 
% Get latency information from the elecResp files. 
for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
    humanLatencies(ii,4) = mean(humanLats); 
end

humanLatencies= humanLatencies/elecResp.details.sample_rate*1000;
figure; 
dots = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,1),'x','markersize',15);
dias = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,4),'d','MarkerFaceColor','red','Markersize',6); 

xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (ms)'); 
legend([dots,crosses,dias],'medians','quartiles','means'); % of responding trials
title(sprintf('Human analyzed latencies for n%0.0f',elecResp.cells.main)); 

%% Display electrical stim eis. 
'/Volumes/Analysis/2015-04-09-2/data003/elecResp_n4984_p2.mat'
'/Volumes/Analysis/2012-09-24-3/data006/elecResp_n5748_p441.mat'
'/Volumes/Analysis/2014-09-10-0/data003/elecResp_n4998_p334.mat'
'/Volumes/Analysis/2015-04-14-0/data001/elecResp_n3947_p264.mat'
generateEiFromStimPattern('//Volumes/Analysis/2015-04-09-2/data003/', 2,'movieNo',443)
generateEiFromStimPattern('/Volumes/Analysis/2012-09-24-3/data006/', 441,'movieNo',443)
generateEiFromStimPattern('/Volumes/Analysis/2012-09-24-3/data006/', 334,'movieNo',443)
generateEiFromStimPattern('/Volumes/Analysis/2012-09-24-3/data006/', 264,'movieNo',443)
