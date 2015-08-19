function displayResults_compareOutput(Output,n,varargin)
% Displays the result of Gonzalo's spike sorting and compares it to the
% result from human analysis. There must exist an elecResp file in order to
% run this function
% if the optional argument is one, then the x axis show amplitude condition
% instead of stimulus amplitude



figure; set(gcf,'Position', [387         221        1521         877]); 

neuronIds = Output.neuronInfo.neuronIds;
neuronId = neuronIds(n);
e = Output.neuronInfo.prefElectrodes{n}(1);
recElecs     = Output.tracesInfo.recElecs;
breakAxon    =  Output.tracesInfo.breakAxon{e};
breakRecElec = Output.tracesInfo.breakRecElecs{e};
colors       = colormap(jet);
breakPointsAll = sort(unique([breakAxon breakRecElec]));
J = Output.tracesInfo.J;
I = Output.tracesInfo.I;
T = Output.tracesInfo.T;

spikeProbs    = nansum(Output.spikes{n}')./I;
spikeLogProbs = Output.LogisticReg(n,:);
% Define optional variable
if(nargin == 3)
    if(varargin{1} == 0)
        amps = abs(Output.stimInfo.listAmps(:,1))';
    else
        amps = [1:J];
    end
else
     amps=abs(Output.stimInfo.listAmps(:,1))';
end
clear latencies
for j=1:J
    lats = Output.latencies{n}(j,1:I(j));
    lats = lats(lats>0);
    if(isempty(lats))
        latencies(j,1:3)=NaN;
    else
        
        latencies(j,1) = prctile(lats,25);
        latencies(j,2) = prctile(lats,50);
        latencies(j,3) = prctile(lats,75);
    end
end

latencies = latencies/Output.params.sampRate*1000;

pathname = Output.path; 
fname = ['elecResp_n' num2str(neuronId) '_p' num2str(Output.stimInfo.patternNo) '.mat']; 
filename = fullfile(pathname,fname); 
temp = load(filename); 
elecResp = temp.elecResp; 

subplot(2,3,1); 
hum = plot(abs(elecResp.stimInfo.stimAmps),elecResp.analysis.successRates,...
    'o-','LineWidth',3); grid on; 
hold on; alg = plot(amps,spikeProbs,'^-','linewidth',3); 
logRegAlg = plot(amps,spikeLogProbs,'--','linewidth',2); 
legend([hum alg(1,1) logRegAlg],'Human','Algorithm','Logistic Regression (algorithm)'); 
xlabel('stimulation amplitude (\muA)'); 
title(sprintf('Activation curves for n%0.0f',neuronId')); 
ylim([0 1]); 

xlabel('stimulation amplitude (\muA)'); 
ylabel('spike probability'); 
%title(sprintf('activation curve for n%0.0f',neuronId)); 

for b = breakAxon
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','yellow')
end
axis([amps(1) amps(end) 0 1])



subplot(2,3,2);

e=Output.neuronInfo.prefElectrodes{n}(1);
elec = Output.tracesInfo.recElecs(e);
sigma  = Output.sigma(e,:);
try
    plot(abs(elecResp.stimInfo.stimAmps),sigma,...
    'o-','LineWidth',3); grid on; 
catch
     plot(abs(elecResp.stimInfo.stimAmps(1:2:end)),sigma,...
    'o-','LineWidth',3); grid on; 
end

xlabel('stimulation amplitude (\muA)'); 
title(sprintf('Standard deviation of residuals for electrode %0.0f',elec)); 
 
ylabel('recorded daqs'); 
%title(sprintf('activation curve for n%0.0f',neuronId)); 
hold on
for b = breakAxon
    
    plot([amps(b) amps(b)], [min(sigma)*0.8 max(sigma)*1.2],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)], [min(sigma)*0.8 max(sigma)*1.2],'linewidth',1.5,'color','yellow')
end
grid('on')

    axis([amps(1) amps(J) min(sigma)*0.8 max(sigma)*1.2])


subplot(2,3,3)
Artifact=Output.Artifact{e};

Time = [Output.tracesInfo.Trange(1) : Output.tracesInfo.Trange(2)]/Output.params.sampRate*1000;

for j = 1:J
    plot(Time,Artifact(j,:),'color',colors(floor(64/J)*j,:),'linewidth',2);
    hold on
    grid('on')
end
title(sprintf('Artifact estimages for electrode %0.0f',elec)); 
 
xlabel('time (ms)'); 
ylabel('recorded daqs');

% Plot the cell template
subplot(2,3,4); 
xlabel('time (ms)'); 
[thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(Output);
ActivationThreshold = thresholdAlg(n);
index = find(amps >= ActivationThreshold,1,'first');
if isempty(index)
    index = length(amps);
end
ph = plot(Output.ResidualArtifact{index}','r'); 
% times = (1:length(Output.neuronInfo.templates{1}))/20; 
hold on; th = plot(Output.neuronInfo.templates{1},'k','LineWidth',3);
xlim([0 40]); 
legend([ph(1) th],'trials (w artifact subtracted)',...
    sprintf('n%0.0f template', neuronId)); 
title(sprintf(['first condition after human / algorithm disagreement\n'...
    'condition %0.0f , stimAmp %0.3f uA'], index, amps(index)));

subplot(2,3,5);
dots = plot(amps,latencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(amps,latencies(:,3),'x','markersize',15); 
plot(amps,latencies(:,1),'x','markersize',15);
xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (ms)'); 
% legend([dots, crosses],'medians','quartiles'); % of responding trials
title(sprintf('Algorithm  latencies for n%0.0f',neuronId)); 

for b = breakAxon
    
    plot([amps(b) amps(b)], [nanmin(nanmin(latencies))*0.8 nanmax(nanmax(latencies))*1.2],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)], [nanmin(nanmin(latencies))*0.8 nanmax(nanmax(latencies))*1.2],'linewidth',1.5,'color','yellow')
end
grid('on')
if ~all(all(isnan(latencies)))
    axis([amps(1) amps(J) nanmin(nanmin(latencies))*0.8 nanmax(nanmax(latencies))*1.2])
end
% Load elecResp file





% Get latency information from the elecResp files. 
for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
end

 humanLatencies= humanLatencies/Output.params.sampRate*1000;
subplot(2,3,6); 
dots = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,1),'x','markersize',15);
xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (ms)'); 
legend([dots,crosses],'medians','quartiles'); % of responding trials
title(sprintf('Human analyzed latencies for n%0.0f',neuronId)); 

for b = breakAxon
    
    plot([amps(b) amps(b)], [nanmin(nanmin(latencies))*0.8 nanmax(nanmax(latencies))*1.2],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)], [nanmin(nanmin(latencies))*0.8 nanmax(nanmax(latencies))*1.2],'linewidth',1.5,'color','yellow')
end
grid('on')
if ~all(all(isnan(latencies)))
    axis([amps(1) amps(J) nanmin(nanmin(latencies))*0.8 nanmax(nanmax(latencies))*1.2])
end

