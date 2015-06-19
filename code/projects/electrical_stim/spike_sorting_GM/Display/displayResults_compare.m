function displayResults_compare(input,Gibbs,Log,n,varargin)
% Displays the result of Gonzalo's spike sorting and compares it to the
% result from human analysis. There must exist an elecResp file in order to
% run this function
% if the optional argument is one, then the x axis show amplitude condition
% instead of stimulus amplitude



figure; set(gcf,'Position', [387	340     1411    758]); 

neuronIds = input.neuronInfo.neuronIds;
neuronId = neuronIds(n);
e = input.neuronInfo.prefElectrodes{n}(1);
recElecs     = input.tracesInfo.recElecs;
breakAxon    =  input.tracesInfo.breakAxon{e};
breakRecElec = input.tracesInfo.breakRecElecs{e};
colors       = colormap(jet);
breakPointsAll = sort(unique([breakAxon breakRecElec]));
J = input.tracesInfo.J;
I = input.tracesInfo.I;
T = input.tracesInfo.T;

spikeProbs    = nansum(Gibbs.variables.spikes{n}')./I;
spikeLogProbs = Gibbs.variables.Probs(n,:);
% Define optional variable
if(nargin == 5)
    if(varargin{1} == 0)
        amps = abs(input.stimInfo.listAmps(:,1))';
    else
        amps = [1:J];
    end
else
     amps=abs(input.stimInfo.listAmps(:,1))';
end
clear latencies
for j=1:J
    lats = Gibbs.variables.latencies{n}(j,1:I(j));
    lats = lats(lats>0);
    if(isempty(lats))
        latencies(j,1:3)=NaN;
    else
        
        latencies(j,1) = prctile(lats,25);
        latencies(j,2) = prctile(lats,50);
        latencies(j,3) = prctile(lats,75);
    end
end
clear sigma
sigma  = Gibbs.variables.sigma(e,:);

subplot(2,2,1)
plot(amps,spikeProbs,'.-','linewidth',3); 
grid('on'); 
hold on; plot(amps,spikeLogProbs,'--','linewidth',2); 
xlabel('stimulation amplitude (\muA)'); 
ylabel('spike probability'); 
title(sprintf('activation curve for n%0.0f',neuronId)); 

for b = breakAxon
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','red')
end
axis([amps(1) amps(end) 0 1])


subplot(2,2,2);
dots = plot(amps,latencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on
grid('on')
crosses = plot(amps,latencies(:,3),'x','markersize',15); 
plot(amps,latencies(:,1),'x','markersize',15);
xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (sample points)'); 
% legend([dots, crosses],'medians','quartiles'); % of responding trials
title(sprintf('latencies for n%0.0f',neuronId)); 

for b = breakAxon
    
    plot([amps(b) amps(b)], [nanmin(nanmin(latencies))-3 nanmax(nanmax(latencies))+3],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)], [nanmin(nanmin(latencies))-3 nanmax(nanmax(latencies))+3],'linewidth',1.5,'color','red')
end
grid('on')
if ~all(all(isnan(latencies)))
    axis([amps(1) amps(J) nanmin(nanmin(latencies))-3 nanmax(nanmax(latencies))+3])
end
% Load elecResp file
pathname = input.names.path; 
fname = ['elecResp_n' num2str(neuronId) '_p' num2str(input.stimInfo.patternNo) '.mat']; 
filename = fullfile(pathname,fname); 
temp = load(filename); 
elecResp = temp.elecResp; 

subplot(2,2,3); 
hum = plot(abs(elecResp.stimInfo.stimAmps),elecResp.analysis.successRates,...
    'o-','LineWidth',3); grid on; 
hold on; alg = plot(amps,spikeProbs,'^-','linewidth',3); 
legend([hum alg(1,1)],'human','algorithm'); 
xlabel('stimulation amplitude (\muA)'); 
title('human analysis results'); 
ylim([0 1]); 

% Get latency information from the elecResp files. 
for ii = 1:size(elecResp.analysis.latencies,1)
    humanLats = elecResp.analysis.latencies{ii}; 
    humanLats = humanLats(humanLats>0); 
    humanLatencies(ii,1) = prctile(humanLats,25);
    humanLatencies(ii,2) = prctile(humanLats,50);
    humanLatencies(ii,3) = prctile(humanLats,75);
end
subplot(2,2,4); 
dots = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,2),'o','MarkerFaceColor','blue','Markersize',12); 
hold on; grid on; 
crosses = plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,3),'x','markersize',15); 
plot(abs(elecResp.stimInfo.stimAmps),humanLatencies(:,1),'x','markersize',15);
xlabel('stimulation amplitude (\muA)'); 
ylabel('latencies (sample points)'); 
legend([dots,crosses],'medians','quartiles'); % of responding trials
title(sprintf('human analyzed latencies for n%0.0f',neuronId)); 

