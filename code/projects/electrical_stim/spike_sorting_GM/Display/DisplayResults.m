function DisplayResults(input,Gibbs,Log,n)


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
amps=abs(input.stimInfo.listAmps)';

spikeProbs    = nansum(Gibbs.variables.spikes{n}')./I;
spikeLogProbs = Gibbs.variables.Probs(n,:);

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
plot(amps,spikeProbs,'linewidth',3)
grid('on')
hold on
plot(amps,spikeLogProbs,'--','linewidth',2)
for b = breakAxon
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)],[0 1],'linewidth',1.5,'color','red')
end
axis([amps(1) amps(end) 0 1])


subplot(2,2,2)
plot(amps,latencies(:,2),'o','MarkerFaceColor','blue','Markersize',12)
hold on
grid('on')
plot(amps,latencies(:,3),'x','markersize',15)
plot(amps,latencies(:,1),'x','markersize',15)
for b = breakAxon
    
    plot([amps(b) amps(b)], [nanmin(nanmin(lats))-3 nanmax(nanmax(lats))+3],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)], [nanmin(nanmin(lats))-3 nanmax(nanmax(lats))+3],'linewidth',1.5,'color','red')
end
grid('on')
axis([amps(1) amps(J) nanmin(nanmin(lats))-3 nanmax(nanmax(lats))+3])


subplot(2,2,3)

plot(amps,sigma,'linewidth',2)
grid('on')
hold on
for b = breakAxon
    plot([amps(b) amps(b)], [nanmin(sigma) nanmax(sigma)],'linewidth',1.5,'color','green')
end
for b = breakRecElec
    plot([amps(b) amps(b)], [nanmin(sigma) nanmax(sigma)],'linewidth',1.5,'color','red')
end

axis([amps(1) amps(J) nanmin(sigma) nanmax(sigma)])

subplot(2,2,4)

Artifact=Gibbs.variables.ArtifactE{e};

Time = [input.tracesInfo.Trange(1) : input.tracesInfo.Trange(2)]/input.params.sampRate*1000;

for j = 1:J
    plot(Time,Artifact(j,:),'color',colors(floor(1.5*j),:),'linewidth',1.5);
    hold on
    grid('on')
end

