PulseTimes=unique(MovieData(1:3:length(MovieData)));
DataNew=Data;
for i=1:length(PulseTimes)
    DataNew(:,:,PulseTimes(i):PulseTimes(i)+40)=0;
end

N=20;
%figure(121);
%clf
%hold on;
t=[1:40000];
SpikesThresholds=[55 40 25 80 17.02 60 100]
for i=3%Elektrody
    DataForNeuron=reshape(DataNew(:,i,:),N,Length)+SpikesThresholds(i);
    znak=sign(DataForNeuron);
    Spikes=diff(znak')';
    [traces,samples]=find(Spikes>0);
    
    diffsamples=diff(samples);
    falsesamples=find(diffsamples<40); % jesli spikes wykryte blizej niz 2 ms od siebie - false
    
    %if i==5
        
    for j=1:length(traces)
        h=plot(samples(j)/20,700-ElectrodeOrder(i)*100+traces(j)*3.5,'bd');
        set(h,'MarkerSize',5);
    end    
end
%plot(DataForNeuron')