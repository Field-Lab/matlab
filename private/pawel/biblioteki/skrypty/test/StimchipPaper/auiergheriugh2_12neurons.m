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
for i=Elektrody
    DataForNeuron=reshape(DataNew(:,i,:),N,Length)+SpikesThresholds(i);
    znak=sign(DataForNeuron);
    Spikes=diff(znak')';
    [traces,samples]=find(Spikes>0);
    
    diffsamples=diff(samples);
    falsesamples=find(diffsamples<40); % jesli spikes wykryte blizej niz 2 ms od siebie - false
    
      samples(falsesamples)=-10;
    
    if i==1
        j=find(samples>8330 & samples<8350);
        samples(j)=-10;
        j=find(samples>13530 & samples<13550);
        samples(j)=-10;
        j=find(samples>10280 & samples<10300);
        samples(j)=-10;
        j=find(samples>10450 & samples<10470);
        samples(j)=-10;
    end
    
    if i==3
        j=find(samples>11860 & samples<11880);
        samples(j)=-10;
        j=find(samples>6810 & samples<6830);
        samples(j)=-10;
        j=find(samples>12410 & samples<12430);
        samples(j)=-10;
        j=find(samples>18430 & samples<18450);
        samples(j)=-10;
    end    
    
    if i==5
        j=find(samples>6210 & samples<6230);
        samples(j)=-10;       
        j=find(samples>18830 & samples<18850);
        samples(j)=-10;       
    end    
    
    if i==6
        j=find(samples>8330 & samples<8350);
        samples(j)=-10;        
    end 
    
    %if i==5
        
    for j=1:length(traces)
        h=plot(samples(j)/20,700-ElectrodeOrder(i)*50+traces(j)*1.8,'bd');
        set(h,'MarkerSize',5);
    end    
end
%plot(DataForNeuron')