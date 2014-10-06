function PrimaryStimulatedElectrode=NS_FindPrimaryStimulatedElectrode(mins,Threshold);

[Events,Electrodes]=find(mins<-Threshold);
Counts=histc(Electrodes,[1:512]); %how many times the threshold was crossed on every electrode
BestElectrodes=find(Counts=max(Counts)); %usually this gives just one electrode that shows more threshold crossings than other electrodes. However, sometimes several electrodes can give the same number of threshold crossings.
PrimaryStimulatedElectrode=BestElectrodes(1);

if length(BestElectrodes)>1
    
    for i=2:length(BestElectrodes)
        Electrode=BestElectrodes(i);
        EventsForElectrode=find(Electrodes==Electrode);
        Amplitudes=sum(mins(Events(EventsForElectrode),EventsForElectrode);
        
        
    