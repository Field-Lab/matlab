load C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\NumberOfStimulatedElectrodesMouse
load C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\NumberOfStimulatedElectrodesRat

figure(1)

subplot(2,2,1)
plot(sum(NumberOfStimulatedElectrodesMouse,2)/128,'bd-');
subplot(2,2,3);
[c,v]=find(diff(sign(NumberOfStimulatedElectrodesMouse))==1);
hist(c,max(c)-min(c)+1);

subplot(2,2,2)
plot(sum(NumberOfStimulatedElectrodesRat,2)/64,'bd-');
subplot(2,2,4);
[c,v]=find(diff(sign(NumberOfStimulatedElectrodesRat))==1);
hist(c,max(c)-min(c)+1);