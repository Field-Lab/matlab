load C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\NumberOfStimulatedElectrodesMouse
NumberOfStimulatedElectrodesMouse=NumberOfStimulatedElectrodes;
load C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\NumberOfStimulatedElectrodesRat

figure(10)

DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
NS_AmplitudesForPattern_512_1el(DataPath,[1:512],13,16,NS_GlobalConstants);
subplot(1,4,1)
plot(sum(NumberOfStimulatedElectrodesMouse,2)/128,'bd-');
subplot(1,4,2);
[c,v]=find(diff(sign(NumberOfStimulatedElectrodesMouse))==1);
hist(c,max(c)-min(c)+1);

subplot(1,4,3)
plot(sum(NumberOfStimulatedElectrodesRat,2)/64,'bd-');
subplot(1,4,4);
[c,v]=find(diff(sign(NumberOfStimulatedElectrodesRat))==1);
hist(c,max(c)-min(c)+1);

FigurePath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\SpikeAnalysis\2013-12-12-3-PH-2014-08-16';
FullName=[FigurePath '\histogramy.tif']; 
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 4]);o
    set(h,'PaperPosition',[0 0 16 4]); 
    print(h, '-dtiff', '-r120', FullName);