load C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\NumberOfStimulatedElectrodesMouse
NumberOfStimulatedElectrodesMouse=NumberOfStimulatedElectrodes;

DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
movies=[16:18:448]
for i=1:length(movies)
    amp1(i)=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],13,movies(i),NS_GlobalConstants);
end

subplot('position',[0.78 0.57 0.19 0.38]);
loglog(amp1,sum(NumberOfStimulatedElectrodesMouse,2)/128,'bd-');
grid on