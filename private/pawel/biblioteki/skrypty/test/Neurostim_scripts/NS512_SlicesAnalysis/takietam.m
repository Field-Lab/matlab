DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';

Electrodes=[236 244 245 412];
Electrodes=[333 334 327 328]
Amplitudes=[1:2:25];
Pattern=61%151;
MovieStart=7;

figure(101)
clf
for i=1:length(Amplitudes)
    Amplitude=Amplitudes(i);
    Movie=MovieStart+(Amplitude-1)*18;
    [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movie,0,0);
    
    for j=1:length(Electrodes)
        data=reshape(DataTraces(:,Electrodes(j),:),50,400);
        subplot(length(Electrodes),length(Amplitudes),length(Amplitudes)*(j-1)+i);
        h1=plot(data')
        set(h1,'Color','b')
        axis([0 100 -500 150])
        grid on
    end
end
    