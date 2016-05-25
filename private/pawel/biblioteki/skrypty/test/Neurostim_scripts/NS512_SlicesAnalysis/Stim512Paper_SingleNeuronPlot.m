DataPath='G:\analysis\2010-09-14-0\data002_preproc';
ArtifactDataPath=DataPath;
PatternNumber=339;
MovieNumber=130;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,PatternNumber,MovieNumber,0,0);

figure(101);
plot(reshape(DataTraces(:,467,:),50,600)');

Patterns=[339 371 409 395 463 491 32 64 47 30 45 356];
for p=8%1:length(Patterns)
    Pattern=Patterns(p);
    movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,Pattern);
    for Amplitude=2:17
        Movie=movieIDList(Amplitude);
        [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,Pattern,Movie,0,0);
        subplot(4,4,Amplitude-1);
        plot(reshape(DataTraces(:,467,:),50,600)');
        axis([0 600 -200 100]);
        grid on;
    end
end
        
    