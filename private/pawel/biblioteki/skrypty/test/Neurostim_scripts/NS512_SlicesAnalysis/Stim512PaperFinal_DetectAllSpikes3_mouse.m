NS_GlobalConstants=NS_GenerateGlobalConstants(500);

DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2013-12-12-3\movie001';
SPfilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\sp_files-2015-05-24';
ArtifactDataPath=DataPath;

[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
Movies=2:450;

for Movie=1
    PatternsForMovie=[64   128   129   193   260   324   441   505];
    for p=1:length(PatternsForMovie)
        Pattern=PatternsForMovie(p)
        spikes=[];
        fid=fopen([SPfilesPath '\sp_p' num2str(Pattern) 'm' num2str(Movie)],'wb','ieee-le'); 
        fwrite(fid,spikes,'int32');
        fclose(fid); 
    end
end
break

for Movie=Movies
    PatternsForMovie=find(Patterns(Movie,:)>0);
    for p=1:length(PatternsForMovie)
        p
        Pattern=PatternsForMovie(p)
        [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movie,0,0);
        spikes=NS512_DetectAllSpikes(DataTraces(:,:,:),[1:512],-40.1,-20.1);
        
        fid=fopen([SPfilesPath '\sp_p' num2str(Pattern) 'm' num2str(Movie)],'wb','ieee-le'); 
        fwrite(fid,spikes,'int32');
        fclose(fid);        
    end
end      