NS_GlobalConstants=NS_GenerateGlobalConstants(500);

switch DataSet
    case 1
        DataPath='J:\analysis\2010-09-14-0\data002_preproc';
        MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie005';
        SPfilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data002\sp_files_2015-05-24';
        ArtifactDataPath=DataPath;
    case 2
        DataPath='G:\analysis\slices\2010-09-14-0\data005_preproc';
        MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie005';
        SPfilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\sp_files-2015-05-24';
        ArtifactDataPath=DataPath;
    otherwise
        error('glaglaglabulbulbul!!');
end

[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,Patterns]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
Movies=1:137;

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

%odczyt:
%fid=fopen(['C:\pawel\nauka\512paper\SpikesAnalysis\sp_p64m1'],'r','ieee-le');
%a=fread(fid,'int32');
%b=reshape(a,length(a)/3,3);
%fclose(fid)