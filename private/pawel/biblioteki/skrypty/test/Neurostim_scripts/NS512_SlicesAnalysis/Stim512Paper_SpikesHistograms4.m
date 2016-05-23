NS_GlobalConstants=NS_GenerateGlobalConstants(500);
GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\rat005_GaussesFiles600points\';

AllGausses2=[];

MovieFilePath='G:\backups_2013_12\cortex_backup_3\2010-09-14-0\movie005'; %'C:\pawel\nauka\512paper\movie002';
[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,PatternsInfo]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);

for Pattern=AllPatternsUsed(3:length(AllPatternsUsed))
    tic
    Movies=find(PatternsInfo(:,Pattern)>0)
    AllGausses=[];
    InterestingElectrodes=NS512_SpikesAnalysis_FindElectrodesWithTimeLockedResponses('C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\rat005_sp_files',Pattern,Movies([1:17]))
    %tutaj znalexc opodiwednie moviesy!

    for Electrode=InterestingElectrodes        
        Electrode
        AllDelays=[];
        for m=1:17
            Movie=Movies(m);
            fid=fopen(['C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\rat005_sp_files\sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
            a=fread(fid,'int32');
            b=reshape(a,length(a)/3,3);
            fclose(fid);    
        
            SpikesIDs=find(b(:,1)==Electrode);

            Delays=b(SpikesIDs,3);
            AllDelays=[AllDelays' Delays']';
                        
        end
        HistogramsAll=hist(AllDelays,[11:1:600]);
        CFs=NS512_FitWithMultiGauss([11:1:600],HistogramsAll);
        for g=1:length(CFs)
            AllGausses=[AllGausses' [Pattern Electrode CFs{g}.A CFs{g}.tau CFs{g}.sigma]']';
        end
    end
    %AllGausses2=[AllGausses2 AllGausses];
    fid=fopen([GaussesFilesPath 'AllGausses_p' num2str(Pattern)],'wb');
    fwrite(fid,AllGausses,'double');
    fclose(fid);
    toc
end 