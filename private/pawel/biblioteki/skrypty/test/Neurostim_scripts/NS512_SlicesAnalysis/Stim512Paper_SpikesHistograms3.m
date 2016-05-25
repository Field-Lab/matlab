NS_GlobalConstants=NS_GenerateGlobalConstants(500);
%AllTraces=zeros(400,600);

%Patterns=[13];

Movies=[8:8:136];

AllGausses=[];

MovieFilePath='C:\pawel\nauka\512paper\movie002';
[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,PatternsInfo]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);
break
for Pattern=AllPatternsUsed
    
    tic
    MoviesForPattern=find(PatternsInfo(:,275)>0));
    AllGausses=[];
    InterestingElectrodes=NS512_SpikesAnalysis_FindElectrodesWithTimeLockedResponses('C:\pawel\nauka\512paper\SpikesAnalysis',13,[8:8:136])
    %tutaj znalexc opodiwednie moviesy!

    for Electrode=InterestingElectrodes
        Electrode
        AllDelays=[];
        for m=1:length(Movies)
            Movie=Movies(m);
            fid=fopen(['C:\pawel\nauka\512paper\SpikesAnalysis\sp_p13m' num2str(Movie)],'r','ieee-le');
            a=fread(fid,'int32');
            b=reshape(a,length(a)/3,3);
            fclose(fid);    
        
            SpikesIDs=find(b(:,1)==Electrode);

            Delays=b(SpikesIDs,3);
            AllDelays=[AllDelays' Delays']';
            
            HistogramsAll=hist(AllDelays,[10:10:600]);
            CFs=NS512_FitWithMultiGauss([1:60],HistogramsAll);
            
            for g=1:length(CFs)
                AllGausses=[AllGausses' [Pattern Electrode CFs{g}.A CFs{g}.tau CFs{g}.sigma]']';
            end
        end
    end
    toc
end 
