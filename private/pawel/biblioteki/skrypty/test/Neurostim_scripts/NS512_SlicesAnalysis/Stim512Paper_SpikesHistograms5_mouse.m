NS_GlobalConstants=NS_GenerateGlobalConstants(500);

%SpikesAnalysisPath='C:\pawel\nauka\512paper\SpikesAnalysis';
%SPFilesPath='C:\pawel\nauka\512paper\SpikesAnalysis\';
%GaussesFilesPath='C:\pawel\nauka\512paper\SpikesAnalysis\GaussesFiles2\';
%FiguresPath='C:\pawel\nauka\512paper\SpikesAnalysis\GaussesFigures\';
%MovieFilePath='C:\pawel\nauka\512paper\movie002';

%2013-12-12-0, data001:

GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\GaussesFiles\';
SPFilesPath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\sp_files\';
FiguresPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\GaussesFigures\';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\MovieFiles\2013-12-12-3\movie001';

[NumberOfMovies,NumberOfPatternsPerMovie,AllPatternsUsed,PatternsInfo]=NS512_MoviePatterns(MovieFilePath,NS_GlobalConstants);

NumberOfAmplitudes=25; % 17 is for rat data, 25 for mouse data!
for Pattern=AllPatternsUsed%(1:7)
    tic
    Movies0=find(PatternsInfo(:,Pattern)>0)
    Movies=Movies0([length(Movies0)-NumberOfAmplitudes+1:length(Movies0)])
    AllGausses=[];
    %InterestingElectrodes=NS512_SpikesAnalysis_FindElectrodesWithTimeLockedResponses(SpikesAnalysisPath,Pattern,Movies([1:17]))
    InterestingElectrodes=NS512_SpikesAnalysis_FindElectrodesWithTimeLockedResponses(SPFilesPath,Pattern,Movies,40);
    %tutaj znalexc opodiwednie moviesy!

    for Electrode=InterestingElectrodes%(5)
        Electrode
        AllDelays=[];
        for m=1:NumberOfAmplitudes %%% czasem od 2!!!
            Movie=Movies(m);
            fid=fopen([SPFilesPath 'sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
            a=fread(fid,'int32');
            b=reshape(a,length(a)/3,3);
            fclose(fid);    
        
            SpikesIDs=find(b(:,1)==Electrode);

            Delays=b(SpikesIDs,3);
            AllDelays=[AllDelays' Delays']';           
        end
        HistogramsAll=hist(AllDelays,[1:600]);
        CFs=NS512_FitWithMultiGauss([1:600],HistogramsAll);
        for g=1:length(CFs)
            AllGausses=[AllGausses' [Pattern Electrode CFs{g}.A CFs{g}.tau CFs{g}.sigma]']';
        end
        figure(100)
        clf
        %subplot(2,1,1)
        hold on
        hist(AllDelays,[1:600]);
        for i=1:length(CFs)
            tau(i)=CFs{i}.tau;
            sigma(i)=CFs{i}.sigma;
            h1=plot(CFs{i});
            set(h1,'LineWidth',2);
        end
        grid on
        h=gca
        set(h,'XLim',[0 max(tau)*1.2+5*max(sigma)])
        set(h,'FontSize',14)
        xlabel('Delay [samples]');
        ylabel('hits');
        
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]); 
        FullFigurePath=[FiguresPath 'p' num2str(Pattern) 'e' num2str(Electrode)];
        print(h, '-dtiff', '-r100', FullFigurePath);         
    end
    %break
    fid=fopen([GaussesFilesPath 'AllGausses_p' num2str(Pattern)],'wb');
    fwrite(fid,AllGausses,'double');
    fclose(fid);
    toc
end