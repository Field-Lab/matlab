NS_GlobalConstants=NS_GenerateGlobalConstants(500);

NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie002';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
OutputPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\SpikeTimesCombined';

MoviesBegins=NS512_MoviesBegins(MovieFilePath,NS_GlobalConstants);

GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFiles\';
%NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);
FiguresPath='D:\Home\Pawel\analysis\slices\2010-09-14-0\analysis_2015_06_02\figures_alt\';

figure(1)
clf
%ciekawe przypadki: 232, 
StrangeStim=[];
for n=180%length(PrimaryNeurons)
    clf
    fid=fopen([OutputPath '\ID=' num2str(PrimaryNeurons(n)) 'c'],'r','ieee-le'); 
    a=fread(fid,'int32');
    fclose(fid);
    l=length(a);
    b0=reshape(a,6,l/6);
    pst=b0(4,:);
    b=b0(:,find(pst>20)); %remove very early spikes which are likely artifacts    
    
    p1=b(3,:); % ostatnia stymuluj?ca elektroda
    PatternsStim=unique(p1);               
    AllGausses=double([]);
    for i=1:min(length(unique(p1)),25)
        %subplot(5,5,i);
        SpikesForPattern=find(b(3,:)==PatternsStim(i));    
        ns2(i)=length(SpikesForPattern);
        if ns2(i)>50
            PreviousPatterns=unique(b(6,SpikesForPattern));
            HistogramsForPattern=hist(b(4,SpikesForPattern),[1:600]);
            plot(HistogramsForPattern,'bd-');
            hold on;
        
            for j=1:length(PreviousPatterns)
                SpikesforPreviousPattern=find(b(6,:)==PreviousPatterns(j));
                HistogramsForPreviousPattern=hist(b(4,SpikesforPreviousPattern),[1:600]);
                h2=plot(HistogramsForPreviousPattern,'rd-');                
            end
        end
    end
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    FullFigurePath=[FiguresPath 'n' num2str(PrimaryNeurons(n)) '_v2'];
    print(h, '-dtiff', '-r120', FullFigurePath);  
end
