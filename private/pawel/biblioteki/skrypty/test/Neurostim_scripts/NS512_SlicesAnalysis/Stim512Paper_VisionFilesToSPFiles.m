NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie002';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
OutputPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\SpikeTimesCombined';

MoviesBegins=NS512_MoviesBegins(MovieFilePath,NS_GlobalConstants);

SpikeTimesCombined=NS512_SpikeTimesToStimulationParameters_v2b(NeuronFilePath,MovieFilePath,DuplicatesFilePath,OutputPath);

break
NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);
FiguresPath='D:\Home\Pawel\analysis\slices\2010-09-14-0\analysis_2015_06_02\figures\';

figure(1)
clf
%ciekawe przypadki: 232, 
for n=132%1:length(PrimaryNeurons)
    clf
    fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_Matlab\SpikeTimesCombined\ID=' num2str(PrimaryNeurons(n)) 'c'],'r','ieee-le'); 
    a=fread(fid,'int32');
    fclose(fid);
    l=length(a);
    b0=reshape(a,5,l/5);
    pst=b0(4,:);
    b=b0(:,find(pst>20)); %remove very early spikes which are likely artifacts    
        
    p1=b(3,:);
    PatternsStim=unique(p1);
    for i=2:min(length(unique(p1)),25)
        ns2(i)=length(find(b(3,:)==PatternsStim(i)));
        if ns2(i)>=50
            subplot(5,5,i)
            hist(b(4,find(b(3,:)==PatternsStim(i))),[1:600]);
            hold on
            HistogramsAll=hist(b(4,find(b(3,:)==PatternsStim(i))),[1:600]);
            [CFs,leg,blad]=NS512_FitWithMultiGauss_2015_05_26([1:600],HistogramsAll); 
            
            for j=1:length(CFs)
            tau(j)=CFs{j}.tau;
            sigma(j)=CFs{j}.sigma;            
            %h1=plot(CFs{j}([1:600]),'rd-');
            %set(h1,'LineWidth',2); 
            %set(h1,'Color','r'); 
            end
            
            h=gca;
            set(h,'XLim',[1 max(tau)*1.3]);
            %set(h,'XLim',[1 600]);
        end
    end
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    FullFigurePath=[FiguresPath 'n' num2str(PrimaryNeurons(n))];
    print(h, '-dtiff', '-r120', FullFigurePath);     
    
    length(b)/64
    %plot(ns2)
    ile(n)=length(find(ns2>=50));
end