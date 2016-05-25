NS_GlobalConstants=NS_GenerateGlobalConstants(500);

NeuronFilePath = 'D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002\data002.neurons';
MovieFilePath='C:\home\Pawel\nauka\512stim_paper\MovieFiles\2010-09-14-0\movie002';
OutputPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\SpikeTimesCombined';

DuplicatesFilePath='D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_duplicates.txt';
PrimaryNeurons=NS512_IdentifyPrimaryNeurons_v2(NeuronFilePath,DuplicatesFilePath);

%GaussesFilesPath='C:\home\Pawel\nauka\512stim_paper\NeuronAnalysis\2010-09-14-0\data002\GaussesFiles\';
%FiguresPath='D:\Home\Pawel\analysis\slices\2010-09-14-0\analysis_2015_06_02\figures2\';

MoviesBegins=NS512_MoviesBegins(MovieFilePath,NS_GlobalConstants);

figure(1)
clf
StrangeStim=[];
for n=1:2%length(PrimaryNeurons)
    clf
    %fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2010-09-14-0\data002_Matlab\SpikeTimesCombined\ID=' num2str(PrimaryNeurons(n)) 'c'],'r','ieee-le'); 
    fid=fopen([OutputPath '\ID=' num2str(PrimaryNeurons(n)) 'c'],'r','ieee-le'); 
    a=fread(fid,'int32');
    fclose(fid);
    l=length(a);
    b0=reshape(a,6,l/6);
    pst=b0(4,:);
    b=b0(:,find(pst>20)); %remove very early spikes which are likely artifacts    
        
    p1=b(3,:);
    PatternsStim=unique(p1);
    AllGausses=double([]);
    for i=1:min(length(unique(p1)),64)
        SpikesForPattern=find(b(3,:)==PatternsStim(i));
        PreviousPatterns=b(6,SpikesForPattern);
        p1=find(PreviousPatterns==PreviousPatterns(1));
        ns2(i)=length(SpikesForPattern);
        if ns2(i)>=50
            subplot(8,8,i);           
            HistogramsAll=hist(b(4,SpikesForPattern),[1:600]);
            HistogramsSome=hist(b(4,SpikesForPattern(p1)),[1:600]);
            h1=bar(HistogramsAll);
            hold on
            h2=bar(HistogramsSome);
            set(h2,'edgecolor','g')
            set(h2,'facecolor','g')
            if WersjaFitowania==1
                [CFs,leg,blad]=NS512_FitWithMultiGauss_2015_06_11([1:600],HistogramsAll); 
            else
                [CFs,leg,blad]=NS512_FitWithMultiGauss_2015_06_11b([1:600],HistogramsAll); 
            end
            for g=1:length(CFs)
                AllGausses=[AllGausses' [double(PatternsStim(i)) double(PrimaryNeurons(n)) double(CFs{g}.A) double(CFs{g}.tau) double(CFs{g}.sigma) blad(g)]']';
            end            
            
            for j=1:length(CFs)
                tau(j)=CFs{j}.tau;
                sigma(j)=CFs{j}.sigma;            
                h1=plot(CFs{j}([1:600]),'rd-');
                set(h1,'LineWidth',2); 
                set(h1,'Color','r');
            end
            
            if length(CFs)
                %legend(leg);            
                h=gca;           
                set(h,'XLim',[1 max(tau)*1.3]);
                set(h,'YLim',[0 max(HistogramsAll)*1.3]);
                %yl=get(h,'YLim');
                text(max(tau)*0.5,max(HistogramsAll)*1.37,[num2str(length(p1)) '/' num2str(length(SpikesForPattern)-length(p1))]);
            end
            StrangeStim=[StrangeStim' [PrimaryNeurons(n) PatternsStim(i) PreviousPatterns(1) length(p1) PreviousPatterns(2) length(SpikesForPattern)-length(p1)]']';
            %set(h,'XLim',[1 600]);
        end
    end
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[16 9]);
    set(h,'PaperPosition',[0 0 16 9]); 
    FullFigurePath=[FiguresPath 'n' num2str(PrimaryNeurons(n))];
    %print(h, '-dtiff', '-r120', FullFigurePath);     
    
    fid=fopen([GaussesFilesPath 'AllGausses_n' num2str(PrimaryNeurons(n))],'wb');
    fwrite(fid,AllGausses,'double');
    fclose(fid);        
    
    length(b)/64
    %plot(ns2)
    ile(n)=length(find(ns2>=50));
end

fid=fopen([GaussesFilesPath 'StrangeStim'],'wb');
fwrite(fid,StrangeStim,'double');
fclose(fid);