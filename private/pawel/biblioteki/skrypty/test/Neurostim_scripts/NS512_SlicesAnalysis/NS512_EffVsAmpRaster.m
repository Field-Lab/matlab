MovieFilePath='D:\Home\Data\slices\2010-09-14-0\movie002';

Patterns=Electrodes(i12);
Patterns=[30    32    45    47    64   339   339   356   371   395   409   463   491   491];

figure(11)
clf
%subplot('position',[0.05 0.05 0.6 0.92])
%hold on;
LC=[0 0.25 0.5 0.75]
BC=[0 0.25 0.5 0.75];

LeftCoordinates=[LC(1) LC(2) LC(3) LC(4) LC(4) LC(4) LC(4) LC(3) LC(2) LC(1) LC(1) LC(1)]
BottomCoordinates=[BC(4) BC(4) BC(4) BC(4) BC(3) BC(2) BC(1) BC(1) BC(1) BC(1) BC(2) BC(3)]
SubplotWidth=0.2;
SubplotHeigth=0.2;

for p=10%:length(Patterns)
    
    subplot('position',[0.05 0.05 0.56 0.92])
    hold on;
    movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,Patterns(p));
    NumberOfSpikes=zeros(1,length(movieIDList));
    for m=1:length(movieIDList)
        SpikesIDs=find(dane(1,:)==movieIDList(m) & dane(3,:)==Patterns(p));
        SpikeTimes=dane(4,SpikesIDs);  
        NumberOfSpikes(m)=length(SpikeTimes);
        for sp=1:length(SpikeTimes)
            h1=plot([SpikeTimes(sp)/20 SpikeTimes(sp)/20],[-m -m-0.5],'b-');
            set(h1,'LineWidth',1);
        end
        h1=gca;
        set(h1,'XLim',[0 30]);
    end
    subplot('position',[0.7,0.6,0.28,0.38])
    plot(int32(NumberOfSpikes/51*100),'bd-')
    axis([0 20 0 120]);
    grid on
    
    subplot('position',[0.7,0.05,0.28,0.48])
    LatenciesForGivenElectrode
    p=hist(LatenciesForGivenElectrode,[5:10:600]);
    hist(LatenciesForGivenElectrode,[5:10:600]);
    CFs=NS512_FitWithMultiGauss([1:60],p);
    hold on
    FL=CFs{1}
    h=plot([5:10:600],FL([1:60]),'r-')
    set(h,'LineWidth',2)
    h=gca;
    set(h,'XLim',[0 600]);
    xlabel('Time [ms]');
    ylabel('N');
    %text(450,40,['A=' num2str(CFs{1}.A,'%8.2f')]);
    text(400,35,['\tau = ' num2str(CFs{1}.tau*10/20,'%8.1f') ' ms']);
    text(400,30,['\sigma = ' num2str(CFs{1}.sigma*10/20,'%8.1f') ' ms']);
    
end
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\figures\fig1_inset.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[8 4.5]);
set(h,'PaperPosition',[0 0 8 4.5]); 
print(h, '-dtiff', '-r120', FullImageName);