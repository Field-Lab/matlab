MovieFilePath='D:\Home\Data\slices\2010-09-14-0\movie002';

FullName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\histograms5\GaussParameters.bin']; % pokoj 109
fid=fopen(FullName,'r','ieee-le');
a=fread(fid,'double');
GaussParameters=reshape(a,length(a)/5,5); %NeuronID, interesting electrode, Amplitude, tau, sigma

Patterns=Electrodes(i12);
%Patterns=[30    32    45    47    64   339   339   356   371   395   409   463   491   491];
Patterns=[30    32    45    47    64   339   356   371   395   409   463   491];

figure(11)
clf

LC=[0 0.25 0.5 0.75]+0.02;
BC=[0 0.25 0.5 0.75]+0.02;

LeftCoordinates=[LC(1) LC(2) LC(3) LC(4) LC(4) LC(4) LC(4) LC(3) LC(2) LC(1) LC(1) LC(1)];
BottomCoordinates=[BC(4) BC(4) BC(4) BC(4) BC(3) BC(2) BC(1) BC(1) BC(1) BC(1) BC(2) BC(3)];
SubplotWidth=0.18;
SubplotHeight=0.18;

for p=1%:length(Patterns)
    WholeAreaCoordinates(1)=LeftCoordinates(p);
    WholeAreaCoordinates(2)=BottomCoordinates(p);
    WholeAreaCoordinates(3)=SubplotWidth;
    WholeAreaCoordinates(4)=SubplotHeight;
    
    %subplot('position',[0.05 0.05 0.56 0.92])
    Subsubplot(WholeAreaCoordinates,[0.05 0.05 0.56 0.92]);
    hold on;
    movieIDList = NS512_GetMoviesWithPatternIDNew(MovieFilePath,Patterns(p));
    NumberOfSpikes=zeros(1,length(movieIDList));
    AllSpikeTimes=[];
    for m=1:length(movieIDList)
        SpikesIDs=find(dane(1,:)==movieIDList(m) & dane(3,:)==Patterns(p));
        SpikeTimes=dane(4,SpikesIDs);  
        AllSpikeTimes=[AllSpikeTimes SpikeTimes];
        NumberOfSpikes(m)=length(SpikeTimes);
        for sp=1:length(SpikeTimes)
            h1=plot([SpikeTimes(sp)/20 SpikeTimes(sp)/20],[-m -m-0.5],'b-');
            set(h1,'LineWidth',1);
        end
        h1=gca;
        set(h1,'XLim',[0 30]);
    end
    Subsubplot(WholeAreaCoordinates,[0.7,0.6,0.28,0.38]);
    plot(int32(NumberOfSpikes/51*100),'bd-')
    axis([0 20 0 120]);
    grid on
    
    Subsubplot(WholeAreaCoordinates,[0.7,0.05,0.28,0.48])
    LatenciesForGivenElectrode=AllSpikeTimes%SpikeTimes
    p=hist(LatenciesForGivenElectrode,[5:10:600]);
    hist(LatenciesForGivenElectrode/20,[5:10:600]/20);
    CFs=NS512_FitWithMultiGauss([1:60],p);
    hold on
    FL=CFs{1}
    h=plot([5:10:600]/20,FL([1:60]),'r-')
    set(h,'LineWidth',2)
    h=gca;
    set(h,'XLim',[0 600]/20);
    xlabel('Time [ms]');
    ylabel('N');
    %text(450,40,['A=' num2str(CFs{1}.A,'%8.2f')]);
    text(20,35,['\tau = ' num2str(CFs{1}.tau*10/20,'%8.1f') ' ms']);
    text(20,30,['\sigma = ' num2str(CFs{1}.sigma*10/20,'%8.1f') ' ms']);
    
end
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\figures\fig1_inset.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 5]*2);
set(h,'PaperPosition',[0 0 10 5]*2); 
print(h, '-dtiff', '-r240', FullImageName);