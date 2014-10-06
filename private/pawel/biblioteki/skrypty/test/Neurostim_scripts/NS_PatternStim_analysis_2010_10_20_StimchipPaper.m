ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);

Length=40000;

%For neuron 6 (el. 61) must look for spikes smaller than 60!!; for
%electrode 37: make sure we do not capture some other spikes
StimElectrodes=[16 27 28 37 45 51 61];
RecElectrodes=[16 27 28 37 45 50 61];
MovieNumbers=[10 15 21 6 10 ];
Amplitudes=[1 7 4 8 9 3 1]; %[1 6 3 7 8 3 1];
SpikesThresholds=[100 45 40 28 50 30 40]; 
SpikeWidths=[4 3 3 2 4 4 2];
NumberOfElectrodes=length(StimElectrodes);

NumbersOfElectrodes=length(StimElectrodes);

DataPath='F:\2010-09-21-0\data017';
DataMovieFile='F:\2010-09-21-0\movie017';
ArtifactDataPath='F:\2010-09-21-0\data020';
ArtifactMovieFile='F:\2010-09-21-0\movie020';

N=20;
Data=int16(zeros(N,length(RecElectrodes),Length));
ArtifactData=int16(zeros(N,Length));

Times=[];
NumberOfPulses=[];
StimElectrodesIndexes=[];

TimesForHist=zeros(NumbersOfElectrodes,1000);
IloscImpulsow=zeros(1,NumbersOfElectrodes);
Elektrody=1:NumberOfElectrodes;
figure(1)
clf;
subplot('position',[0.15 0.1 0.7 0.8]);
for i=Elektrody %1:NumberOfElectrodes
    StimElectrode=StimElectrodes(i)
    RecElectrode=RecElectrodes(i);
    MovieNumber=2+(Amplitudes(i)-1)*5;
    % 1. Read movie file for RAW data    
    MovieData1=NS512_MovieData2(DataMovieFile,MovieNumber,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData1);
    DataMovieBegin=MovieBegin;
    rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(DataPath);
    
    TimesIndexes=find(MovieData(2:3:length(MovieData))==StimElectrode);
    Times0=MovieData(TimesIndexes*3-2);
    IloscImpulsow(i)=length(Times0);
    Times=[Times' MovieData(TimesIndexes*3-2)']';
    
    % 2. Read movie file for artifact data   
    MovieData=NS512_MovieData2(ArtifactMovieFile,MovieNumber,NS_GlobalConstants);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(MovieData);
    ArtifactMovieBegin=MovieBegin;
    rawArtFile=edu.ucsc.neurobiology.vision.io.RawDataFile(ArtifactDataPath);
    
    for j=1:N
        Start=ArtifactMovieBegin+(j-1)*Length;
        D=int16(rawArtFile.getData(Start,Length)');
        ArtifactData(j,:)=D(RecElectrode+1,:);    
    end
    l=int16(mean(ArtifactData,1));
    
    for j=1:N
        Start=DataMovieBegin+(j-1)*Length;
        D=int16(rawFile.getData(Start,Length)');
        Data(j,i,:)=D(RecElectrode+1,:)-l;                
    end
    
    SpikeIndex=0;
    Dane=reshape(Data(:,i,:),N,Length)+SpikesThresholds(i);
    for j=1:N
        znak=sign(Dane(j,:));
        SpikeTimes=find(diff(znak)>0);
        SpikeTimesNeg=find(diff(znak)<0);
        
        if i==4
            FalseTimes=[690 2610 4130 8250 8410 8530 10610 13130 16370];
            for w=1:length(FalseTimes)
                l=find(SpikeTimes>FalseTimes(w) & SpikeTimes<FalseTimes(w)+5);
                SpikeTimes(l)=-100;
            end                        
        end
        
        if i==7
            znak2=sign(Dane(j,:)-40+55);
            SpikeTimes2=find(diff(znak2)<0);
            
            znak3=sign(Dane(j,:)-40+100);
            SpikeTimes3=find(diff(znak3)<0);
            
            for k=1:length(SpikeTimesNeg)
                ST=SpikeTimesNeg(k);
                st1=find(SpikeTimes2>ST-2 & SpikeTimes2<ST+5);
                st2=find(SpikeTimes3>ST-2 & SpikeTimes3<ST+5);
                
                if (length(st1)>0 & length(st2)==0)
                    in=k;
                    g=find(abs(SpikeTimes-ST)<20); % znajdz w tablicy SpikeTimes (ktora zawiera przekroczenia progu w gore) wartosc najblizsza tej wartosci przekroczenia progu w dol, dla ktorej chcemy usunac spike  zlisty do rysowania
                    SpikeTimes(g)=-100;
                end
            end
        end
        
        for k=1:length(SpikeTimes)-1
            ST=SpikeTimes(k);
            if SpikeTimes(k+1)-SpikeTimes(k)<20
                SpikeTimes(k+1)=-100;
            end
        end
        
        kropki=ones(1,length(SpikeTimes))*(i+j*0.03)*(-1);
        h=plot(SpikeTimes/20,kropki,'bd');
        set(h,'MarkerSize',3);
        hold on;        
        
        for k=1:length(Times0)  
            delays=SpikeTimes-Times0(k);
            l=find(delays>2 & delays<30); % ktore spiki sa opoznione w stosunku do tego konkretnego impulsu o maksymalnie 1.5 ms
            if l
                SpikeIndex=SpikeIndex+1;
                TimesForHist(i,SpikeIndex)=delays(l(1));  
                if delays(l(1))<8
                    [i j k delays(l(1))]                                 
                end
            end      
        end                
    end
    TimesForHist(i,1)=SpikeIndex;
    for j=1:length(Times0)        
        h=plot([Times0(j),Times0(j)]/20,[-i+0.2 -i+0.05]);
        set(h,'Color','r');
        set(h,'LineWidth',2);              
    end
        
end
h=gca;
set(h,'XLim',[0 1000]);
set(h,'FontSize',14);
h=xlabel('Time [ms]');
set(h,'FontSize',14);
h=ylabel('Neuron number');
set(h,'FontSize',14);

% * * * * dodajemy histogramy
for i=1:NumbersOfElectrodes
    subplot('position',[0.88,1-i*0.12,0.10,0.1]);
    l=find(TimesForHist(i,:)>0);
    p=TimesForHist(i,2:TimesForHist(i,1));
    hist(p/20,[1.5:1:25.5]/20);
    h=gca;
    set(h,'XLim',[0.1 1.2]);
    set(h,'YLim',[0 300]);
    m=mean(p);
    text(0.5,250,['delay=' num2str(m/20) 'ms']);
    jit=std(p);
    text(0.8,200,['\sigma=' num2str(round(jit*50),'%d') '\mus']);
    text(0.8,150,['efficacy=' num2str(TimesForHist(i,1)/IloscImpulsow(i)/N*100, '%10.1f')]);
end

FullName=['C:\home\pawel\nauka\analiza\retina\ScatterPlot'];            
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[12 10]);
set(h,'PaperPosition',[0 0 12 10]); 
print(h, '-dtiff', '-r250', FullName);