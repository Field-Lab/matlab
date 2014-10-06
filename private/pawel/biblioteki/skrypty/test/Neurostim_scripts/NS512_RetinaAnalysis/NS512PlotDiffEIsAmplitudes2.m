Pattern=495;
PrimaryEl=[510];
PrimaryElAplitudes=[40 80];
EIsDataPath=['D:\Home\Pawel\analysis\retina\2012-09-27-4\analysis_2013_08_06\EI_p' num2str(Pattern)];
%EIsDataPath='D:\Home\Pawel\analysis\retina\2012-09-27-4\analysis_2013_08_06\EI_p360';

%pattern 139, neuron 3037

TraceLength=140;

fid2=fopen(EIsDataPath,'r');
a=fread(fid2,'integer*2');
b=reshape(a,32,512,TraceLength);

SamplesOffset=7;

TopHalf=[193:448];
LowHalf=[1:192 449:512];
%ElectrodesToShow=TopHalf;
ElectrodesToShow0=[1:512];

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
ElectrodesToExclude=electrodeMap.getAdjacentsTo(Pattern,Radius)';
ElectrodesToShow=NS_RemoveBadChannels(ElectrodesToShow0,ElectrodesToExclude);

c0=reshape(b(32,:,SamplesOffset+1:TraceLength),512,TraceLength-SamplesOffset);
MinEI=min(c0(ElectrodesToShow,:)')-c0(ElectrodesToShow,133)';
ElectrodesToPlot0=find(MinEI<-50); %indexes within the ElectrodesToShow array
ElectrodesToPlot=ElectrodesToShow(ElectrodesToPlot0); %now we are back to real electrode numbers

%ElectrodesToPlot=[290 294 295 297 298]

AllEIs=zeros(31,512);
AllDiffEIs=zeros(31,512);


%  * * * * wyswietlic wzgledna zmiane amplitudy w funkcji pradu stymulacji * * * * 


for i=23:28
    c1=reshape(b(i,:,SamplesOffset+1:TraceLength),512,TraceLength-SamplesOffset);
    c2=reshape(b(i+1,:,SamplesOffset+1:TraceLength),512,TraceLength-SamplesOffset);    
    EIdiff=c2-c1;
    figure(1);
    MinEI=min(c1')-c1(:,133)';
    AllEIs(i,:)=abs(MinEI);
    MinDiffEI=min(EIdiff');
    AllDiffEIs(i,:)=abs(MinDiffEI);
    plot(abs(MinEI),abs(MinDiffEI),'bd');
    axis([0 400 0 250])
    figure(11)
    clf
    h=plot(abs(MinDiffEI(ElectrodesToPlot)),abs(MinDiffEI(ElectrodesToPlot))./abs(MinEI(ElectrodesToPlot))*100,'bd');
    hold on
    h=plot(abs(MinDiffEI(1,PrimaryEl)),abs(MinDiffEI(1,PrimaryEl))./abs(MinEI(1,PrimaryEl))*100,'rd');
    set(h,'MarkerSize',10)
    set(h,'MarkerFaceColor','r')    
    axis([0 200 0 400])
    grid on
    figure(12)
    clf
    h=semilogy(abs(MinDiffEI(ElectrodesToPlot)),abs(MinDiffEI(ElectrodesToPlot))./abs(MinEI(ElectrodesToPlot))*100,'bd');
    hold on
    h=semilogy(abs(MinDiffEI(1,PrimaryEl)),abs(MinDiffEI(1,PrimaryEl))./abs(MinEI(1,PrimaryEl))*100,'rd');
    set(h,'MarkerSize',10)
    set(h,'MarkerFaceColor','r')    
    axis([0 200 0 2000])
    grid on
    %plot(MinEI,MinDiffEI,'bd');
    pause(1)
    %for j=1:length(ElectrodesToPlot)
    %    HalfWidths(j)=FindNegativeHalfAmplitudeWidth(EIdiff(ElectrodesToShow(ElectrodesToPlot(j)),:));
    %end    
end
break
figure(3)
plot(AllDiffEIs(:,ElectrodesToPlot))
grid on
h=gca
set(h,'XLim',[12 32])

figure(4)
plot(diff(AllDiffEIs(:,ElectrodesToPlot)))

RelativeEIsChanges=AllDiffEIs(2:31,:)./AllDiffEIs(1:30,:);
figure(5)
plot(RelativeEIsChanges(:,ElectrodesToPlot))
h=gca
set(h,'XLim',[12 32])