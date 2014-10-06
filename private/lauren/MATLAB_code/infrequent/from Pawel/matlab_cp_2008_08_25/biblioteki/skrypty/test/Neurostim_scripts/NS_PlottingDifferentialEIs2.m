ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

BadElectrodes=[1 9 25 28 31 33 37 41 57 64];
Channel=17;
Fn=19;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
Electrodes = electrodeMap.getAdjacentsTo(Channel,1);
GoodElectrodes=NS_RemoveBadElectrodes(Electrodes,BadElectrodes);

MovieNumber=22;
Diff=4;
nr_el=7;
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
fid=fopen(name,'r');
a=fread(fid,'double');
b=reshape(a,6,nr_el,51);
EI1=-reshape(b(Diff,:,:),nr_el,51);
FigureProperties=struct('FigureNumber',Fn,'Subplot',[2 3 3],'TimeRange',[18 50],'AmplitudeRange',[-100 100],'FontSize',10,'Colors',['r' 'b' 'm' 'k' 'g']);
%WaveformTypes=ones(length(GoodElectrodes));
%y=NS_PlotSignatureOnArrayLayout(EI1,GoodElectrodes,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);
%hj=gcf;
%set(hj, 'PaperOrientation', 'portrait');
%name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\processed\' num2str(Channel) 'Movie' num2str(MovieNumber) 'Difference' num2str(Diff)];
%print(hj, '-dtiff', name);

MovieNumber=27;
Diff=5;
%nr_el=7;
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\PCA_EI\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)]
fid=fopen(name,'r');
a=fread(fid,'double');
b=reshape(a,6,nr_el,51);
EI2=-reshape(b(Diff,:,:),nr_el,51);

MovieNumber=28;
Diff=5;
%nr_el=7;
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\PCA_EI\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)]
fid=fopen(name,'r');
a=fread(fid,'double');
b=reshape(a,6,nr_el,51);
EI3=-reshape(b(Diff,:,:),nr_el,51);

MovieNumber=29;
Diff=5;
%nr_el=7;
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
fid=fopen(name,'r');
a=fread(fid,'double');
b=reshape(a,6,nr_el,51);
EI4=-reshape(b(Diff,:,:),nr_el,51);

MovieNumber=30;
Diff=3;
%nr_el=7;
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
fid=fopen(name,'r');
a=fread(fid,'double');
b=reshape(a,6,nr_el,51);
EI5=reshape(b(Diff,:,:),nr_el,51);

EIs=zeros(2,nr_el,51);
%EIs(1,:,:)=-EI1;
%EIs(1,:,:)=EI2;        
EIs(1,:,:)=EI3;
EIs(2,:,:)=EI4; 
EIs(3,:,:)=-EI5; 

FigureProperties=struct('FigureNumber',Fn+1,'Subplot',[2 3 3],'TimeRange',[18 50],'AmplitudeRange',[-200 200],'FontSize',10,'Colors',['r' 'b' 'm' 'k' 'g']);
y18=NS_PlotManySignaturesOnArrayLayout(EIs,GoodElectrodes,[1:5],ArrayID,FigureProperties,NS_GlobalConstants);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\processed\' num2str(Channel) 'neuron1517_el17_tri_mobies28_30'];
%name=['E:\analysis\2008-03-21-0\PCA_EI\raport\biphasic\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber) 'DifferenceMovies' num2str(MovieNumber) 'and' num2str(MovieNumber2)];
print(hj, '-dtiff', name);