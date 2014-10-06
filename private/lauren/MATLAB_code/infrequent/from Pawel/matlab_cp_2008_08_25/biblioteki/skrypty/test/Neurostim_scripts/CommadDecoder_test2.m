clear;
Fs=20000;

ChipAddresses=[31 30];
NumberOfChannelsPerChip=32;

cd C:/praca;
filename='movie003';
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];

time_start=1540000;
time_length=200000;
time=[time_start+1 time_start+time_length];
channel=6;
unit=1;
%output1=ReadMovieDataChunk2(filename,channel,time,ChipAddresses,NumberOfChannelsPerChip);
output1=ReadMovieDataChunk4(filename,channel,time,unit,ChipAddresses,NumberOfChannelsPerChip,CurrentRanges);

t=[time_start+1:time_start+time_length]/Fs;
%t=[1:26994];
figure(2)
%plotPH(t,output1(1,:));
%plotPH(t,output1(1,:).*output1(3,:),100);
%clear output1;

%gca
%xlabel('t')
%ylabel('I')

'sgsgesrgherh';
javaaddpath C:\praca\vision5-2006-11-27-updated-filter\vision5/vision.jar;
%cd C:\praca;
full_path=[pwd '\' 'data' '003'];
a='C:\praca\data003';
a=[pwd '\' 'data' '003']
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

%rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile('data003');
data = rawFile.getData(time_start+1,time_length);
size(data);

subplot(2,1,1);
plotPH(t,output1(1,:).*output1(3,:),100);
clear output1;
gca;
xlabel('t')
ylabel('I')
subplot(2,1,2);
plot(t,data(:,channel+1));
clear data;