function y=NS_ShowChannelData(FileName,Channel,TimeStart,TimeLength,Unit,FigureNumber,NS_GlobalConstants);
%A simple function that shows the raw data and the stimulation current
%waveforms on the saeme figure.
%FileName - format like: '003';
Fs=NS_GlobalConstants.SamplingFrequency;
time=[TimeStart+1 TimeStart+TimeLength];

filename_movie=['movie' FileName];
output1=ReadMovieDataChunk4(filename_movie,Channel,time,Unit,NS_GlobalConstants);

t=[TimeStart+1:TimeStart+TimeLength]/Fs*1000;
t=[1:TimeLength]/Fs*1000;

figure(FigureNumber);

full_path=[pwd '\' 'data' FileName]
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
data = rawFile.getData(TimeStart+1,TimeLength);

%subplot(2,1,1);
plotPH(t,output1(1,:).*output1(3,:),100);
clear output1;
y=gca;
xlabel('t');
ylabel('I');
%subplot(2,1,2);
%figure(25);
plot(t,data(:,Channel+1));
clear data;