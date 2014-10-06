%clear;
Fs=20000;
%javaaddpath C:\praca\vision5-2006-11-27-updated-filter\vision5/vision.jar;
%javaaddpath C:\praca\vision5-2006-11-27-updated-filter\vision5\vision.jar;

ChipAddresses=[31 30];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

cd E:\2008-02-19-1;
filename='013';
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
unit=1;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

%[CI,MI]=NS_Report(filename,unit,NS_GlobalConstants);

TimeStart=14773000;
TimeLength=10000;
%time=[time_start+1 time_start+time_length];
Channel=13;
Unit=1;

FigureNumber=30;
%y=NS_ShowChannelDataExtended(CI,MI,filename,Channel,[],TimeStart,TimeLength,Unit,FigureNumber,NS_GlobalConstants);

y=NS_ShowChannelData(filename,Channel,TimeStart,TimeLength,Unit,FigureNumber,NS_GlobalConstants)
h=gca;
%set(h,'FontSize',13);
%xlabel('t [ms]');
%ylabel('ADC units');
%axis([3257 3260 -800 0]);
%grid on;
