ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;

%cd C:\praca\data\2008-02-04-0;
%cd C:\praca\data\2008-02-19\2008-02-19-1;
%cd E:\2008-03-18saline;
cd C:\praca\data\2008-03-18saline;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

FileName='000';
%[CI,MI]=NS_Report(FileName,1,NS_GlobalConstants);

StimulationElectrodeNumber=5;
MovieNumber=24;
Fn=2;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);

RecordingElectrodes=[StimulationElectrodeNumber StimulationElectrodeNumber+1];
RecordingElectrodes=[5 7];
RecordingElectrodes = electrodeMap.getAdjacentsTo(StimulationElectrodeNumber,1);
[Pulse,Status,patterns,PatternsIndexes]=NS_FindPulseShapeForMovie(FileName,StimulationElectrodeNumber,MovieNumber,NS_GlobalConstants);
%figure(1);
%clf;
%subplot(2,3,1);
%[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,StimulationElectrodeNumber,0,1,'b-',NS_GlobalConstants);

TimeStart=-10;
NumberOfSamples=70;
BadElectrodes=[9 25 31 33 41 57 64];

%Movies=[34:46];
Movies=[1:33];
figure(2);
for i=[58:64]
    StimulationElectrodeNumber=i;
    RecordingElectrodes = electrodeMap.getAdjacentsTo(StimulationElectrodeNumber,1)
    mask=ones(1,length(RecordingElectrodes));
    for i=1:length(RecordingElectrodes)
        if find(BadElectrodes==RecordingElectrodes(i))
            mask(1,i)=0;
        end
    end
    mask
    RE=RecordingElectrodes(find(mask==1))
    [Ind,Traces]=NS_ReportForOneElectrode(FileName,Movies,StimulationElectrodeNumber,RE,NS_GlobalConstants);
end