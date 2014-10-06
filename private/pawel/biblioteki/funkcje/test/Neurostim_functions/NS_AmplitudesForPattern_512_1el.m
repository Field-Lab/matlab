function Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,Channels,PatternNumber,MovieNumber,NS_GlobalConstants);
%This specific version of the function includes correction for
%non-clear problem with chip adresses in the log files forthe 512 system.
%This correction works fine for 1-electrode scan data, when all the
%stimulating electrodes in given movie have the same amplitudes.
%The function gives back the amplitude value of stimulation current for
%each electrode.
%DataPath - with preprocessed data

FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
for i=1:512
    r(i)=Status.ChannelsStatus(i).range;
end
for i=1:512
    Status.ChannelsStatus(i).range=max(r);
end
range=Status.ChannelsStatus(PatternNumber).range;
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
A=max(max(abs(Pattern.data)));
Amplitude=A/127*NS_GlobalConstants.CurrentRanges(range+1);