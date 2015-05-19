function NS_GlobalConstants=NS_GenerateGlobalConstants(electrodeNumber)
% Electrode number must be 512, 61 (64)
switch electrodeNumber
    case 512
        ChipAddresses=31:-1:24;
        NumberOfChannelsPerChip=64;
        ArrayID=500;
    case {61,64}
        ChipAddresses=[30 31];
        NumberOfChannelsPerChip=32;
        ArrayID=1;
    otherwise
        err = MException('MATLAB:specifyElecNo',...
            'Must specify the number of electrodes in the stim system');
        throw(err);
end
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040]; % units are in uA, i.e., 66nA, 266nA, 1.07uA, 4.25uA,....
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
