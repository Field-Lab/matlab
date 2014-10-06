function Command=NS_RangeCommand(GlobalChipAddress,ChipAddress,GlobalChannelAddress,ChannelAddress,Range);
Header='1010';
ChipAddress=[dec2bin(GlobalChipAddress,1) dec2bin(ChipAddress,5)];
OC='101';
ChannelAddress=[dec2bin(GlobalChannelAddress,1) dec2bin(ChannelAddress,6)];
Range2=dec2bin(Range,3);
[Header ChipAddress OC ChannelAddress Range2 '00000'];
Command=bin2dec([Header ChipAddress OC ChannelAddress Range2 '00000']);