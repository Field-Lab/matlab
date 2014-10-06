ChipAddresses=[31 30];
NumberOfChannelsPerChip=32;

NumberOfChips=length(ChipAddresses);
NumberOfChannels=NumberOfChips*NumberOfChannelsPerChip;

ChipStatus=struct('StimActive',0','ConnectDelay',100,'RefreshDelay',100,'RecordDelay',110,'DischargeDelay',105,'HoldDelay',115);
ChipsStatus(1:NumberOfChips)=ChipStatus;

ChannelStatus=struct('active',0,'mode',0,'range',3);
ChannelsStatus(1:NumberOfChannels)=ChannelStatus;

Status=struct('ChipsStatus',ChipsStatus,'ChannelsStatus',ChannelsStatus);

command=176065920;
NewStatus=CommandDecoder(ChipAddresses,NumberOfChannelsPerChip,Status,command)

NewStatus.ChannelsStatus(13)