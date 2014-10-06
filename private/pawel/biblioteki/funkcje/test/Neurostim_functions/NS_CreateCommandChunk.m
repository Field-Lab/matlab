function Commands=NS_CreateCommandChunk(Status);
%Input status includes the status of each channel!!
NS_GlobalConstants=NS_GenerateGlobalConstants(61);

ChannelsStatus=Status.ChannelsStatus;

NumberOfCommands=length(ChannelsStatus);
for i=1:NumberOfCommands
    range=ChannelsStatus(i).range;
    ChipNumber=floor((i-1)/NS_GlobalConstants.NumberOfChannelsPerChip)
    ChipAddress=NS_GlobalConstants.ChipAddresses(ChipNumber+1);
    Channel=i-1-ChipNumber*NS_GlobalConstants.NumberOfChannelsPerChip
    Commands(i)=NS_RangeCommand(0,ChipAddress,0,Channel,range);
end