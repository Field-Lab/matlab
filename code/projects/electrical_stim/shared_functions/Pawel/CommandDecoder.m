function NewStatus=CommandDecoder(ChipAddresses,NumberOfChannelsPerChip,status,command)
%
% updates values in NewStatus.ChannelsStatus.range (all other fields taken
% directly from 'status')
%

%command - I32 value
%channels are numbered from 1 here!!!
NewStatus=status;

%1) I32 value to bit stream:
bit=zeros(1,32); %value of 32 comes from the fact that each command is encoded in the file in single I32 value
for i=1:32
    bit(i)=floor(command/(2^(32-i)));
    command=command-bit(i)*2^(32-i);
end
clear command;

%2) Looking for header - for now works onlu for 4-bit margin at the
%begining! (2007-12-01)
margin=4; %this value comes from the Labview software!!! (STimbuild_stream.vi)
header=bit(margin+1:margin+4);

%3) Looking for chips with physical address identical to the address in the
%command
ChipAddress=bit(margin+5:margin+10);
ChipAddressGlobal=ChipAddress(1);
ChipAddressLocal=ChipAddress(2)*16+ChipAddress(3)*8+ChipAddress(4)*4+ChipAddress(5)*2+ChipAddress(6);
ChipNumbers=find(ChipAddresses==ChipAddressLocal);

if ChipAddressGlobal==1
    ChipNumbers=[1:1:length(ChipAddresses)]
end
if length(ChipNumbers)==0
    error ('Wrong definition of ChipAddresses - nono of the chips has the correct address');
end

%4) Realizing the command
%4a) Operation Code:
OC=bit(margin+11:margin+13);
code=OC(1)*4+OC(2)*2+OC(3);

for i=ChipNumbers       
    switch code %depending on the Operation Code... NOW WORKING WELL ONLY FOR SOME COMMANDS
        case 0      
            disp('StimActive');
        case 1
            disp('Test');
        case 2
            disp('StimStart');
        case 3
            disp('Stop');
        case 4
            disp('SoftReset');
        case 5
            %disp('StimRange');
            ChannelAddress=bit(margin+14:margin+20); %channel address in command
            
            ChannelAddressGlobal=ChannelAddress(1); %if global (for this chip)
            ChannelAddressLocal=ChannelAddress(2)*32+ChannelAddress(3)*16+ChannelAddress(4)*8+ChannelAddress(5)*4+ChannelAddress(6)*2+ChannelAddress(7)+1; %local for this chip - numbered from 1!!!!
            if ChannelAddressGlobal==1
                ChannelAddressLocal=[0:NumberOfChannelsPerChip-1];
            end                            
            ChannelNumbers=NumberOfChannelsPerChip*(i-1)+ChannelAddressLocal; %addresses of all the channels 
            %(in given chip) that will be affected - in fact, there may be
            %one or 64 such channels
                      
            Range=bit(margin+21)*4+bit(margin+22)*2+bit(margin+23); %reading in the range value in the command        
            for j=ChannelNumbers                
                NewStatus.ChannelsStatus(j).range=Range;
            end                            
        case 6
            disp('StimDelay');
        case 7
            disp('SwitchDelay');
    end
end