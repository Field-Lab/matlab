function patterns_out=readPDchunk(fid);
%Loads the whole PD chunk data and returns them in 1-dimensional array.
%Each element of the array is a structure with two fields:
%1. channel - channel number
%2. data - which is a structure with elements identical to total length of the
%pulse in given channel, where the unit is one sampling period. Each
%element has five fields:
%2a) DAC - DAC value
%2b) record - state of "record" switch, may be 0 or 1
%2b) connect - state of "connect" switch, may be 0 or 1
%2b) discharge - state of "discharge" switch, may be 0 or 1
%2b) hold - state of "hold" switch, may be 0 or 1;

ID=fread(fid,8,'int8')'

%if (ID~=[106 -23 100 -113 -38 79 0 -93])
%   error('chunk ID not correct');
%end

size=fread(fid,1,'int64')

patterns=fread(fid,size,'int32') %loads all the data from the chunk
i=1;

patterns_out=[];

while i<size % for every patter...
    pattern=[]
    number_of_channels=patterns(i);
    i=i+1;
    for j=1:number_of_channels % for every channel in the pattern...
        channel=patterns(i)
        i=i+1;
        number_of_phases=patterns(i)
        i=i+1;
        scaling_factor=patterns(i)
        i=i+1;
        index=1;
        
        durations=zeros(1,number_of_phases);      
        DAC_and_switches=zeros(1,number_of_phases);
        
        index0=1;
        for k=1:number_of_phases % for each phase of the pulse in this channel...
            DAC_and_switches=patterns(i)     ;       
            i=i+1;            
            duration=patterns(i);
            i=i+1;
            index1=index0+duration;            
         
            pattern_channel(index0:index1-1)=DAC_and_switches_decode(DAC_and_switches);
            channel_data=struct('channel',channel,'data',pattern_channel); 
          
            index0=index1;                                   
        end
        pattern=[pattern channel_data];
    end
    patterns_out=[patterns_out pattern];
end                                                                