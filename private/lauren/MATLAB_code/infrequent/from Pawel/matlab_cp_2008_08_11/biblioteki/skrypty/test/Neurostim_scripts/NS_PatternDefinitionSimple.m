NumberOfSamples=20000;

Patterns = LaurenelectrodeOrderGenerator(0)
%Patterns=[1:64];
Patterns=[Patterns Patterns Patterns];
Patterns=[1:196];
Patterns=randperm(196);
TimeShiftInMs=0;
DelayInMs=5;

TimeShift=TimeShiftInMs*20; %in sampling periods (50 microseconds)
Delay=DelayInMs*20;
Times=[TimeShift:Delay:TimeShift+Delay*(length(Patterns)-1)];

Array=zeros(1,6+length(Patterns)*3);
Array(1,6)=NumberOfSamples;

for i=1:length(Patterns)
    index=6+(i-1)*3;
    Array(1,index+1)=Times(i); %TimeShift+(i-1)*Delay; %point in time for this pattern
    Array(1,index+2)=Patterns(i);
    Array(1,index+3)=1; %scaling factor - for now always equal to 1
end

Array=[length(Array) Array];

cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('Movie','wb')
fwrite(fid,Array,'int32');
fclose(fid);