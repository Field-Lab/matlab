NumberOfSamples=20000;

%Patterns=[1:64];
%Patterns = LaurenelectrodeOrderGenerator(0)
%Patterns=[Patterns 1:60];
Patterns = [1:108];
TimeShiftInMs=0;
DelayInMs=15;

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
fid = fopen('Movie4','wb')
fwrite(fid,Array,'int32');
fclose(fid);