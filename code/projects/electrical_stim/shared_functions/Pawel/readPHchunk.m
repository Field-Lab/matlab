function header=readPHchunk(fid);
%Reads the header-of-pattern-file chunk and returns the two-field structure:
%header=struct(number_of_chunks,time_unit), where 'number_of_chunks' is the
%total number of Pattern Data (PD) chunks in the file and the 'time_unit' value
%is the time unit (in nanoseconds) for each pattern in the file as well as
%for each movie in the corresponding movie file. The time unit is by
%default equal to the sampling frequency of the DACs in Stimchips, which is
%typically equal to 50 microseconds (then the 'time_unit' value would be
%50000).

%fid0=fopen(filename,'r','b');
ID=fread(fid,8,'int8')';

if (ID~=[-15 110 -9 -21 88 6 -76 43])
   error('chunk ID not correct');
end

header=ID;
size=fread(fid,1,'int64');

number_of_chunks=fread(fid,1,'int32');
time_unit=fread(fid,1,'int32');

header=struct('number_of_chunks',number_of_chunks,'time_unit',time_unit);