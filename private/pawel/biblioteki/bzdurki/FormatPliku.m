cd C:\home\Pawel\nauka\BTrzpil\FormatyPlikowTesty; %path
fid=fopen('file.bin','w','l');
fwrite(fid,PulseLibrary,'int16'); % PulseLibrary is a Matlab array of I16 type
fclose(fid);

fid1=fopen('file.bin','r','l');
PL=fread(fid1,'int16')';
fclose(fid1);

EventSequence=[1 2 3; 4 5 6; 7 8 9; 10 11 12]; % example for 4 events
cd C:\home\Pawel\nauka\BTrzpil\FormatyPlikowTesty;
fid=fopen('file2.bin','w','l');
fwrite(fid,EventSequence,'int32'); % EventSequence is a Nx3 Matlab array of I32 type (N - number of events)
fclose(fid);

fid1=fopen('file2.bin','r','l');
ESvector=fread(fid1,'int32')'; % here the data is just 1-D vector with N*3 elements
fclose(fid1);
ES=reshape(ESvector,length(ESvector)/3,3); % ES is the recostructed EventSequence array