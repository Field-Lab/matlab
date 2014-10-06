cd C:\praca\data\;
filename='pattern044';

fid0=fopen(filename,'r','b');
header=readPHchunk(fid0);

%ble=fread(fid0,20,'int8')';
%patterns1=readPDchunk(fid0);

%ble=fread(fid0,16,'int8')';
%patterns2=readPDchunk(fid0);

ble=fread(fid0,16,'int8')';
patterns3=readPDchunk(fid0);
