%test parameters
currVec = [3 .1 4];
iter = 27;
ratios = [2 -3 1];
time = [50 50 50];
randflg = 1;
gapParam = 7500; %in microseconds
saveFile = 'test';
saveLocation = '/Users/sasidharmadugula/Desktop';

processGuiOutput(currVec,iter,ratios,time, randflg, gapParam, saveFile, saveLocation)

%load files again
fid = fopen([saveLocation filesep saveFile '.slf'],'r','l'); 
grot1 = fread(fid,'int16');
fclose(fid);

fid = fopen([saveLocation filesep saveFile '.sif'],'r','l'); 
grot2 = fread(fid,'int32');
fclose(fid);

fid = fopen([saveLocation filesep saveFile '.sef'],'r','l'); 
grot3 = fread(fid,'int32');
fclose(fid);
grot3 = reshape(grot3,length(grot3)/3,3);
