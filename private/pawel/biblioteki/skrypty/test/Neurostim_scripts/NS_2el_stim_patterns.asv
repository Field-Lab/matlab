electrodes=[1:8 10 11];
electrodes=[1:64];
no_of_patterns=64;

cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('electrodes','wb')
fwrite(fid,electrodes,'integer*4');
fclose(fid);

Array=(rand(length(electrodes),no_of_patterns)-0.5)*2;
%Array=eye(no_of_patterns);

fid = fopen('patterns_random','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);