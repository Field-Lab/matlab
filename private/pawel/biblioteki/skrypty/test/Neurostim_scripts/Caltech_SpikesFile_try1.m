% 1. Generation of artificial spike times
spikes=[1 1 10 25 40
    1 1 10 28 35
    1 1 20 28 45
    2 1 5 55 35
    2 1 7 55 38
    2 1 11 55 31
    2 1 15 55 36
    2 1 15 57 42
    15 1 10 55 30
    15 1 15 28 40]

% 2. Writing the file
fid1=fopen('C:\pawel\nauka\Caltech\spikes_files\SpikesFile1','w+');
fwrite(fid1,spikes,'integer*2');
fclose(fid1);

% 3. Loading the data from the file
fid2=fopen('C:\pawel\nauka\Caltech\spikes_files\SpikesFile1','r');
a=fread(fid2,'integer*2'); %dane 1-8: pierwszy przypadek 
fclose(fid2);
b=reshape(a,length(a)/5,5);