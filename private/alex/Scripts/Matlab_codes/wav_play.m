
[Y,FS,NBITS]=wavread('C:\Documents and Settings\atikidzhi\Desktop\Stratovarius - Phoenix_0.wav');
Y=Y(1:500000,:);
disko=Y;
save('C:\Documents and Settings\atikidzhi\Desktop\MEA_bin_matlab\Matlab\Matlab-Import-Filter\Matlab_Interface\music.mat','disko','FS')
flag=1;
a=dir('\\134.2.181.58\ag_muench\user\alexandra\20110518\MEA_mcd');
while flag
    NoF=length(a);
    pause(120)
    a=dir('\\134.2.181.58\ag_muench\user\alexandra\20110518\MEA_mcd');
    flag=NoF<length(a); 
    display('All OK')
end

while 1
     wavplay(Y,FS)
end