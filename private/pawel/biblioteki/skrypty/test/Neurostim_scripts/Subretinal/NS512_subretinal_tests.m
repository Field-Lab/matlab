TimeShiftInMs=50;
DelayInMs=50;

TimeShift=TimeShiftInMs*20;
Delay=DelayInMs*20;

NumberOfSamples=10000;
electrodes=[1:512];
Array=eye(512,512); 
%Array=zeros(512,64); 
Times=[];

El=[10:64:512];

Patterns=[];
Times=[];
Disconnecting=zeros(1,8);
TTLs=ones(1,6)*(-1);

addTTL=[0 0 1 1 0 0 1 1]*0;
%addTTL=zeros(1,8);
addDisc=[0 1 0 1 0 1 0 1]*0;
AddElectrodes=[0 0 0 0 1 1 1 1]*1;

%addDisc=zeros(1,8);

T=[50:50:400];

for i=1:8
    if addTTL(i)==1
        p1=TTLs;
        t1=[T(i)+1.4:0.05:T(i)+1.65];
    else
        p1=[];
        t1=[];
    end
    
    if addDisc(i)==1
        p2=Disconnecting;
        t2=[T(i)+1.35:0.05:T(i)+1.7];
    else
        p2=[];
        t2=[];
    end
    
    if AddElectrodes(i)==1
        p3=El(i);
        t3=T(i);
    else
        p3=[];
        t3=[];
    end
    
    Patterns=[Patterns p1 p2 p3];
    Times=[Times t1 t2 t3]; %the patterns for generating the stimulation current will overwrite the disconnstion signals for the stimulating channel, only if they follow the disconnection signal in the stim files
end
Patterns
Times=Times*20

%Patterns=El;
Chunk=NS_MovieChunkGenerationForExperiment(Times,NumberOfSamples,Patterns);
MovieChunksFile=[1 Chunk]; %only one movie

cd C:\home\Pawel\nauka\StimFiles\subretinal\stimfiles; 
fid = fopen('SubTest_el','wb','l')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('SubTest_pt','wb','l')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('SubTest_mv','wb','l')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid); 

break
h2=fopen('ble1','r','l');
s1=fread(h2,'int32');