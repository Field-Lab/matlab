%full field
mu=30;
sigma=9;
stim_time_min=1/3;
max_frames=round(stim_time_min*3600);
a=sigma*randn(max_frames,1)+mu;
a(a<0)=0;
a(a>60)=60;
b=a;
a=round(a);
a=int8(a);
hist(double(a),59)
mean(a)
std(double(a))


repeated=a;

stim_time_min=1/6;
max_frames=round(stim_time_min*3600);
for i=1:50    
    a=sigma*randn(max_frames,1)+mu;
    a(a<0)=0;
    a(a>60)=60;
    seqB(1:1800,i)=[b; a];
    a=round(a);
    a=int8(a);
%     mean(a);
%     std(double(a));
    seq(1:1800,i)=[repeated; a];
end


sigmaLow=1.8;
seqLow=(seqB-mu)/sigma*sigmaLow+mu;
seqLow=round(seqLow);
seqLow(seqLow<0)=0;
seqLow(seqLow>60)=60;
seqLow=int8(seqLow);
std(double(seqLow(:)))
mean(seqLow(:))



cd('S:\user\alexandra\stimfiles\_Stefano_\');
for i=1:50
    fid=fopen(['S:\user\alexandra\stimfiles\_Stefano_\LCStefano_',int2str(i),'.stim.txt'],'w');
    fprintf(fid,'Parameters:\r\nID\tduration\tposx\tpoy\r\n-1\t1000\t0\t0\r\n');
    for j=1:size(seqLow,1)
        fprintf(fid,['-1\t1\t2\t',num2str(seqLow(j,i)),'\r\n']);
        
    end
    fprintf(fid,'-1\t1000\t0\t0');
    fclose(fid);
end



% create batch file

cd('S:\user\alexandra\stimfiles\_Stefano_\');
nm='NatMov1.stim.txt';

fid_batch=fopen('flickerNM1.batch.txt','w');
fprintf(fid_batch,'nd\tnone\t0\r\n');

for k=1:50

    fprintf(fid_batch,['HCStefano_',int2str(k),'\tnone\t0\r\n']);
    fprintf(fid_batch,['LCStefano_',int2str(k),'\tnone\t0\r\n']);
    fprintf(fid_batch,'NatMov1\tnone\t0\r\n');
end
for k=1:5
    fprintf(fid_batch,'chirp\tnone\t0\r\n');
end
fprintf(fid_batch,'##BeginStimParam##\r\nCounter\r\n')
for i=1:50
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'\r\n');
    fprintf(fid_batch,'1\r\n');
end
for k=1:5
    fprintf(fid_batch,'30\t30\t500\r\n');
end
fprintf(fid_batch,'##BeginVoltParam##')
for i=1:156
    fprintf(fid_batch,'\r\n');
end
fprintf(fid_batch,'##BeginCounterSettings##\r\n');
fclose(fid_batch);

