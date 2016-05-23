%fid=fopen('C:\home\Pawel\nauka\Stanford\StreamExamples\stim1el.bin','r','l');
%a=fread(fid,'uint32');
%fclose(fid);

%break
fid=fopen('C:\home\Pawel\nauka\Stanford\StreamExamples\neutral.bin','r','l');
a=fread(fid,'uint32');
fclose(fid);

fid=fopen('C:\home\Pawel\nauka\Stanford\StreamExamples\stim64el.bin','r','l');
b=fread(fid,'uint32');
fclose(fid);


break
fin=a(2:length(a));
fin(150000:160000)=b(150001:160001);
fid=fopen('C:\home\Pawel\nauka\Stanford\StreamExamples\stim1el.bin','wb','l');
b=fwrite(fid,fin,'uint32');
fclose(fid);