function result=convert_data(name,header,channels,record);
%Converts the data file "name" to output file "nameconv".
%The input file contains the data in specified format - from
%the "eye" experiment; The output file contains the data
%in the "Matlab-readable" format 'int16' (or 'short'); The
%header is rewritten from the input file;
%The out file has the same size as the input file only
%if the number of samples in every channels is the same;

fid0 = fopen(name,'r');

%2)plik wyjsciowy:
namesize=size(name);
fileoutname(1,1:namesize(2))=name;
fileoutname(1,(namesize(2)+1):(namesize(2)+4))='conv'
fid1=fopen(fileoutname,'r+');

a=fseek(fid0,0,1); %skok na koniec pliku
p=ftell(fid0)

%z=zeros(1,p);
%fwrite(fid1,z,'short');
%3. Naglowek:
a=fseek(fid0,0,-1);
s=fread(fid0,header,'int8');
fwrite(fid1,s,'int8');
clear s;

a=fseek(fid0,0,1); %skok na koniec pliku
p=ftell(fid0)

nrsamples=floor((p-header)/2/channels)
iterations=floor(nrsamples/record)  %number of copy-block-operations

record_2=p-record*iterations*channels*2-header %for the special iteration 
%(the last one, if the number of samples are not equal N*record);

tabin=zeros(2,channels*record); %only for 2-byte sample representation!
tabout=zeros(2*record,channels); %look above

fseek(fid0,header,-1);
for i=1:iterations
   i
   tabin=fread(fid0,[2,channels*record],'int8');
   for k=1:record
      rownr=2*(k-1)+1;
      start=channels*(k-1);
      tabout(rownr,:)=tabin(2,(start+1):(start+channels));
      tabout((rownr+1),:)=tabin(1,(start+1):(start+channels));
   end
   for k=1:channels
      point=header+nrsamples*2*(k-1)+record*2*(i-1);
      fseek(fid1,point,-1);
      fwrite(fid1,tabout(:,k),'int8');
   end
   p=ftell(fid1)
end
clear tabin;
clear tabout;

if (record_2~=0)
   tabin=zeros(2,channels*record_2);
   tabout=zeros(2*record_2,channels);
   
   tabin=fread(fid0,[2,channels*record_2],'int8');
   for k=1:record_2
      rownr=2*(k-1)+1;
      start=channels*(k-1);
      tabout(rownr,:)=tabin(2,(start+1):(start+channels));
      tabout((rownr+1),:)=tabin(1,(start+1):(start+channels));
   end
   for k=1:channels
      point=header+nrsamples*2*(k-1)+record*2*iterations;
      fseek(fid1,point,-1);
      fwrite(fid1,tabout(:,k),'int8');
   end
   p=ftell(fid1)
   clear tabin;
   clear tabout;
end

fclose(fid0);
fclose(fid1);

result='ok';
      
   