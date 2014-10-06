function r=readconv(name,header,nrchns,channels,samples);
%r=readconv(name,header,nrchns,channels,samples);
%
%Read the values from the "name" file, created using the "convert_data"
%function. Other inputs:
%header (in bytes) - usually 206;
%nrchns - number of channels in file - usually 65;
%channels - vector contains numbers of channels to read the data;
%samples - vector contains the numbers of first and last sample to read
%from each channel specified in "channel" argument;
%The output table contains the blocks of samples (one row is one
%channel).

as=size(samples);
if (as(1)>1)
   samples=samples';
   as=as';
end

if (as(1)>1)
   error('The input argument "samples" has contains exactly two values');
end

if (samples(1,1)>=samples(1,2))
   error('The second value in "samples" has to be greater than the first');
end

if (samples(1,1)<1)
   error('The samples are indexed by the positive numbers');
end
clear as;

if (header<0)
   error('The "header" value cannot be negative');
end

if (nrchns<1)
   error('Number of channels myst be positive');
end

chsize=size(channels);
if (chsize(1)>1)
   samples=samples';
   chsize=chsize';
end

for i=1:chsize(2)
   if (channels(1,i)>nrchns)
      error('The input file does not contain samples for specified channels');
   end
end

nrsamples=samples(1,2)-samples(1,1)+1;

r=zeros(chsize(2),nrsamples);
%size(r)
fid0=fopen(name,'r');

a=fseek(fid0,0,1); %skok na koniec pliku
p=ftell(fid0);
size=p;

chlength=floor((p-header)/2/nrchns); %samples per channel

if (samples(1,2)>chlength)
   fclose(fid0);
   error('Too large values in "samples" argument');
end

for i=1:chsize(2)
   channel=channels(1,i);
   start=header+(channel-1)*chlength*2+(samples(1,1)-1)*2;
   fseek(fid0,start,-1);
   r(i,:)=fread(fid0,nrsamples,'int16')';
end

fclose(fid0);
