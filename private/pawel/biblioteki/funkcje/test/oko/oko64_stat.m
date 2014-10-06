function w=oko64_stat(signal,clock,nrchns);

ls=length(signal);
lc=length(clock);

if ls~=lc
    error('Length of signal must be equal to length of clock');
end

disc=mean(clock) %na razie !!!!
disc=2.5;

sigcl=(sign(clock-disc));

%s=find(sigcl==1);
s=sigcl+abs(sigcl);
ds=diff(s);

p=find(ds==2);
size(p)
dlugosc=floor(length(p)/nrchns)

chns=zeros(nrchns,dlugosc);

for i=1:nrchns
    %i:nrchns:dlugosc*nrchns
    chns(i,1:dlugosc)=signal(p(i:nrchns:dlugosc*nrchns));
end

w=chns;