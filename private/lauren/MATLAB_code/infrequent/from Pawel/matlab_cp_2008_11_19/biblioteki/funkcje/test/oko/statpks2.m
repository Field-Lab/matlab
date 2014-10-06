function h=statpks2(name,channels,block,samples,prog,histereza,figura);
%block - rozmiar maksymalnego bloku probek, na ktorym nalezy 
%operowac - np. przy funkcji readcnst - dla unikniecia 
%"swapowania" po dysku. Typowo: 2000000;
%figure(figura);
%clf(figura);
samples;

nrchannels=length(channels);

p0=sqrt(nrchannels);

a2=floor(p0);
%a2=20;

if a2==0
   a2=1;
end

a1=ceil(nrchannels/a2);

chlength=samples(2)-samples(1)+1;
%data=zeros(1,chlength);

nrloops=floor(chlength/block);
w=0;

for i=1:nrchannels
    i;
    w=0;
    pks=zeros(3,0); %inicjacja
    for j=1:nrloops
        %data(1,((j-1)*chlength+1):(j*chlength));
        start=(j-1)*block+samples(1)
        stop=j*block+samples(1)-1
        s=readcnst3(name,206,65,channels(i),[start stop]);
        p=detekcja6(s,prog,histereza);
        if p~=0
            l='grrrr';
            p(1,:)=p(1,:)+start-1;
            p(2,:)=p(2,:)+start-1;
            pks=[pks p];
        end
        %pks;
    end
    
    if nrloops*block<chlength
        s=readcnst3(name,206,65,channels(i),[nrloops*block+samples(1) samples(2)]);
        %size(s)
        p=detekcja6(s,prog,histereza);
        if p~=0
            l='uama';
            p(1,:)=p(1,:)+nrloops*block+samples(1)-1;
            p(2,:)=p(2,:)+nrloops*block+samples(1)-1;
            pks=[pks p];
        end
        %whos
    end
    size(pks);
    
    %w;
    %subplot(a1,a2,i);
    %hist(pks(3,:),100);
end

h=pks;
%h=pks(1,:);

    
        
        




