function h=statpks(name,channels,block,samples,prog,histereza,figura);
%block - rozmiar maksymalnego bloku probek, na ktorym nalezy 
%operowac - np. przy funkcji readcnst - dla unikniecia 
%"swapowania" po dysku. Typowo: 2000000;
figure(figura);
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
    i
    w=0;
    pks=zeros(2,0); %inicjacja
    for j=1:nrloops
        %data(1,((j-1)*chlength+1):(j*chlength));
        start=(j-1)*block+samples(1);
        stop=j*block+samples(1)-1;
        s=readcnst3(name,206,65,channels(i),[start stop]);
        p=detekcja5(s,prog,histereza);
        if p~=0
            l='grrrr';
            p=p+start-1;
            pks=[pks p];
        end
        %pks;
    end
    
    if nrloops*block<chlength
        s=readcnst3(name,206,65,channels(i),[nrloops*block+samples(1) samples(2)]);
        %size(s)
        p=detekcja5(s,prog,histereza);
        if p~=0
            l='uama';
            p=p+nrloops*block+samples(1)-1;
            pks=[pks p];
        end
        %whos
    end
    if pks~=0 
        spks=size(pks)
    else
        pks;
        spks=[0 0];
    end
    
    %size(spks)
    
    for j=1:spks(2)
        %j;
        %a=pks(:,j);
        s=readcnst3(name,206,65,channels(i),[pks(1,j) pks(2,j)]);
        w(j)=max(abs(s))*sign(s(1,1));
    end
    
    w;
    subplot(3,4,8+i);
    hist(w,100);
end

h=w;
%h=pks(1,:);

    
        
        




