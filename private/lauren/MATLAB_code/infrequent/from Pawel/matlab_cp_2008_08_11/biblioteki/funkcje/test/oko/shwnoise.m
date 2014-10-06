function w=shwnoise(name,channels,dlugosc,start,N,M,fp,lsb,figura1,figura2,figura3);


header=206;
nrchns=65;

w=zeros(length(start),N);

%figure(figura1);
%clf(figura1);
%figure(figura2);
%clf(figura2);
figure(figura3);
%clf(figura3);

kolory=['m','b','g','r','k'];
lkolory=length(kolory)

lstart=length(channels);

p=sqrt(lstart);

a2=ceil(p);
%a2=20;

if a2==0
   a2=1;
end

a1=ceil(lstart/a2);

dl=ceil(N/2);
%szum=rand(1,dlugosc);

for i=1:lstart
    
    %subplot(a1,a2,i);
    for j=1:length(start)
        r=readcnst2(name,header,nrchns,channels(i),[start(j) start(j)+dlugosc-1]);
        [f,w(j,:)]=spdf(r,N,M,fp,lsb);

        f0=f(1,1:dl);
        
        s=j-floor(j/lkolory)*lkolory+1
        kolor=kolory(j-floor(j/lkolory)*lkolory+1)
        
        %figure(figura1);
        %subplot(a1,a2,i); %wlasnie i!
        %hold on;
        %plot(r,kolor);
        %grid on;
        %mr(j)=max(abs(r));
        
        %figure(figura2);
        %hold on;
        %subplot(a1,a2,i);
        %plot(f0,w(j,1:dl),kolor);
        %axis([0 f0(1,dl) 0 0.5]);
        %grid on;
        %mw=max(abs(w(j,dl)));
        
        figure(figura3);
        hold on;
        subplot(3,2,i+4);
        loglog(f0,w(j,1:dl),kolor);
        hand=gca;
        axis([0 f0(1,dl) 1e-4 1]);
        set(hand,'XTick',[100 1000]);
        %set(hand,'XTickLabel',[50] 5000]);
        grid on;
        
        %plot(f,w(j,:));
    end
    %mr0=max(mr);
    %mw0=max(mw);
    %figure(1);
    %plot(w(1,:)
    %figure(figure2);
    
    %plot(f0,w(1,1:dl),f0,w(2,1:dl),f0,w(3,1:dl));
    %grid on;
    %hold off;
end 