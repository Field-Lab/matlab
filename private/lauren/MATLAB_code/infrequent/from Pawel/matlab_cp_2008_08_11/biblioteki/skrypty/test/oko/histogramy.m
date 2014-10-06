cd /mnt/win3/oko/2000-12-12/;

clear;
name='Data013';
figure(4);
%chns=[6 17 28 29 31 34 59 65];
%channels=[10 12 14 17 20:34 36:46 51 55 58:65];
%chns=[2:65];
%chns=22;
%chns=[6 10 12 14 17 20]
%chns=[26 27];
%channels=[24 28 29 31];
%channels=[55 59 62 65];
channels=[55 59 62 64];
%channels=64;
%channels=[34 38 40 42];
%channels=55;
%channels=[29];
%channels=[10 12 14 17];
ampl=[150 150 150];
start=1;
%for i=1:3
%25000000
    name='Data009';
    %w0=statpks2(name,chns,2000000,[start start+10000000-1],50,40,1);
    name='Data011';
    %w0=statpks2(name,chns,2000000,[start start+10000000-1],50,40,2);
    name='Data013';
    %w0=statpks2(name,chns,2000000,[start start+10000000-1],50,40,3);
%chns=[17 26 31 51 65];
%chns=[17 65];
osie=[1000 1000 500 500 200 2000];
rozdz=200;
prog=50;
histereza=15;

%for j1=1:5
 %   j1    
  %  chns=channels(1,(8*(j1-1)+1):(8*j1))
   % size(chns);
    %figure(j1);
    chns=channels;
lchns=length(channels)
for i=4:lchns    
    i
    name='Data009';
    w1=statpks2(name,chns(i),2000000,[start start+10000000-1],prog,histereza,4);
    %subplot(3,lchns,i);
    min1=min(w1(3,:));
    max1=max(w1(3,:));
    roz1=size(w1);
    h1=hist(w1(3,:),max1-min1);
    
    name='Data011';
    w2=statpks2(name,chns(i),2000000,[start start+10000000-1],prog,histereza,5);
    %subplot(3,lchns,i+lchns);
    %hist(w0(3,:));
    min2=min(w2(3,:));
    max2=max(w2(3,:));
    roz2=size(w2);
    h2=hist(w2(3,:),max2-min2);
    
    name='Data013';
    w3=statpks2(name,chns(i),2000000,[start start+10000000-1],prog,histereza,6);
    %subplot(3,lchns,i+2*lchns);
    %hist(w0(3,:));
    min3=min(w3(3,:));
    max3=max(w3(3,:));
    roz3=size(w3);
    h3=hist(w3(3,:),max3-min3);
    
    min0=min([min1 min2 min3])-200;
    max0=max([max1 max2 max3])+200;
    h0=max([h1 h2 h3])*1.2;            
    
    subplot(3,lchns,i);
    %plot(h1);
    hist(w1(3,:),max1-min1);
    axis([min0 max0 0 h0]);
    grid on;
    
    subplot(3,lchns,i+lchns);
    %plot(h2);
    hist(w2(3,:),max2-min2);
    axis([min0 max0 0 h0]);
    grid on;
    
    subplot(3,lchns,i+2*lchns);
    %plot(h3);
    hist(w3(3,:),max3-min3);
    axis([min0 max0 0 h0]);
    grid on;
      
end
end
    a='blabla'