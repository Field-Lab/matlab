function y=stim_analyze(filename,dlugosc,okres,stim_channels,figura,marg,zakres,plots);

wspolrzedne=zeros(61,2);

wspolrzedne(1:10,:)=[[5,1]',[4,2]',[6,2]',[3,3]',[5,3]',[7,3]',[2,4]',[4,4]',[6,4]',[8,4]']';
wspolrzedne(11:15,:)=[[1,5]',[3,5]',[5,5]',[7,5]',[9,5]']';
wspolrzedne(16:19,:)=[[2,6]',[4,6]',[6,6]',[8,6]']';

wspolrzedne(20:24,1)=wspolrzedne(11:15,1);
wspolrzedne(20:24,2)=wspolrzedne(11:15,2)+2;
wspolrzedne(25:28,1)=wspolrzedne(16:19,1);
wspolrzedne(25:28,2)=wspolrzedne(16:19,2)+2;

wspolrzedne(29:33,1)=wspolrzedne(11:15,1);
wspolrzedne(29:33,2)=wspolrzedne(11:15,2)+4;
wspolrzedne(34:37,1)=wspolrzedne(16:19,1);
wspolrzedne(34:37,2)=wspolrzedne(16:19,2)+4;

wspolrzedne(38:42,1)=wspolrzedne(11:15,1);
wspolrzedne(38:42,2)=wspolrzedne(11:15,2)+6;
wspolrzedne(43:46,1)=wspolrzedne(16:19,1);
wspolrzedne(43:46,2)=wspolrzedne(16:19,2)+6;

wspolrzedne(47:51,1)=wspolrzedne(11:15,1);
wspolrzedne(47:51,2)=wspolrzedne(11:15,2)+8;
wspolrzedne(52:55,1)=wspolrzedne(16:19,1);
wspolrzedne(52:55,2)=wspolrzedne(16:19,2)+8;

wspolrzedne(56:61,:)=[[3,15]',[5,15]',[7,15]',[4,16]',[6,16]',[5,17]']';
    
numery(1,1:19)=[30 31 27 34 29 23 37 32 24 20 40 35 28 21 19 39 33 22 18];
numery(1,20:37)=[42 38 26 17 16 43 36 15 14 45 44 41 12 13 46 47 4 11];
numery(1,38:61)=[48 49 58 6 10 50 54 1 7 51 53 60 3 8 52 56 64 5 55 61 2 59 63 62];

xstart=0.05;
ystart=0.02;
xkrok=0.09;
ykrok=0.05;
xszer=0.075;
yszer=0.075;

start=1;
margines=1;

%time=[1:okres]/20000*1000;
time=[1:dlugosc]/20000*1000;
marg=marg*20000/1000;

h=figure(figura);
clf;

%s=readconv(filename,206,65,1,[(start) (start+dlugosc-1)]);
%smin=min(min(s));
%smax=max(max(s));
%detect=find(s<(smin+smax)/2-200);
%poczatek=detect(1,1);
%if poczatek<marg+1
%    poczatek=poczatek+okres;
%end

%s=readconv(filename,206,65,1,[(start) (start+period-1)]);
%(det0,det1)=min(s)
 
poczatek=0;

margines=poczatek-marg

margines=marg;
            
x0=0;
x1=x0+zakres;
y0=-1200;
y1=1200;

for i=0:61
    if i~=0
        s=readconv(filename,206,65,numery(1,i)+1,[(start) (start+dlugosc-1)]);
    else
        s=readconv(filename,206,65,1,[(start) (start+dlugosc-1)]);
    end
        
    
    s=s(margines:dlugosc);
    ss=size(s);
    s_dl=ss(2)
    time=[1:s_dl]/20000*1000;
    
    n_okr=floor(s_dl/okres);
    s_usr=zeros(1,okres);
    
    for j=1:n_okr
        s0=s(1+(j-1)*okres:j*okres);
        s_usr=s_usr+s0;
    end
    s_usr=s_usr/n_okr;
    
    if i~=0
        subplot('position',[wspolrzedne(i,1)*xkrok+xstart,wspolrzedne(i,2)*ykrok+ystart,xszer,yszer]);
        if numery(1,i)<33                       
            q=plot(time,s,'b-');
        else
            q=plot(time,s,'r-');
        end
        %size(s)
        
        if find(stim_channels==numery(1,i))
            set(q,'LineWidth',2);
        end
        
        axis([x0 x1 y0 y1]);
        h=gca;
        grid on;
        num2str(numery(1,i));
        q=title(num2str(numery(1,i)));
        set(q,'Position',[x0+(x1-x0)*0.84,y0+(y1-y0)*0.7]);
        if numery(1,i)~=40
            set(h,'XTickLabel',[]);
            set(h,'YTickLabel',[]);
        else
            xlabel('time [msec]');
            ylabel('[DAC units]');
        end
        
    else
        subplot('position',[0.14,0.87,xszer,yszer]);
        size(s)
        plot(time,s,'k-');
        
        axis([x0 x1 -2500 -1000]);
        grid on;
        xlabel('time [msec]');
        title('Vstim')
    end
   
end

if plots
    cd obrazki;
    namesize=size(filename);
    obr_name=[filename(1,1:(namesize(2)-4)) '_' num2str(zakres) 'ms']
    print(figura,'-dbmp',obr_name);
    print(figura,'-depsc',obr_name);
    print(figura,'-dpng',obr_name);
    cd .. ;
    y=obr_name;
else
    y='no graphic files';
end