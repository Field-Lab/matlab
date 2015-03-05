nds='76543212345'
figure

cd('/mnt/muench_data/user/alexandra/scripts')
pp=1;
for j=1:3
    date=['20120829_',int2str(j)];
    
    path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
    keyWord='ffflicker';
    lfpFilters=dir([path2data,'*',keyWord,'*.mat']);
    
    all_ch=0;
    for ch=1:59
        load([path2data,lfpFilters(ch).name])
        all_ch=all_ch+HighFilter;
    end
    
    cnt=1;
    for i=1:4:25
        subplot(3,7,pp)
        plot(all_ch(1:500,i:i+3)/59,'LineWidth',2)
        axis([0,500,-15,10])
        title(['ND',nds(cnt),'  ',date,'  OPN4 LFP'])
        cnt=cnt+1;
        pp=pp+1;
    end
end



nds='76543212345'
figure

cd('/mnt/muench_data/user/alexandra/scripts')
pp=1;
for j=1:2
    if j==1
        date='20121002';
        l=12;
    else
        date='20120928';     
        l=0;
    end
    
    path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
    keyWord='HCseq';
    lfpFilters=dir([path2data,'*',keyWord,'*.mat']);
    
    all_ch=0;
    for ch=1:59
        load([path2data,lfpFilters(ch).name])
        all_ch=all_ch+HighFilter;
    end
    
    cnt=1;
    for i=1:12:84
        subplot(2,7,pp)
        plot(all_ch(1:500,i+l:2:i+11+l)/59,'LineWidth',2)
        axis([0,500,-60,30])
        title(['ND',nds(cnt),'  ',date,'  CNGa3 LFP'])
        cnt=cnt+1;
        pp=pp+1;
    end
end



date='20121004';

path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
keyWord='HCseq';
lfpFilters=dir([path2data,'*',keyWord,'*.mat']);

all_ch=0;
for ch=1:59
    load([path2data,lfpFilters(ch).name])
    all_ch=all_ch+HighFilter;
end

t=[12 30:9:71 71 72:12:194];
tot_high=[1 25 9 9 9 9 9 1 12 12 12 12 12 9 9 9 9 9 9 9];
a=cumsum(tot_high);
a(end)=a(end)-1;
nds='9876543432123456783'
figure
cnt=1;
for i=1:19
    subplot(4,5,cnt)
    plot(all_ch(1:500,a(i):a(i+1)-1)/59,'LineWidth',2)

    axis([0 500 -60 20])
title(['ND',nds(cnt)])
cnt=cnt+1;
end









cd('/mnt/muench_data/user/alexandra/scripts')
figure

for j=1:5
    if j==1
        date='20121002';
        tt=49:60;
        n=30;
        kk=tt-24;
    elseif j==2
        date='20120928';     
        tt=37:48;
        n=18;
        kk=tt-24;
    elseif j==3
        date='20121017';
        tt=37:48;
        n=18;
        kk=tt-24;
    elseif j==4
        date='20121018';
        tt=49:60;
        n=30;
        kk=tt-24;
    else
        date='20121004';
        tt=62:70;
        n=45;
        kk=tt-18;
    end
    
    path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
    keyWord='HCseq';
    lfpFilters=dir([path2data,'*',keyWord,'*.mat']); 

    all_ch=0;
    for ch=1:59
        load([path2data,lfpFilters(ch).name])
        all_ch=all_ch+HighFilter;
    end
    a=min(all_ch(:,n));
    
    subplot(5,2,j*2-1)
    plot(-all_ch(:,kk)/a,'LineWidth',2)
    title(['ND6', date, ' normalized to ND6 middle min'])
    axis([0,500,-1.5,0.2])    
    
    subplot(5,2,j*2)
    plot(-all_ch(:,tt)/a,'LineWidth',2)
    title(['ND4', date, ' normalized to ND6 middle min'])
    axis([0,500,-1,0.2])
end






nds='765432123'
cd('/mnt/muench_data/user/alexandra/scripts')
pp=1;
figure
for j=1:5
    if j==1
        date='20121002';
        n=12;
        kk='cng3a';
    elseif j==2
        date='20120928';     
        n=0;
        kk='cng3a';
    elseif j==3
        date='20121017';
        n=0;
        kk='cpfl1';
    elseif j==4
        date='20121018';
        n=12;
        kk='cpfl1';
    else
        date='20121004';
        n=45;
        kk='cng3a';
    end
    
    path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
    keyWord='HCseq';
    lfpFilters=dir([path2data,'*',keyWord,'*.mat']);
    
    all_ch=0;
    for ch=1:59
        load([path2data,lfpFilters(ch).name])
        all_ch=all_ch+HighFilter;
    end
    all_ch=-all_ch/min(all_ch(:));
    a(1)=max(all_ch(:));
    a(2)=min(all_ch(:));
    
    
    if j<5
        cnt=1;
        for i=1:12:108
            subplot(5,9,pp)
            plot(all_ch(1:500,i+n:i+11+n),'LineWidth',2)
            %         axis([0,500,-15,10])
            title(['ND',nds(cnt),' ',date,' ',kk,' LFP'])
            cnt=cnt+1;
            pp=pp+1;
            axis([0 500 a(2) a(1)])
        end
    else
        cnt=1;
        for i=1:9:36
            subplot(5,9,pp)
            plot(all_ch(1:500,34+i:34+i+8),'LineWidth',2)
            %         axis([0,500,-15,10])
            title(['ND',nds(cnt),' ',date,' ',kk,' LFP'])
            cnt=cnt+1;
            pp=pp+1;
            axis([0 500 a(2) a(1)])
        end
    end
end