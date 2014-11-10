clear


piece = '2010-09-24-1';%'2008-08-27-5';

run = 'data003';%'data003';
path2data=['/Volumes/Analysis/' piece '/' run '/' run];
% path2data=['/Volumes/Analysis/' piece '/' run '/' run '/' run];

datarun = load_data(path2data);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);

datarun=load_sta(datarun,'load_sta','all');
datarun = load_cones(datarun);%,'standard');
% datarun = load_cones(datarun,magic_number);
datarun = make_mosaic_struct(datarun);



sta=datarun.stas.stas{2};
figure
imagesc(sta(:,:,2,6))
[a,b]=find(sta(:,:,2,6)>0.02)

m=[];
for i=1:length(a)
    m=[m squeeze(sta(a(i),b(i),:,6))];
end
m=m';
[~,ic]=sort(m(:,3))
m=m(ic,:);
figure
for i=1:3
    subplot(3,1,i)
    bar(diff(m(:,i)))
end



m=sort(m);
figure
for i=1:3
    subplot(3,1,i)
    bar(m(:,i))
end

m=sort(m);
figure
for i=1:3
    subplot(3,1,i)
    bar(diff(m(:,i)))
    axis([0 45 0 0.015])
end




figure
for i=1:3
    subplot(3,1,i)
    bar(m(:,i))
end

sum(m)

myR=sta(:,:,1,5);
myR=myR(:);
myG=sta(:,:,2,5);
myG=myG(:);
myB=sta(:,:,3,5);
myB=myB(:);
a=corr(myR,myG);
corr(myG,myB)
corr(m)

[a,b]=find(sta(:,:,2,5)>0.015);
m=[];
for i=1:length(a)
    m=[m squeeze(sta(a(i),b(i),:,5))];
end
figure
bar(m')

sum(m')


sta=datarun.stas.stas{9};

figure

for j=1:3
    
    [a,b]=find(sta(:,:,j,5)>0.032);
    m=[];
    for i=1:length(a)
        m=[m squeeze(sta(a(i),b(i),:,5))];
    end
    m=m';
    m=sort(m);
    subplot(3,1,j)
    bar(1:length(m),m(:,1)/max(m(:,j)),'r','barwidth',0.2)
    hold on
    bar(1.3:length(m)+1,m(:,2)/max(m(:,j)),'g','barwidth',0.2)
    bar(1.6:length(m)+1,m(:,3)/max(m(:,j)),'b','barwidth',0.2)
end


sta=datarun.stas.stas{1};
figure

for j=1:3
    
    [a,b]=find(sta(:,:,j,5)<-0.02);
    m=[];
    for i=1:length(a)
        m=[m squeeze(sta(a(i),b(i),:,5))];
    end
    m=m';
    m=sort(m,'descend');
    subplot(3,1,j)
    bar(1:length(m),m(:,1)/min(m(:,j)),'r','barwidth',0.2)
    hold on
    bar(1.3:length(m)+1,m(:,2)/min(m(:,j)),'g','barwidth',0.2)
    bar(1.6:length(m)+1,m(:,3)/min(m(:,j)),'b','barwidth',0.2)
end

figure
cnt=1;p1B=[];p2B=[];p1G=[];p2G=[]; p1R=[];p2R=[]; cellCl=[];
for j=get_cell_indices(datarun,{3})
    sta=datarun.stas.stas{j};
    
    subplot(3,3,cnt)
    tmp=squeeze(sta(:,:,1,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp0=tmp;
    ttt=cumsum(tmp(102300:end));
    
    forfit=double(tmp(100000:102200));
    k=fit([1:length(forfit)]',forfit,'poly1');
    p1R=[p1R k.p1];
    p2R=[p2R k.p2];
    
    plot(tmp,'.-r')
    hold on
    
    x=102200:102400;
    plot(x,k.p1*(x-100000)+k.p2,'r')
    
    tmp=squeeze(sta(:,:,2,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    ttt=cumsum(tmp(102300:end));
    
    forfit=double(tmp(100000:102200));
    k=fit([1:length(forfit)]',forfit,'poly1');
    p1G=[p1G k.p1];
    p2G=[p2G k.p2];
    
    plot(tmp,'.-g')
    
    plot(x,k.p1*(x-100000)+k.p2,'g')
    
    tmp=squeeze(sta(:,:,3,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    ttt=cumsum(tmp(102300:end));
    
    forfit=double(tmp(100000:102200));
    k=fit([1:length(forfit)]',forfit,'poly1');
    p1B=[p1B k.p1];
    p2B=[p2B k.p2];
    
    plot(tmp,'.-b')
    
    plot(x,k.p1*(x-100000)+k.p2,'b')
    axis([102300,102400,0,Inf])
    a=datarun.cell_ids(j);
    for i=1:length(datarun.cell_types)
        if sum(datarun.cell_types{1, i}.cell_ids==a)~=0
            tt=datarun.cell_types{1, i}.name;
            break
        end
    end
    

    cellCl=[cellCl i];
    title(['cell ID ',int2str(a),'  ',tt])
    cnt=cnt+1;
end

figure
plot(cellCl+0.4,p1B,'*b')
hold on
plot(cellCl,p1R,'*r')
plot(cellCl+0.2,p1G,'*g')


figure
plot(cellCl,(p1G-p1B)./(p1G+p1B),'*b')


figure
loglog(a(end:-1:1))

figure
a=diff(tmp);
a=cumsum(a);
plot(a)


figure
plot(diff(tmp1-tmp2));
axis([99800,102400,-Inf,Inf])
tt=tmp1-tmp2;
find(tt<0,1,'last')

figure
hist(diff(tmp1-tmp2),1000)



figure
cnt=1;pp=[];pm=[];myTMP=[];
for k=[1:5, 11]
    subplot(2,3,cnt)
    myTMP=[];myNND=[];
    if k==2 || k==4 
        ss=-1;
    else
        ss=1;
    end
    for j=get_cell_indices(datarun,{k})
        sta=datarun.stas.stas{j};
        
        
        tmp=squeeze(sta(:,:,1,6));
        tmp=sort(ss*tmp(:));%/max(tmp(:)));
        tmp0=tmp;
        
        tmp=squeeze(sta(:,:,2,6));
        [tmp,ic]=sort(ss*tmp(:));%/max(tmp(:)));
        tmp1=tmp;
        
        tmp=squeeze(sta(:,:,3,6));
        tmp=sort(ss*tmp(:));%/max(tmp(:)));
        tmp2=tmp;
        
%         pp=[pp max(tmp1-tmp2)];
%         pm=[pm min(tmp1-tmp2)];


        tt=tmp1-tmp2;
        myTMP=[myTMP tt];
        icc=ic(end-10:end);
        [rows,cols]=ind2sub([320,320],icc);
        a=pdist([rows, cols]);
        myNND=[myNND mean(a)];
%         plot(cumsum(tt(end-50:end)));
        hold on
    end
    plot(myNND)
%     plot(mean(myTMP,2),'r','linewidth',2)
    title(datarun.cell_types{k}.name)
    axis([0 Inf 0 25])
%     axis([-10000 112400 -0.3 0.2])
    cnt=cnt+1;
end


%% best version so far
figure
tt=[];
for k=[1,3]
    if k==1
        col='r';
    else col='b';
    end
    mySTA=zeros(320,320);
    thresh=20;
    cnt=1;
    % figure
    for j=get_cell_indices(datarun,{k})
        sta=datarun.stas.stas{j};
        
        
        tmp=squeeze(sta(:,:,1,6));
        tmp=sort(tmp(:));%/max(tmp(:)));
        tmp0=tmp;
        
        tmp=squeeze(sta(:,:,2,6));
        [tmp, ic]=sort(tmp(:));%/max(tmp(:)));
        tmp1=tmp;
        
        tmp=squeeze(sta(:,:,3,6));
        tmp=sort(tmp(:));%/max(tmp(:)));
        tmp2=tmp;
        
        my_tmp=tmp1-tmp2;
        
        %         thr=min(max([tmp(end-5) tmp2(end-5)]));
        thr=tmp2(end-10);
        myValues=find(tmp1<=thr);
        tmp=squeeze(sta(:,:,2,6));
        tmp(ic(myValues))=0;
        mySTA=mySTA+tmp;
        
        %         subplot(5,5,cnt)
        %         imagesc(tmp)
        %         cnt=cnt+1;
        
        %         plot(my_tmp,col);
        
        %         tt=[tt find(my_tmp<0,1,'last')];
        
        %         m=robust_std(my_tmp(100:end-100))*thresh;
        %         m=find(my_tmp>m);
        %         [row,column]=ind2sub([320,320],ic(m));
        %         for h=1:length(row)
        %             sta(row(h),column(h),:,6)=sta(row(h),column(h),:,6)*1000;
        %         end
        %          mySTA=mySTA+squeeze(sta(:,:,2,6));
        %          subplot(5,5,cnt)
        %          imagesc(squeeze(sta(:,:,2,6)))
        %          cnt=cnt+1;
        %          figure
        %          imagesc(tmp)
        
        %         plot(my_tmp(m))
        %         [row,col]=ind2sub([320,320],m);
        %         sta(sta<)
        %
        %         hold on
        %         line(get(gca,'XLim'),[m,m],'color','k')
    end
    
    figure
    imagesc(mySTA)
end

axis([150, 200, 100,150]*3)
%% single cone template

k=3;
halfSize=25;
cnt=1;
figure
for color_cone=1:3
    
    myCone=zeros(halfSize*2+1);
    side_count=0;    
    for j=get_cell_indices(datarun,{k})
        sta=datarun.stas.stas{j};
        tmp=squeeze(sta(:,:,color_cone,6));
        [a,b]=find(tmp==min(tmp(:)),1);
        % [x,y]=find(tmp<0);
        % m=find(x>(halfSize+100)&y>(halfSize+100),1);
        % a=x(m);b=y(m);
        
        try myCone=myCone+tmp(a-halfSize:a+halfSize,b-halfSize:b+halfSize);
        catch err
            side_count=side_count+1;
        end
    end
    myCone=myCone/(length(get_cell_indices(datarun,{k}))-side_count)*100;
    subplot(3,4,cnt)
    imagesc(myCone)
    subplot(3,4,cnt+1)
    plot(myCone(halfSize+1,:),'linewidth',2)
    hold on
    plot(myCone(:,halfSize+1),'r','linewidth',2)
    
    myCone=zeros(halfSize*2+1);
    side_count=0;
    for j=get_cell_indices(datarun,{k})
        sta=datarun.stas.stas{j};
        tmp=squeeze(sta(:,:,color_cone,6));
        [a,b]=find(tmp==max(tmp(:)),1);
        
        try myCone=myCone+tmp(a-halfSize:a+halfSize,b-halfSize:b+halfSize);
        catch err
            side_count=side_count+1;
        end
    end
    myCone=myCone/(length(get_cell_indices(datarun,{k}))-side_count)*100;
    subplot(3,4,cnt+2)
    imagesc(myCone)
    subplot(3,4,cnt+3)
    plot(myCone(halfSize+1,:),'linewidth',2)
    hold on
    plot(myCone(:,halfSize+1),'r','linewidth',2)
    
    cnt=cnt+4;
end
  
%% template matching
figure
imagesc(myCone)
figure
imagesc(mySTA)


subSize=1;
mySubCone=myCone((26-subSize):(26+subSize),(26-subSize):(26+subSize));
figure;imagesc(mySubCone)

tmp=1;
myConeList=[];
while tmp>0
    tmp=max(mySTA(:));
    [a,b]=find(mySTA==tmp);
    mySTA(a-subSize:a+subSize,b-subSize:b+subSize)=...
        mySTA(a-subSize:a+subSize,b-subSize:b+subSize)-...
        mySubCone/(min(mySubCone(:))/tmp);
%     mySubCone/(mySubCone(subSize+1,subSize+1)/tmp);
    
    myConeList=[myConeList; [a b]];
end


mySTA=zeros(320);
k=3;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    
    
    tmp=squeeze(sta(:,:,1,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp0=tmp;
    
    tmp=squeeze(sta(:,:,2,6));
    [tmp, ic]=sort(tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    
    tmp=squeeze(sta(:,:,3,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    
    my_tmp=tmp1-tmp2;
    
    %         thr=min(max([tmp(end-5) tmp2(end-5)]));
    thr=tmp2(end-10);
    myValues=find(tmp1<=thr);
    tmp=squeeze(sta(:,:,2,6));
    tmp(ic(myValues))=0;
    mySTA=mySTA+tmp;
    
end
mySTAOff=zeros(320);
k=4;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    
    
    tmp=squeeze(sta(:,:,1,6));
    tmp=sort(-tmp(:));%/max(tmp(:)));
    tmp0=tmp;
    
    tmp=squeeze(sta(:,:,2,6));
    [tmp, ic]=sort(-tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    
    tmp=squeeze(sta(:,:,3,6));
    tmp=sort(-tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    
    my_tmp=tmp1-tmp2;
    
    %         thr=min(max([tmp(end-5) tmp2(end-5)]));
    thr=tmp2(end-10);
    myValues=find(tmp1<=thr);
    tmp=-squeeze(sta(:,:,2,6));
    tmp(ic(myValues))=0;
    mySTAOff=mySTAOff+tmp;
    
end
k=3;


figure
imagesc(mySTA)
hold on
plot(myConeList(:,2),myConeList(:,1),'x','color',[1 1 1])
plot(datarun.cones.centers(:,1),datarun.cones.centers(:,2),'+y')
%individual cells
figure
colormap gray
cnt=1;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    tmp=squeeze(sta(:,:,2,6));
    subplot(3,3,cnt)
    imagesc(tmp)
    [a,b]=find(tmp==max(tmp(:)));
    hold on    
    plot(myConeList(:,2),myConeList(:,1),'xr')
    plot(datarun.cones.centers(:,1),datarun.cones.centers(:,2),'+b')
    axis([b-10 b+10 a-10 a+10])
    cnt=cnt+1;
end
   

%% scaled version (x3)

% get sta
mySTA=zeros(320*3);
k=3;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    
    
    tmp=squeeze(sta(:,:,1,6));
    tmp=imresize(tmp,3);    
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp0=tmp;
    
    tmp=squeeze(sta(:,:,2,6));
    tmp=imresize(tmp,3);
    [tmp, ic]=sort(tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    
    tmp=squeeze(sta(:,:,3,6));
    tmp=imresize(tmp,3);
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    
    my_tmp=tmp1-tmp2;
    
    %         thr=min(max([tmp(end-5) tmp2(end-5)]));
    thr=tmp2(end-30);
    myValues=find(tmp1<=thr);
    tmp=squeeze(sta(:,:,2,6));
    tmp=imresize(tmp,3);
    tmp(ic(myValues))=0;
    mySTA=mySTA+tmp;
    
end

k=4;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    
    
    tmp=-squeeze(sta(:,:,1,6));
    tmp=imresize(tmp,3);    
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp0=tmp;
    
    tmp=-squeeze(sta(:,:,2,6));
    tmp=imresize(tmp,3);
    [tmp, ic]=sort(tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    
    tmp=-squeeze(sta(:,:,3,6));
    tmp=imresize(tmp,3);
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    
    my_tmp=tmp1-tmp2;
    
    %         thr=min(max([tmp(end-5) tmp2(end-5)]));
    thr=tmp2(end-30);
    myValues=find(tmp1<=thr);
    tmp=-squeeze(sta(:,:,2,6));
    tmp=imresize(tmp,3);
    tmp(ic(myValues))=0;
    mySTA=mySTA+tmp;
    
end

k=5;
mySTA=zeros(320*3);
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    
    
    tmp=squeeze(sta(:,:,2,6));
    tmp=imresize(tmp,3);
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    
    tmp=squeeze(sta(:,:,3,6));
    tmp=imresize(tmp,3);
    [tmp, ic]=sort(tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    
    my_tmp=tmp2-tmp1;
    
    %         thr=min(max([tmp(end-5) tmp2(end-5)]));
    thr=tmp1(end-50);
    myValues=find(tmp2<=thr);
    tmp=squeeze(sta(:,:,3,6));
    tmp=imresize(tmp,3);
    tmp(ic(myValues))=0;
    mySTA=mySTA+tmp;
    
end


figure
imagesc(mySTA)

% get cone template
k=3;
halfSize=25*3;
color_cone=3;

myCone=zeros(halfSize*2+1);
side_count=0;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    tmp=squeeze(sta(:,:,color_cone,6));
    tmp=imresize(tmp,3);
    [a,b]=find(tmp==max(tmp(:)),1);
    
    try myCone=myCone+tmp(a-halfSize:a+halfSize,b-halfSize:b+halfSize);
    catch err
        side_count=side_count+1;
    end
end
myCone=myCone/(length(get_cell_indices(datarun,{k}))-side_count)*100;
figure
imagesc(myCone)

subSize=4;
mySubCone=myCone((76-subSize):(76+subSize),(76-subSize):(76+subSize));
figure;imagesc(mySubCone)

% find cones
tmp=1;
myConeList=[];
mySTAsaved=mySTA;
% mySTA=mySTAsaved;
while tmp>0
    tmp=max(mySTA(:));
    [a,b]=find(mySTA==tmp);
    mySTA(a-subSize:a+subSize,b-subSize:b+subSize)=...
        mySTA(a-subSize:a+subSize,b-subSize:b+subSize)-...
        mySubCone/(min(mySubCone(:))/tmp);
%     mySubCone/(mySubCone(subSize+1,subSize+1)/tmp);
    
    myConeList=[myConeList; [a b]];
end


figure
imagesc(mySTAsaved)
hold on
plot(myConeList(:,2),myConeList(:,1),'x','color',[1 1 1])

%individual cells
figure
colormap gray
cnt=1;
for j=get_cell_indices(datarun,{k})
    sta=datarun.stas.stas{j};
    tmp=squeeze(sta(:,:,3,6));
    tmp=imresize(tmp,3);
    subplot(3,3,cnt)
    imagesc(tmp)
    [a,b]=find(tmp==max(tmp(:)));
    hold on    
    plot(myConeList(:,2),myConeList(:,1),'xr')
    plot(datarun.cones.centers(:,1)*3,datarun.cones.centers(:,2)*3,'+b')
    axis([b-30 b+30 a-30 a+30])
    cnt=cnt+1;
end
   
%% 
mySTA=zeros(320,320);
thresh=20;
cnt=1;
myCells=get_cell_indices(datarun,{k});
figure
for j=myCells(1:5:35)
    sta=datarun.stas.stas{j};
    
    
    tmp=squeeze(sta(:,:,1,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp0=tmp;
    
    tmp=squeeze(sta(:,:,2,6));
    [tmp, ic]=sort(tmp(:));%/max(tmp(:)));
    tmp1=tmp;
    
    tmp=squeeze(sta(:,:,3,6));
    tmp=sort(tmp(:));%/max(tmp(:)));
    tmp2=tmp;
    
    my_tmp=tmp1-tmp2;
    
    m=robust_std(my_tmp(100:end-100))*thresh;
    mind=find(my_tmp>m);
    [row,column]=ind2sub([320,320],ic(mind));
    for h=1:length(row)
        sta(row(h),column(h),:,6)=sta(row(h),column(h),:,6)*1000;
    end
    
    m=-robust_std(my_tmp(100:end-100))*thresh;
    mind=find(my_tmp<m);
    [row,column]=ind2sub([320,320],ic(mind));
    for h=1:length(row)
        sta(row(h),column(h),:,6)=-sta(row(h),column(h),:,6)*1000;
    end
    
    mySTA=mySTA+squeeze(sta(:,:,2,6));
    subplot(5,5,cnt)
    imagesc(squeeze(sta(:,:,2,6)))
    cnt=cnt+1;
end
figure
imagesc(mySTA)


figure
plot(my_tmp)
