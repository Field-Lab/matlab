clear

%% prepare parameters and HEKA files
date='20140121b';
codeWord='cbSparse';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
heka=dir([mainpath,'HEKA/*.phys']);

%find nd changes
nd_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'nd%nd', 'once'))
        nd_list(1,end+1)=i;        
        a=regexp(heka(i).name,'%nd_');
        a=heka(i).name(a+4);
        nd_list(2,end)=str2num(a);
    end
end

% select heka files
ndFullList=[];
file_list=[];
cnt=1;
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];        
        ndFullList=[ndFullList nd_list(2,find(nd_list(1,:)<file_list(end),1,'last'))];        
    end
end


%% extract responses to patterns

units=dir([mainpath,'units/*.mat']);

cbSparse=zeros(4500,length(file_list)*8,length(units)); % response matrix: (response, patterns, units)
names=cell(length(units),1);

for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    trialNumber=1;
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,41500);
        conv=conv(121:end-120);
        
        protocol=read_header_field_heka([mainpath,'HEKA/'], heka(file_list(i)).name, 'Stimulus Protocol');
        protocol=round(protocol(2:end-1,[1,2]));

        for j=1:2:size(protocol,1)-1
            cbSparse(:,trialNumber,cnt)=conv(protocol(j,1)-499:protocol(j,1)+4000);
            trialNumber=trialNumber+1;
        end
        
  
    end  
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end



%% save stuff out

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/analysis/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

save([path2save,'cbSparse'],'ndFullList','cbSparse','names')


%% find cells locations from sta

stapath='/Users/alexth/Desktop/old_stuff/ONresp/20140121b/analysis/sta/';
stas=dir([stapath,'sta_*.mat']);
sta_lists{1}=[1 2]; %ND8
sta_lists{2}=[3 4];%ND7
sta_lists{3}=[5 6 7];%ND6
sta_lists{4}=[8 9 10];%ND5
sta_lists{5}=[11];%ND4
sta_all=zeros(500,1600,14,5);
for i=1:length(sta_lists)
    tmp=0;
    for j=1:numel(sta_lists{i})
        load([stapath,stas(sta_lists{i}(j)).name])
        tmp=tmp+sta;
    end
    sta_all(:,:,:,i)=tmp;
end
sta=sta_all; clear sta_all stas tmp sta_lists stapath

% plot all RFs at ND6 and find 'center' stixel
nd=6; % nd you want
nds=[8 7 6 5 4]; % all nd list

figure 
colormap gray
locs=zeros(14,2);
for j=1:14
    subplot(3,5,j)
    m=reshape(std(sta(:,:,j,nds==nd)),40,40); 
%     m=rot90(m,2); %flipped twice to align with electrodes
    imagesc(m)
    
    [~,loc]=max(m(:));

    [x,y]=ind2sub([40,40],loc);
    hold on
    plot(y,x,'or')
    locs(j,:)=[y,x];
    title(names{j},'interpreter','none')
end


%check threshold for RF (u can't run this part)
figure % ND6
for j=1:14
    subplot(3,5,j)
    m=std(sta(:,:,j,3)); 
    [ai,ic]=sort(m,'descend');
    plot(ai)
    r=robust_mean(m)+5*robust_std(m);% threshold! 
    line([0 1600],[r r],'color','r')
    axis tight
end



%% plot both grids

% locations of binary checkerboards!
bin_beg=((1:40)-21)*20;
bin_end=bin_beg+19;

% locations of sparse checkerboards!
sp_beg=((1:20)-11)*35;
sp_end=sp_beg+34;

figure
% plot binary
for j=1:40
    line([-400,400],[bin_beg(j),bin_beg(j)])   
end
line([-400,400],[bin_end(j),bin_end(j)])
for j=1:40
    line([bin_beg(j),bin_beg(j)],[-400,400])
end
line([bin_end(j),bin_end(j)],[-400,400])

%plot sparse
for j=1:20
    line([-350,350],[sp_beg(j),sp_beg(j)],'color','m')
end
line([-350,350],[sp_end(j),sp_end(j)],'color','m')
for j=1:20
    line([sp_beg(j),sp_beg(j)],[-350,350],'color','m')
end
line([sp_end(j),sp_end(j)],[-350,350],'color','m')


%% RF in full screen version - unitwise

unit=6;

m=reshape(std(sta(:,:,unit,3)),40,40);

bin_full=zeros(800,800);
beg_coord=(0:39)*20+1;
for i=1:40
    for j=1:40
        bin_full(beg_coord(i):beg_coord(i)+19,beg_coord(j):beg_coord(j)+19)=m(i,j);
    end
end
%normalize
bin_full=bin_full/(max(bin_full(:)));

figure;
colormap gray
imagesc(bin_full)



%% all sparse full matrices

% make full sparse CB matrices (equal in size to binary checkerboard in full)

load('/Users/alexth/Desktop/old_stuff/stim_new/matr_assigned_35_4_700.mat','matr')
matr(:,:,1:4)=[];
matr_tmp=zeros(20,20,357);
matr_tmp(:,:,1:34)=matr(:,:,1:34);
matr_tmp(:,:,35)=zeros(20,20);
matr_tmp(:,:,36:end)=matr(:,:,35:end);
matr=matr_tmp;
clear matr_tmp


all_sp_full=zeros(800,800,size(matr,3));
beg_coord=(0:19)*35+1+50;
for k=1:size(matr,3)
    for i=1:size(matr,1)
        for j=1:size(matr,2)
            all_sp_full(beg_coord(i):beg_coord(i)+34,beg_coord(j):beg_coord(j)+34,k)=matr(i,j,k);
        end
    end
end

% making positive and negative matrices
all_pos_sp_full=all_sp_full;
all_pos_sp_full(all_sp_full<0)=0;
all_neg_sp_full=-all_sp_full;
all_neg_sp_full(all_sp_full>0)=0;

save('/Users/alexth/Desktop/old_stuff/ONresp/20140121b/analysis/all_sparse_full.mat','all_sp_full','all_pos_sp_full','all_neg_sp_full')

%% make combinations

comb=zeros(800,800,3);
comb(:,:,1)=bin_full;
comb(:,:,2)=all_pos_sp_full(:,:,100);
comb(:,:,3)=all_neg_sp_full(:,:,100);

figure
imagesc(comb)

%% auxillary debugging
cbSparse=cbSparse(:,361:720,6);
figure
cnt=1;
cor=8*4;
for i=1:8
    subplot(8,2,cnt)
    comb=zeros(800,800,3);
    comb(:,:,1)=bin_full;
    comb(:,:,2)=all_pos_sp_full(:,:,i+cor);
    comb(:,:,3)=all_neg_sp_full(:,:,i+cor);
    imagesc(comb)
    
     subplot(8,2,cnt+1)
     plot(cbSparse(:,i+cor))
     title(['pattern',int2str(i+cor)])
     axis([0 4500 0 70])
     cnt=cnt+2;
end


[ai,ic]=sort(max(cbSparse(500:820,:)));
maxpos=ic(ai>80);

figure
plot(cbSparse(:,maxpos))
hold on
plot(mean(cbSparse(:,maxpos),2),'k','linewidth',3)


[ai,ic]=sort(max(cbSparse((500:820)+2000,:)));
maxpos=ic(ai>80);

figure
plot(cbSparse(:,maxpos))
hold on
plot(mean(cbSparse(:,maxpos),2),'k','linewidth',3)


figure
cnt=1;
for i=1:6
    subplot(6,2,cnt)
    comb=zeros(800,800,3);
    comb(:,:,1)=bin_full;
    comb(:,:,2)=all_pos_sp_full(:,:,maxpos(i));
    comb(:,:,3)=all_neg_sp_full(:,:,maxpos(i));
    imagesc(comb)
    subplot(6,2,cnt+1)
    plot(cbSparse(:,maxpos(i)))
    axis([0 4500 0 100])
    cnt=cnt+2;
end

%% loop through cells
load('/Users/alexth/Desktop/old_stuff/ONresp/20140121b/analysis/quick.mat')
white_flash=squeeze(mean(white_flash(:,ndFullList==6,:),2));
black_flash=squeeze(mean(black_flash(:,ndFullList==6,:),2));

thresh=5;

for unit=1:14;
    
    m=reshape(std(sta(:,:,unit,3)),40,40);
    
    bin_full=zeros(800,800);
    beg_coord=(0:39)*20+1;
    for i=1:40
        for j=1:40
            bin_full(beg_coord(i):beg_coord(i)+19,beg_coord(j):beg_coord(j)+19)=m(i,j);
        end
    end
    
    bin_full=bin_full/(max(bin_full(:)));
    threshold=robust_mean(bin_full(:))+thresh*robust_std(bin_full(:));
    
    rf_parts=bin_full;
    rf_parts(rf_parts<threshold)=0;
    
    
    % positive and negative trials
    pos_trials=[];
    neg_trials=[];
    c=zeros(800,800,2);
    d=c;
    for i=1:size(all_pos_sp_full,3)
        tmp=rf_parts.*all_pos_sp_full(:,:,i);
        if nnz(tmp)>0
            pos_trials=[pos_trials i];
            c(:,:,1)=c(:,:,1)+all_pos_sp_full(:,:,i);
            c(:,:,2)=c(:,:,2)+all_neg_sp_full(:,:,i);
        end
        
        tmp=rf_parts.*all_neg_sp_full(:,:,i);
        if nnz(tmp)>0
            neg_trials=[neg_trials i];
            d(:,:,1)=d(:,:,1)+all_pos_sp_full(:,:,i);
            d(:,:,2)=d(:,:,2)+all_neg_sp_full(:,:,i);
        end        
    end
    
    
    figure
    set(gcf,'position',[1807  1 755 951]);
    
    trials=setdiff(pos_trials,neg_trials);
    data=cbSparse(1:4500,trials+360,unit);
    subplot(4,2,1)
    plot(mean(data,2))
    axis tight
    a=get(gca,'YLim');
    line([500,500],a,'color','k')
    line([2500,2500],a,'color','k')
    title([names{unit},',  positive'],'interpreter','none')
    
    trials=setdiff(neg_trials,pos_trials);
    data=cbSparse(1:4500,trials+360,unit);
    subplot(4,2,3)
    plot(mean(data,2))
    axis tight
    a=get(gca,'YLim');
    line([500,500],a,'color','k')
    line([2500,2500],a,'color','k')
    title([names{unit},',  negative'],'interpreter','none')
    
    trials=1:360;
    trials=setdiff(trials,[pos_trials neg_trials]);
    data=cbSparse(1:4500,trials+360,unit);
    subplot(4,2,2)
    plot(mean(data,2))
    axis tight
    a=get(gca,'YLim');
    line([500,500],a,'color','k')
    line([2500,2500],a,'color','k')
    title([names{unit},',  rest'],'interpreter','none')
    
   
    subplot(4,2,5)
    plot(white_flash(:,unit))
    axis tight
    a=get(gca,'YLim');
    line([500,500],a,'color','k')
    line([2500,2500],a,'color','k')
    title([names{unit},',  white flash FF'],'interpreter','none')
    
    subplot(4,2,7)
    plot(black_flash(:,unit))
    axis tight
    a=get(gca,'YLim');
    line([500,500],a,'color','k')
    line([2500,2500],a,'color','k')
    title([names{unit},',  black flash FF'],'interpreter','none')
    
 
    c=c/max(c(:));
    d=d/max(d(:));
    comb=zeros(800,800,3);
    comb(:,:,1)=bin_full;
    
    subplot(4,2,4)
    imagesc(comb)
    title('original RF')
    
    comb(:,:,2:3)=c;

    subplot(4,2,6)
    imagesc(comb)
    title(['positive stimulus map, ',int2str(length(unique(pos_trials))),' trials'])
    
    comb(:,:,2:3)=d;
    subplot(4,2,8)
    imagesc(comb)
    title(['negative stimulus map, ',int2str(length(unique(neg_trials))),' trials'])
    
    saveas(gcf,['/Users/alexth/Desktop/old_stuff/ONresp/20140121b/cbSparse_plot/',names{unit},'.png'])
    
end

% version with overlay

thresh=20;

for unit=1:14;
    
    m=reshape(std(sta(:,:,unit,3)),40,40);
    
    bin_full=zeros(800,800);
    beg_coord=(0:39)*20+1;
    for i=1:40
        for j=1:40
            bin_full(beg_coord(i):beg_coord(i)+19,beg_coord(j):beg_coord(j)+19)=m(i,j);
        end
    end
    
    bin_full=bin_full/(max(bin_full(:)));
    threshold=robust_mean(bin_full(:))+thresh*robust_std(bin_full(:));
    
    rf_parts=bin_full;
    rf_parts(rf_parts<threshold)=0;
    
    
    % positive and negative trials
    pos_trials=[];
    neg_trials=[];
    c=zeros(800,800,2);
    d=c;
    for i=1:size(all_pos_sp_full,3)
        tmp=rf_parts.*all_pos_sp_full(:,:,i);
        if nnz(tmp)>0
            pos_trials=[pos_trials i];
            c(:,:,1)=c(:,:,1)+all_pos_sp_full(:,:,i);
            c(:,:,2)=c(:,:,2)+all_neg_sp_full(:,:,i);
        end
        
        tmp=rf_parts.*all_neg_sp_full(:,:,i);
        if nnz(tmp)>0
            neg_trials=[neg_trials i];
            d(:,:,1)=d(:,:,1)+all_pos_sp_full(:,:,i);
            d(:,:,2)=d(:,:,2)+all_neg_sp_full(:,:,i);
        end        
    end
    
    
    figure
    set(gcf,'position',[1807  1 755 951]);
    
    trials=setdiff(pos_trials,neg_trials);
    data=cbSparse(1:4500,trials+360,unit);
    
    subplot(3,2,1)
    plot(mean(data,2))
    hold on
    plot(white_flash(:,unit),'r')
    legend('Sp','FF')
    title([names{unit},',  white'],'interpreter','none')
    
    axis tight
    a=get(gca,'YLim');
    
    
    trials=setdiff(neg_trials,pos_trials);
    data=cbSparse(1:4500,trials+360,unit);
    subplot(3,2,3)
    plot(mean(data,2))
    hold on
    plot(black_flash(:,unit),'r')
    legend('Sp','FF')
    title([names{unit},',  black'],'interpreter','none')
    
    axis tight
    a1=get(gca,'YLim');
    
    
    trials=1:360;
    trials=setdiff(trials,[pos_trials neg_trials]);
    data=cbSparse(1:4500,trials+360,unit);
    subplot(3,2,5)
    plot(mean(data,2))
    title([names{unit},',  rest'],'interpreter','none')
    
    axis tight
    a2=get(gca,'YLim');
    
    a(1)=0;
    a(2)=max([a a1 a2]);
    for sbpl=1:2:5
        subplot(3,2,sbpl)
        line([500,500],a,'color','k')
        line([2500,2500],a,'color','k')
        axis([0 4500 0 a(2)])
    end

    c=c/max(c(:));
    d=d/max(d(:));
    comb=zeros(800,800,3);
    comb(:,:,1)=bin_full;
    
    subplot(3,2,2)
    imagesc(comb)
    title('original RF')
    
    comb(:,:,2:3)=c;

    subplot(3,2,4)
    imagesc(comb)
    title(['positive stimulus map, ',int2str(length(unique(pos_trials))),' trials'])
    
    comb(:,:,2:3)=d;
    subplot(3,2,6)
    imagesc(comb)
    title(['negative stimulus map, ',int2str(length(unique(neg_trials))),' trials'])
    
    saveas(gcf,['/Users/alexth/Desktop/old_stuff/ONresp/20140121b/cbSparse_plot/',names{unit},'.png'])
    
end




%% find all trials for given location

load('/Users/alexth/Desktop/old_stuff/stim_new/matr_assigned_35_4_700.mat','matr')
matr(:,:,1:4)=[];
matr_tmp=zeros(20,20,357);
matr_tmp(:,:,1:34)=matr(:,:,1:34);
matr_tmp(:,:,35)=zeros(20,20);
matr_tmp(:,:,36:end)=matr(:,:,35:end);
matr=matr_tmp;
clear matr_tmp

coord_pos=cell(20,20);
coord_neg=cell(20,20);
for i=1:20
    for j=1:20
        coord_pos{i,j}=find(matr(i,j,:)==1);
        coord_neg{i,j}=find(matr(i,j,:)==-1);
    end
end
