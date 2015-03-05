load('/mnt/muench_data/user/alexandra/scripts/single_cone/s')
defaultStream = RandStream.getDefaultStream;
defaultStream.State = s(1).State;
imageArray=zeros(800,800);
b=uint8(randi([0 1],1,160000)*60); % colors 0 and 60
b=repmat(b,2,1);
b=reshape(b,320000,1);
b=reshape(b,800,400);
imageArray(:,1:2:end)=b;
imageArray(:,2:2:end)=b;
linux_Array=imageArray;
load('/mnt/muench_data/user/alexandra/scripts/single_cone/squared')
sum(sum(linux_Array-imageArray))

imageArray=rot90(imageArray,2);

load('/mnt/muench_data/user/alexandra/MEA_data/20121024/coarse_coord/A20121024_CH11_sort1_100_unit_0009_CBFlicker_coars_coord')
coarse_coord=i;
for bigCh=1:1600
    row(bigCh)=floor((bigCh-1)/40)+1;
    col(bigCh)=mod(bigCh-1,40)+1;
end

b=zeros(400,400);
for k=coarse_coord  
    b((col(k)-1)*10+1:(col(k)-1)*10+10,(row(k)-1)*10+1:(row(k)-1)*10+10)=1;
end
fine_coord=b;

a=zeros(40,40);
a(coarse_coord)=1;
figure
subplot(2,1,1)
imagesc(a)
% b=zeros(400,400);
% b(fine_coord)=1;
subplot(2,1,2)
imagesc(b)




path2units='/mnt/muench_data/user/alexandra/MEA_data/20121024/easy_formatted_units/';
units=dir([path2units,'*spike_info_SingleCone.mat']);
load([path2units,units(1).name]);
trial=26;
flips=spike_info.flip_times{trial,1}(:,1);
spikes=cell(length(units),1);
for units_cnt=1:length(units)
    load([path2units,units(units_cnt).name]);
    spikes{units_cnt}=spike_info.spike_times{trial,1};
    spikes{units_cnt}(spikes{units_cnt}<spike_info.flip_times{trial,1}(1,1)+501)=[];
    spikes{units_cnt}(spikes{units_cnt}>flips(end-1))=[];
end

% rngs=spike_info.flip_times{1,1}(:,2);
% rng_order=zeros(flips(end)+1,2);
% rng_order([1; flips(1:end-1)],1)=rngs;
% for i=1:200
%     subscr=flips(1:end-1)+i;
%     subscr(subscr>=flips(end)+1)=[];    
%     subscr(rng_order(subscr,1)~=0)=[];
%     rng_order(subscr,1)=rng_order(subscr-1,1);
% end
% rng_order(spikes,2)=1;

LF=zeros(500,2500);

filt_add=800;

tmp=zeros(filt_add,160000);
load('/mnt/muench_data/user/alexandra/scripts/single_cone/s')
defaultStream = RandStream.getDefaultStream;
defaultStream.State = s(1).State;
% first sequence
cnt=1;Last=1;
tic
while Last<filt_add
    ll=flips(cnt+1)-flips(cnt);
    k=randi([0 1],1,160000);
    k=repmat(k,ll,1);
    tmp(Last:Last+ll-1,:)=k; % considers pixels of size "1", contrast 0/1
    Last=Last+ll;
    cnt=cnt+1;
end
toc
clear k
% fill filter for the first sequence
tic
% for units_cnt=1:length(units)
units_cnt=2
    k=spikes{units_cnt}(spikes{units_cnt}<flips(1)+filt_add)-flips(1)+1;
    length(k)
    for i=k
         LF=LF+tmp(i-499:i,topLeft);
%         LF(:,:,units_cnt)=LF(:,:,units_cnt)+tmp(i-499:i,:);
    end
% end
toc
% fill sequence with shift
Last=800;
tic
while cnt<length(flips)-1
%     tic
    ll=flips(cnt+1)-flips(cnt);
    Last=Last+ll;
    k=randi([0 1],1,160000);
    k=repmat(k,ll,1);
    tmp=[tmp(ll+1:end,:); k];
%     for units_cnt=1:length(units)
        k=Last-(spikes{units_cnt}(spikes{units_cnt}>=flips(cnt)&spikes{units_cnt}<flips(cnt+1))-flips(1)-1);
        for i=k
            LF=LF+tmp(end-i-498:end-i+1,topLeft);
        end
%     end
    cnt=cnt+1
%     toc
end
toc

save('/mnt/muench_data/user/alexandra/MEA_data/20121024/LF_unit0009_s1_nd4','LF')







%% prepared partial sequence

for bigCh=1:1600
    row(bigCh)=floor((bigCh-1)/40)+1;
    col(bigCh)=mod(bigCh-1,40)+1;
end

path2units='/mnt/muench_data/user/alexandra/MEA_data/20121024/CBF_filter/';
units=dir([path2units,'*CBFlicker_LF*.mat']); 
load([path2units,units(1).name]);
fine_coord=zeros(2500,length(units));
findHotSpot=3; % trial number
for units_cnt=1:length(units)
    load([path2units,units(units_cnt).name])    
  
    tmp=CBFfilter(:,:,findHotSpot);
    tmp=tmp';
    [val co]=max(std(tmp));
    centrSize=2;
    centr=[co-centrSize:co+centrSize];
    i=[];
    for j=-centrSize:centrSize
        i=[i centr+40*j];
    end
    i(i>1600)=[];
    i(i<1)=[];
    
    b=zeros(400,400);
    for k=i
        b((col(k)-1)*10+1:(col(k)-1)*10+10,(row(k)-1)*10+1:(row(k)-1)*10+10)=1;
    end
    fine_coord(:,units_cnt)=find(b==1);
end

clear CBFfilter centr centrSize SpikeCount col b row val

% load('/mnt/muench_data/user/alexandra/MEA_data/20121024/coarse_coord/A20121024_CH11_sort1_100_unit_0009_CBFlicker_coars_coord')



load('/mnt/muench_data/user/alexandra/scripts/single_cone/s')
path2units='/mnt/muench_data/user/alexandra/MEA_data/20121024/easy_formatted_units/';
units=dir([path2units,'*spike_info_SingleCone.mat']);
for units_cnt=1:length(units)
    load([path2units,units(units_cnt).name]);
    filter_length=500;
    coneFilter=zeros(2500,filter_length);
    
    part_seq=zeros(300000,2500);
    
    ss=1;
    for trial_counter=14:61
        tic
        trial_counter
        flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
        pxls=spike_info.flip_times{trial_counter}(:,2); % pxs in ms, normalized -1:1
        spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
        defaultStream = RandStream.getDefaultStream;
        defaultStream.State = s(ss).State;
        % first sequence
        cnt=1;Last=1;
        totFlips=length(flips);
        while cnt+1<=totFlips
            ll=flips(cnt+1)-flips(cnt);
            k=randi([0 1],1,160000);
            m=k(fine_coord(:,units_cnt));
            m=repmat(m,ll,1);
            part_seq(Last:Last+ll-1,:)=m; % considers pixels of size "1", contrast 0/1
            Last=Last+ll;
            cnt=cnt+1;
        end
        clear k
        
        
        spikes(spikes<filter_length|spikes>300000)=[];
        spikes(spikes<flips(3)+filter_length)=[];
        for i=1:filter_length
            coneFilter(:,i)=coneFilter(:,i)+sum(part_seq(spikes-i,:))'; % fill actual data
        end
        toc
        ss=ss+1;
        if ss==13
            ss=1;
        end
    end
    save(['/mnt/muench_data/user/alexandra/MEA_data/20121024/cones/',units(units_cnt).name(1:end-28),'_cones'],'coneFilter')
end



figure
colormap('gray')
imagesc(reshape(std(CBFfilter,0,2),50,50))





m=figure(3);

timePoint2show=[30:30:500];
units=dir('/mnt/muench_data/user/alexandra/MEA_data/20121024/cones/*.mat')

for i=1:7
    load(['/mnt/muench_data/user/alexandra/MEA_data/20121024/cones/',units(i).name])
tmp=coneFilter;
tmp=tmp/max(abs(tmp(:)));
b=mean(tmp);
acc=[];
for j=1:length(timePoint2show)
    a=reshape(tmp(:,timePoint2show(j)),50,50);
    acc=[acc a ones(50,1)];
end
acc(:,end-2:end)=[];
subplot(7,1,i)
colormap('gray')
imagesc(acc)
title(units(i).name)
end
