%%
datarun = load_data('/Volumes/Analysis/2012-09-24-5/d03-06-07-norefit/data003/data003');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones_ath(datarun,'d03');

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'date', '2012-09-24-5');

%%
datarun = load_data('/Volumes/Analysis/2012-09-24-5/d00_06-norefit/data003-from-d00_06/data003-from-d00_06');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_cones_ath(datarun,'25.00-all_');
datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'date', '2012-09-24-5');

%%
datarun = load_data('/Volumes/Analysis/2012-09-13-2/d01_09-norefit/data001-from-d01_09/data001-from-d01_09');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);
datarun = load_cones_ath(datarun,'15.00-all_');
datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'date', '2012-09-13-2');

%%

load('/Volumes/Analysis/2012-09-24-5/subunits/d03-06-07-norefit/data003/conepreprocess.mat');


load('/Volumes/Analysis/2012-09-24-5/subunits/data003-from-d00_06/conepreprocess.mat');

%% cone finding locally


if ~isfield(datarun, 'stimulus') || ~isfield(datarun.stimulus, 'java_movie')
    datarun = load_java_movie(datarun);
end
datarun.stimulus.refresh_time = datarun.stimulus.refresh_period/1000;
datarun.stimulus.stimobject = datarun.stimulus.java_movie;
datarun.stimulus.get_movie_frame = @get_movie_frame_from_java_movie;
datarun.stimulus.wc_preprocess = @wc_preprocess_javamovie_rawmovie;

if ~isfield(datarun, 'duration')
    datarun = load_neurons(datarun);
end
datarun.stimulus.duration = min(max(cell2mat(datarun.spikes)), datarun.duration);

start_stim = floor(1+datarun.stimulus.start_time/datarun.stimulus.refresh_time);
end_stim = floor(1+opts.end_time/datarun.stimulus.refresh_time);
stims = start_stim:end_stim;

datarun.cone_inputs = calc_cone_inputs(Wc, stims, datarun.stimulus.stimobject, datarun.stimulus.get_movie_frame);


%%
cellID = 48%7127%4098
ncones = 1779;

datarunID=find(datarun.cell_ids==cellID);
raw_cone_weights=datarun.cones.weights(:,datarunID);
norm_cone_weights=raw_cone_weights/max(raw_cone_weights);
datarun = load_sta(datarun,'load_sta',cellID);
sta = squeeze(datarun.stas.stas{datarunID});
cone_inputs = datarun.cone_inputs;
cone_locs=datarun.cones.centers;
spike_rate=double(datarun.spike_rate);
cell_spike_rate=spike_rate(datarunID,:);

figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(cone_locs(norm_cone_weights>0,:),norm_cone_weights(norm_cone_weights>0),'EdgeColor','r')
viscircles(cone_locs(norm_cone_weights<0,:),-norm_cone_weights(norm_cone_weights<0),'EdgeColor','b')


% plot alt
norm_spike_rate=cell_spike_rate-mean(cell_spike_rate);
[~,sorted_cone_order]=sort(norm_cone_weights,'descend');

alt=zeros(5,ncones);
for j=1:ncones
    cnt=1;
    for i=-0.3:0.2:0.5
        tmp=find(cone_inputs(:,sorted_cone_order(j))<i & cone_inputs(:,sorted_cone_order(j))>i-0.2);
        a=randperm(length(tmp));
        alt(cnt,j)=mean(norm_spike_rate(tmp(a(1:1300))));
%         alt(cnt,j)=mean(norm_spike_rate(tmp));
%         count(cnt,j)=length(tmp);
        cnt=cnt+1;
    end
end

figure
plot(alt(:,1:10))
axis([1 5 -0.2 0.3])

figure
plot(alt(:,500:510))


figure
plot(alt(:,sorted_cone_order(1:10)))

figure
plot(count)

figure
j=1
plot(cone_inputs(:,sorted_cone_order(j)),norm_spike_rate,'*')



bin=0.1;
alt=zeros(length(-0.5:bin:(0.5-bin)),ncones);
sterr=alt;
count=alt;
cnt_cone=1;
binc=zeros(length(-0.5:bin:(0.5-bin)),1);
a=squareform(pdist(cone_locs));
for j=1:ncones
    cnt=1;
%     my_locs=cone_locs(j,:);
%     nncones=find(a(j,:)<10);
%     nncones(nncones==j)=[];
%     t=sum(cone_inputs(:,nncones),2);
%     good_inputs=t<std(t)&t>-std(t);
%     good_inputs=sum(abs(cone_inputs(:,nncones))<0.3,2)==length(nncones);
%     
    
    for i=-0.5:bin:(0.5-bin)
        myData=norm_spike_rate(cone_inputs(:,j)>i & cone_inputs(:,j)<=i+bin);% & good_inputs);
        alt(cnt,cnt_cone)=mean(myData);
        sterr(cnt,cnt_cone)=std(myData)/sqrt(length(myData));
        count(cnt,cnt_cone)=length(myData);
        binc(cnt)=i+bin*0.5;
        cnt=cnt+1;
    end
    max_diff(cnt_cone)=mean(alt(1:2,cnt_cone))-mean(alt(end-1:end,cnt_cone));
    cnt_cone=cnt_cone+1;
end


alt005=alt;

figure
plot(binc,alt(:,sorted_cone_order(1:10)))
b=sum(alt005(1:2,:));


figure
plot(alt(1:3,sorted_cone_order(1500:1600)))


figure
plot(alt005(1:6,sorted_cone_order(1:10)))

figure
clear tmp
excss=zeros(10,ncones);
for i=1:10
    subplot(3,4,i)
    hist(alt(i,:),50)
    tmp(i)=robust_std(alt(i,:));
    title(num2str(tmp(i)))
    if i<6
        excss(i,:)=alt(i,:)>tmp(i)*2;
    elseif i>5
        excss(i,:)=alt(i,:)<tmp(i);
    end
end

sum(excss')
tt=find(sum(excss(1:2,:))==2 & sum(excss(19:20,:))==2)
sorted_cone_order(1:7)
figure
plot(alt005(:,tt))

figure
for i=1:20
    subplot(5,4,i)
    hold on
    plot(alt(:,sorted_cone_order((i-1)*5+1:i*5)));
    plot(tmp,'k', 'linewidth',4)    
end
bl=mean(alt(:,sorted_cone_order(1:5))')
figure
plot(bl)
for i=1:ncones
    myc(i)=corr(alt(:,i),bl');
end

figure
plot(myc)
tt=find(myc>0.95)
sorted_cone_order(1:10)

figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(cone_locs(tt,:),max_diff(tt),'EdgeColor','m')



figure
plot(binc,alt(:,sorted_cone_order(1500:1600)))

figure
plot(max_diff(sorted_cone_order))
figure
plot(alt(:,1328))

figure 
plot(binc,sterr,'r')

figure
plot(count)

figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(cone_locs(max_diff>0.3,:),max_diff(max_diff>0.3),'EdgeColor','m')
viscircles(cone_locs(max_diff>0,:),max_diff(max_diff>0),'EdgeColor','r')
viscircles(cone_locs(max_diff<0,:),-max_diff(max_diff<0),'EdgeColor','b')
viscircles(cone_locs(norm_cone_weights>0,:),norm_cone_weights(norm_cone_weights>0),'EdgeColor','g')
viscircles(cone_locs(norm_cone_weights<0,:),-norm_cone_weights(norm_cone_weights<0),'EdgeColor','y')


alt=zeros(10,ncones);
for j=1:ncones
    cnt=1;
    for i=-0.4:0.1:0.5
        tmp=find(cone_inputs(:,j)<i & cone_inputs(:,j)>i-0.1);        
        a=randperm(length(tmp));
        alt(cnt,j)=mean(norm_spike_rate(tmp(a(1:250))));
        count(cnt,j)=length(tmp);
        stdev(cnt,j)=std(norm_spike_rate(tmp(a(1:250))));
        cnt=cnt+1;
    end
end
c=find(alt(1,:)>0.15);
c=find(alt(1,:)<-0.1);
norm_cone_weights(c)
figure
plot(alt(:,sorted_cone_order(end-4:end)))

figure
plot(count(:,count(5,:)<3000))

figure
plot(count)


figure
plot(stdev)

figure
cnt=1;
for i=c
    subplot(3,3,cnt)
    hist(cone_inputs(:,i),10)
    cnt=cnt+1;
end

figure
b=alt(:,sorted_cone_order(1:100));
plot(b(:,myrfweights<0.5))

% fit weights
mycrs=alt(1:5,sorted_cone_order(1:100))';
mycrsx=repmat([-0.45:bin:0],size(mycrs,1),1);
mycrs(:,5)=max(mean(mycrs(:,5)), 0);

[p resnorm residual] = normcdfxscalesimple(mycrs, mycrsx, 'plot', false, 'title', false);
myrfweights = p(1:end-2)';
myrfweights=myrfweights/max(myrfweights);

figure
plot(norm_cone_weights(sorted_cone_order(1:size(mycrs,1))),myrfweights,'*');
hold on
line([-0.2 1], [-0.2 1], 'color','r')
line([-0.2 1], [0 0], 'color','k')
line([0 0], [-0.2 1], 'color','k')
xlabel('raw cone weights')
ylabel('fitted cone weights')

figure
plot(sort(myrfweights, 'descend'),'b')
hold on
plot(norm_cone_weights(sorted_cone_order(1:size(mycrs,1))),'r')
legend('fitted', 'raw')



locs=datarun.cones.centers(sorted_cone_order(1:size(mycrs,1)),:);
figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(locs,myrfweights,'EdgeColor','b')
hold on
viscircles(locs,norm_cone_weights(sorted_cone_order(1:size(mycrs,1))),'EdgeColor','r')
axis ij



%%
figure
mycones=sorted_cone_order(1:10);
plot(cell_spike_rate)
weighted_inputs=cone_inputs(:,mycones).*repmat(raw_cone_weights(mycones)', size(cone_inputs,1),1);
figure
hist(weighted_inputs(:), 1000)
weighted_inputs(weighted_inputs>0)=0;
summed_input=sum(weighted_inputs');
corr(summed_input',cell_spike_rate')


%%


stimmap='/Volumes/Data/2012-09-24-5/Visual/1234d03/map-0000.txt';
map=load(stimmap);
myboundaries=cell(14,4);
for i=1:4
    BW=map;
    BW(BW~=i)=0;
    boundaries = bwboundaries(BW);
    myboundaries(1:size(boundaries,1),i)=boundaries;
end
figure
imagesc(map)
col='rgby';
hold on
for i=1:4
    for k=1:size(myboundaries,1)
        b = myboundaries{k,i};
        if ~isempty(b)
            plot(b(:,2),b(:,1),col(i));
        end
    end
end



figure
colormap gray
imagesc(sta(:,:,6))
hold on
c=find(a>30);
plot(datarun.cones.centers(c,1),datarun.cones.centers(c,2),'mx')
col='rgby';
for i=1:4
    for k=1:size(myboundaries,1)
        b = myboundaries{k,i};
        if ~isempty(b)
            plot(b(:,2)/2,b(:,1)/2,col(i));
        end
    end
end

c2=find(datarun.cones.centers(:,1)>141.6 & datarun.cones.centers(:,1)<142 &...
    datarun.cones.centers(:,2)==98);

c3=find(datarun.cones.centers(:,2)>92 & datarun.cones.centers(:,2)<92.5 &...
    datarun.cones.centers(:,1)==146);

c4=find(datarun.cones.centers(:,1)>143 & datarun.cones.centers(:,1)<143.5 &...
    datarun.cones.centers(:,2)>89 & datarun.cones.centers(:,2)<89.5);

figure
% plot(my_ci)
plot(my_ci(:,1),cell_spike_rate,'x')
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,1)<i & my_ci(:,1)>i-0.1);
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end
plot(alt,'-x')
hold on
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,2)<i & my_ci(:,2)>i-0.1);
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end
plot(alt,'-xr')
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,3)<i & my_ci(:,3)>i-0.1);
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end
plot(alt,'-xk')

clear alt
cnti=1;cntj=1;
for i=-0.4:0.1:0.5
    cntj=1;
    for j=-0.4:0.1:0
        tmp=find(my_ci(:,2)<i & my_ci(:,2)>i-0.1 & my_ci(:,3)<j & my_ci(:,3)>j-0.1 );
        alt(cnti,cntj)=mean(cell_spike_rate(tmp));
        cntj=cntj+1;
    end
    cnti=cnti+1;
end

figure
plot(alt)
legend('1','2','3','4','5')



clear alt
cnti=1;cntj=1;
for i=-0.4:0.1:0.5
    cntj=1;
    for j=0:0.1:0.5
        tmp=find(my_ci(:,2)<i & my_ci(:,2)>i-0.1 & my_ci(:,1)<j & my_ci(:,1)>j-0.1 );
        alt(cnti,cntj)=mean(cell_spike_rate(tmp));
        cntj=cntj+1;
    end
    cnti=cnti+1;
end

figure
plot(alt)
legend('1','2','3','4','5')


clear alt
cnti=1;cntj=1;
for i=-0.4:0.1:0.5
    cntj=1;
    for j=-0.4:0.1:0.5
        tmp=find(my_ci(:,2)<i & my_ci(:,2)>i-0.1 & my_ci(:,1)<j & my_ci(:,1)>j-0.1 );
        alt(cnti,cntj)=mean(cell_spike_rate(tmp));
        cntj=cntj+1;
    end
    cnti=cnti+1;
end

figure
plot(alt)


clear alt
cnti=1;cntj=1;
for i=-0.4:0.1:0.5
    cntj=1;
    for j=-0.4:0.1:0.5
        tmp=find(my_ci(:,3)<i & my_ci(:,3)>i-0.1 & my_ci(:,1)<j & my_ci(:,1)>j-0.1 );
        alt(cnti,cntj)=mean(cell_spike_rate(tmp));
        cntj=cntj+1;
    end
    cnti=cnti+1;
end

figure
plot(alt)





%% c3 and c4; c3 and c2
clear alt
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,2)<i & my_ci(:,2)>i-0.1 & my_ci(:,3)>-i & my_ci(:,3)<-i+0.1 );
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end

figure
hold on
plot(alt)
clear alt
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,2)<i & my_ci(:,2)>i-0.1 & my_ci(:,1)>-i & my_ci(:,1)<-i+0.1 );
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end
plot(alt,'r')
legend('c3 and c4','c3 and c2')


%% c4 and c3; c4 and c2
clear alt
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,3)<i & my_ci(:,3)>i-0.1 & my_ci(:,2)>-i & my_ci(:,2)<-i+0.1 );
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end

figure
hold on
plot(alt)
clear alt
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,3)<i & my_ci(:,3)>i-0.1 & my_ci(:,1)>-i & my_ci(:,1)<-i+0.1 );
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end
plot(alt,'r')
legend('c4 and c3','c4 and c2')

%% c2 and c3; c2 and c4
clear alt
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,1)<i & my_ci(:,1)>i-0.1 & my_ci(:,2)>-i & my_ci(:,2)<-i+0.1 );
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end

figure
hold on
plot(alt)
clear alt
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(my_ci(:,1)<i & my_ci(:,1)>i-0.1 & my_ci(:,3)>-i & my_ci(:,3)<-i+0.1 );
    alt(cnt)=mean(cell_spike_rate(tmp));
    cnt=cnt+1;
end
plot(alt,'r')
legend('c2 and c3','c2 and c4')

