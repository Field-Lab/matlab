load('/Volumes/Analysis/2012-09-24-5/subunits/data003-from-d00_06/conepreprocess.mat');
load('/Volumes/Analysis/2012-09-13-2/subunits/data001-from-d01_09/conepreprocess.mat');

cellID = 1892 %(14k,3)%7127%4098 2026(18k,4) 6871(17k,8)
ncones = 1495;

datarunID=find(datarun.cell_ids==cellID);
raw_cone_weights=datarun.cones.weights(:,datarunID);
norm_cone_weights=raw_cone_weights/max(raw_cone_weights);
datarun = load_sta(datarun,'load_sta',cellID);
sta = squeeze(datarun.stas.stas{datarunID});
cone_inputs = datarun.cone_inputs;
cone_locs=datarun.cones.centers;
spike_rate=double(datarun.spike_rate);
cell_spike_rate=spike_rate(datarunID,:)/max(spike_rate(datarunID,:));

figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(cone_locs(norm_cone_weights>0,:),norm_cone_weights(norm_cone_weights>0),'EdgeColor','r')
viscircles(cone_locs(norm_cone_weights<0,:),-norm_cone_weights(norm_cone_weights<0),'EdgeColor','b')
% select cones around the cell
max_cone=find(raw_cone_weights==max(raw_cone_weights));
dists=squareform(pdist(cone_locs));

dist_thresh=7;
weight_thresh=0.4;

% my_cones=find(dists(:,max_cone)<dist_thresh);
my_cones=find(norm_cone_weights>weight_thresh);

my_cone_weights=norm_cone_weights(my_cones);
my_cone_locs=cone_locs(my_cones,:);
figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(my_cone_locs(my_cone_weights>0,:),my_cone_weights(my_cone_weights>0),'EdgeColor','r')
viscircles(my_cone_locs(my_cone_weights<0,:),-my_cone_weights(my_cone_weights<0),'EdgeColor','b')

% singleton model

% bin cone_inputs
norm_cone_inputs=double(cone_inputs(:,my_cones)*2);
for i=-0.8:0.2:1
    norm_cone_inputs(norm_cone_inputs<=i & norm_cone_inputs>i-0.2)=i-0.1;
end

% initial BC nonlinearity

alt=zeros(10,1);
mean_cell_spike_rate=cell_spike_rate-mean(cell_spike_rate);
cnt=1;
for i=-0.4:0.1:0.5
    tmp=find(cone_inputs(:,max_cone)<=i & cone_inputs(:,max_cone)>i-0.1);
    alt(cnt)=mean(mean_cell_spike_rate(tmp));
    count(cnt)=length(tmp);
    stdev(cnt)=std(mean_cell_spike_rate(tmp));
    cnt=cnt+1;
end
figure
plot(alt,'*-')



BC=cell(100,1);
BC_weights=cell(100,1);
RGC=cell(100,1);
my_corr=zeros(1,100);
iter=1;
% max negative nonlinearity
BC{1}=alt;
BC_weights{1}=norm_cone_weights(my_cones)';
my_corr(1)=-1;
cnt_iter=0;
while iter<100
    
    cnt_iter=cnt_iter+1;
    % pass through nonlinearity, assuming the same parameters for all inputs
    BC_tmp=(rand(10,1)-0.5)/10 + BC{iter};
%     hold on
%     plot(BC_tmp,'r')
    BC_output=norm_cone_inputs;
    
%     figure
%     plot(norm_cone_inputs(1:5,:))

    cnt=1;
    for i=unique(norm_cone_inputs)'        
        BC_output(norm_cone_inputs==i) = BC_tmp(cnt);
        cnt=cnt+1;
    end
    
%     figure
%     plot(BC_output(1:5,:))
%     
    
    % apply weights to BC outputs
    
    BC_weights_tmp = BC_weights{iter} + (rand(1,size(BC_weights{iter},2))-0.5)/10;
    BC_weights_tmp = repmat(BC_weights_tmp,size(BC_output,1),1);
    BC_weighted_output=BC_weights_tmp.*BC_output;
    
%     figure
%     plot(BC_weighted_output(1:5,:))
    
    % sum BC output, normalize
    RGC_input = sum(BC_weighted_output,2);
    RGC_input = RGC_input/max(RGC_input);
    
%     figure
%     plot(RGC_input(1:5))
%     hold on
%     plot(cell_spike_rate(1:5),'r')
    
    % bin RGC input
    minRGC = min(RGC_input(:));
    binRGC=(1-minRGC)/9;
    for i=minRGC:binRGC:(1-binRGC)
        RGC_input(RGC_input>=i & RGC_input<(i+binRGC))=i+binRGC/2;
    end
    
    % pass through RGC nonlinearity
    RGC_init = minRGC:binRGC:1;
    RGC_tmp = [(rand(1,9)-0.5)/5 + RGC_init(1:9) 1];
%     RGC_tmp = RGC_tmp(end:-1:1);
    
%     figure
%     plot(RGC_tmp)

    
    cnt=1;
    RGC_output=RGC_input;
    for i=unique(RGC_input')
        RGC_output(RGC_input==i) = RGC_tmp(cnt);
        cnt=cnt+1;
    end
    
    RGC_output=RGC_output+min(RGC_output);
%     figure
%     plot(RGC_output(1:5))

sse = sum((RGC_output-cell_spike_rate').^2);
sst = sum((cell_spike_rate-mean(cell_spike_rate)).^2);
r2 = 1 - sse/sst;
    
%     tmp=corr(RGC_output, cell_spike_rate');
    tmp=r2;
    if my_corr(iter)<tmp        
        iter=iter+1;
        my_corr(iter)=tmp;
        BC{iter} = BC_tmp;
        BC_weights{iter} = BC_weights_tmp(1,:);
        RGC{iter} = RGC_tmp;
        my_corr(iter)=tmp;
        disp(['step ', int2str(iter), ' , corr ', num2str(my_corr(iter)), ', iter ', int2str(cnt_iter)]);
    end
    
end

figure
plot(RGC{iter})
hold on
plot(RGC{2},'r')

figure
plot(BC{iter})
hold on
plot(BC{1},'r')

figure
plot(BC_weights{iter})
hold on
plot(BC_weights{1},'r')


figure
colormap gray
imagesc(sta(:,:,4))
hold on
viscircles(my_cone_locs(my_cone_weights>0,:),my_cone_weights(my_cone_weights>0),'EdgeColor','m')
viscircles(my_cone_locs(my_cone_weights<0,:),-my_cone_weights(my_cone_weights<0),'EdgeColor','c')
viscircles(my_cone_locs(BC_weights{iter}>0,:),BC_weights{iter}(BC_weights{iter}>0),'EdgeColor','r')
viscircles(my_cone_locs(BC_weights{iter}<0,:),-BC_weights{iter}(BC_weights{iter}<0),'EdgeColor','b')



cnt1=1;
for weight_thresh=-0.2:0.04:0.5
    
    my_cones=find(norm_cone_weights>weight_thresh);
    my_cone_weights=norm_cone_weights(my_cones);
    a=double(cone_inputs(:,my_cones));
    a=a.*repmat(my_cone_weights',size(a,1),1);
    a=sum(a,2);
    
    bin=(max(a)-min(a))/10;
    alt_gc=zeros(10,1);
    cnt=1;
    for i=min(a):bin:max(a)-bin
        alt_gc(cnt)=mean(cell_spike_rate(a>=i & a < i+bin));
        cnt=cnt+1;
    end
    

    cnt=1;
    myrate=zeros(size(cell_spike_rate));
    for i=min(a):bin:max(a)-bin
        myrate(a>=i & a < i+bin)=alt_gc(cnt);
        cnt=cnt+1;
    end
%     corr(myrate',cell_spike_rate');
    
    
    sse = sum((myrate-cell_spike_rate).^2);
    sst = sum((cell_spike_rate-mean(cell_spike_rate)).^2);
    r2 = 1 - sse/sst;

    
    subplot(5,5,cnt1)
    plot(alt_gc)
    title([int2str(length(my_cones)), '  r2 ',num2str(r2)])
    cnt1=cnt1+1;
end

weight_thresh=0.3
my_cones=find(norm_cone_weights>weight_thresh);
my_cone_weights=norm_cone_weights(my_cones);
a=double(cone_inputs(:,my_cones));
a=a.*repmat(my_cone_weights',size(a,1),1);
a=sum(a,2);

bin=(max(a)-min(a))/10;
alt_gc=zeros(10,1);
cnt=1;
for i=min(a):bin:max(a)-bin
    alt_gc(cnt)=mean(cell_spike_rate(a>=i & a < i+bin));
    cnt=cnt+1;
end

a=double(cone_inputs(:,my_cones));
a=a.*repmat(my_cone_weights',size(a,1),1);
a=sum(a,2);

bin=(max(a)-min(a))/10;
cnt=1;
myrate=zeros(size(cell_spike_rate));
for i=min(a):bin:max(a)-bin
    myrate(a>=i & a < i+bin)=alt_gc(cnt);
    cnt=cnt+1;
end
corr(myrate',cell_spike_rate');


sse = sum((myrate-cell_spike_rate).^2);
sst = sum((cell_spike_rate-mean(cell_spike_rate)).^2);
r2 = 1 - sse/sst;

% fit LN
my_cone_weights=norm_cone_weights(my_cones);
tt=double(cone_inputs(:,my_cones));
my_r2=zeros(100,1);
my_r2(1)=r2;
fgc=cell(100,1);fw=fgc;
iter=1;
fgc{1}=alt_gc;
fw{1}=my_cone_weights;
cnt_iter=0;
while iter<100
    
    fit_weights=fw{iter}+(rand(length(my_cones),1)-0.5)/10;
    a=tt.*repmat(fit_weights',size(tt,1),1);
    a=sum(a,2);
    gc=fgc{iter}+(rand(10,1)-0.5)/10;
    
    
    bin=(max(a)-min(a))/10;
    cnt=1;
    myrate=zeros(size(cell_spike_rate));
    for i=min(a):bin:max(a)-bin
        myrate(a>=i & a < i+bin)=alt_gc(cnt);
        cnt=cnt+1;
    end
    
    sse = sum((myrate-cell_spike_rate).^2);
    sst = sum((cell_spike_rate-mean(cell_spike_rate)).^2);
    r2 = 1 - sse/sst;
    if my_r2(iter)<r2
        iter=iter+1;
        my_r2(iter)=r2;
        fw{iter} = fit_weights;
        fgc{iter} = gc;
        disp(['step ', int2str(iter), ' , r2 ', num2str(my_r2(iter)), ', iter ', int2str(cnt_iter)]);
    end
    cnt_iter=cnt_iter+1;
end

figure
plot(fw{1})
hold on
plot(fw{iter},'r')

figure
plot(fgc{1})
hold on
plot(fgc{iter},'r')