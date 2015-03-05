clear
load('/mnt/muench_data/user/alexandra/MEA_data/20121017/waveforms/CH1_20121017')

% align by min
[~,pos]=min(spikes(10:14,:));
spikes_tmp=zeros(20,size(spikes,2));
pos=pos+9;
for i=10:14
    first_point=i-9;
    last_point=i+10;
    spikes_tmp(:,pos==i)=spikes(first_point:last_point,pos==i);
end
spikes=spikes_tmp;
clear val pos first_point last_point vals ranges max_val pos_max

figure
plot(min(spikes(:,1:10:end)),1:size(spikes(:,1:10:end),2),'.','MarkerSize',1)
axis tight

%% get template 1

thr=-1050;

val=min(spikes);
spikes_tmp=spikes(:,val<thr);

% normalize amplitude (from 0 to -1)
spikes_tmp=spikes_tmp-repmat(max(spikes_tmp),20,1);
spikes_tmp=-spikes_tmp./repmat(min(spikes_tmp),20,1);


sparseness=ceil(size(spikes_tmp,2)/5000);
tmp=spikes_tmp(:,1:sparseness:end);
% figure
% plot(tmp)

tic
common=pdist(tmp');
common = squareform(common);
[common,IX]=sort(common);
tt=common(1:2500,:);
[~,b]=min(mean(tt));
a=IX(:,b);
template1=mean(tmp(:,a(1:2500)),2);
toc


%% get template 2
thr=-600;
thr2=-950;

val=min(spikes);
spikes_tmp=spikes(:,val<thr);
val=min(spikes_tmp);
spikes_tmp=spikes_tmp(:,val>thr2);

% normalize amplitude (from 0 to -1)
spikes_tmp=spikes_tmp-repmat(max(spikes_tmp),20,1);
spikes_tmp=-spikes_tmp./repmat(min(spikes_tmp),20,1);


sparseness=ceil(size(spikes_tmp,2)/5000);
tmp=spikes_tmp(:,1:sparseness:end);
% figure
% plot(tmp)

tic
common=pdist(tmp');
common = squareform(common);
[common,IX]=sort(common);
tt=common(1:2500,:);
[~,b]=min(mean(tt));
a=IX(:,b);
template2=mean(tmp(:,a(1:2500)),2);
toc

clear sparseness spikes_tmp common IX tt val a b thr thr2 tmp

%% plot templates
figure
plot([template1,template2],'linewidth',3)
title('templates')
legend('template 1','template 2','Location','SouthEast')


%% analysis
figure
i=1;
for thr=[-200 -350 -450 -550 -650 -950]
    subplot(2,3,i)
%     thr=-350; % common threshold
    val=min(spikes);
    spikes_tmp=spikes(:,val<thr);
    % normalize amplitude (from 0 to -1)
    spikes_tmp=spikes_tmp-repmat(max(spikes_tmp),20,1);
    spikes_tmp=-spikes_tmp./repmat(min(spikes_tmp),20,1);
    clear val
%     figure
%     plot([template1,template2,mean(spikes_tmp,2)],'linewidth',3)
%     legend('template 1','template 2','mean of all','Location','SouthEast')
    
%     % distance between templates
%     sqrt(sum((template1-template2).*(template1-template2)))
    
    tt=repmat(template1,1,size(spikes_tmp,2));
    tt=sqrt(sum((tt-spikes_tmp).*(tt-spikes_tmp)));
    
    pp=repmat(template2,1,size(spikes_tmp,2));
    pp=sqrt(sum((pp-spikes_tmp).*(pp-spikes_tmp)));
    
    plot(tt,pp,'.','MarkerSize',1)
    xlabel('template 1')
    ylabel('template 2')
    title(['thr=', int2str(thr)])
    i=i+1;
end






figure
plot(pp,1:length(tt),'.','MarkerSize',1)

[~,a]=sort(tt);
[~,b]=sort(pp);

figure
hold on
plot(spikes_tmp(:,a(1:10000)),'b')
plot(spikes_tmp(:,b(1:10000)),'g')
plot([template1,template2],'k','linewidth',3)

