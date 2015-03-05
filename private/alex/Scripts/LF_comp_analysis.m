clear
load('/mnt/muench_data/user/alexandra/MEA_data/20121023/FFFlicker_LF_raw/A20121023_CH24_sort1_100_unit_0009_FFFlicker_linear_filter_raw.mat')
clear LowFilterRow LowSpikeCount
cnt=1;
nd='87654321234567';
figure
for i=1:12:12*8
    a=HighFilterRow{i};
    a(sum(abs(a),2)<50,:)=[];
    % k=mean(a);
    % figure
    % plot(k)
    
    b=pdist(a(:,50:350));
    b=squareform(b);
    l(cnt)=mean(b(:));
    subplot(2,4,cnt)
    imagesc(b)
    title(['ND',nd(cnt)])
    cnt=cnt+1;
end


cnt=1;
for i=1:12:12*8
    a=HighFilterRow{i};
    a(sum(abs(a),2)<50,:)=[];
    k(cnt,1:500)=mean(a);
    cnt=cnt+1;
end
figure
plot(k','LineWidth',3)
legend(['87654321']')


i=25
a=HighFilterRow{i};
a(sum(abs(a),2)<50,:)=[];

b=pdist(a(:,50:350));
b=squareform(b);
plot(mean(a))


tt=repmat(mean(a),size(a,1),1);
distance2template=sqrt(sum((tt'-a').*(tt'-a')));
plot(distance2template)
[val,ID]=sort(distance2template);
figure
plot(mean(a(ID(1:50),:)))
figure
imagesc(a(ID,:))
figure
plot(mean(a(sum(a(:,120:220)')>10,:)))

for k=1:length(b)
    [val,ID]=sort(b(:,k));
    tt(k,1:length(b))=val;
end
figure
plot(tt')
plot(mean(a(ID(val>0.5),:)))
hold on
m=ID(1);
[val,ID]=sort(b(:,m));
plot(mean(a(ID(val>0.5),:)),'k')

plot(max(a(:,1:350)))
hold on
plot(min(a(:,1:350)))

 figure;
 b=arepmat(mean(a),length(a),1);
 plot(std(b))
 plot(mean(a)-std(a),'r')
 hold on
 plot(mean(a)+std(a),'r')
 plot(mean(a))


 
 
 a=HighFilterRow{25};
 a(sum(abs(a),2)<50,:)=[];
 plot(mean(a))
 hold on
 plot(std(a)*10-3,'r')
 for i=1:25
     subplot(5,5,i)
     plot(a(i+600,:))
     hold on
     plot(mean(a)*2,'r')
 end
 
 
 
 