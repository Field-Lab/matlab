clear
load('S:\data\alexandra\MEA_data\flash_analysis\flash_moment.mat')



c=snr>25;

close all
figure
set(gcf,'position',[ 1 31 1680 946])

cnt=1;
for i=3:8
    subplot(6,4,cnt)
    a=flash(c&onOff<0,:,i)';
    b=repmat(mean(a(1:400,:)),9000,1);
    a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
%     a(:,8)=[];
    plot(a)
        hold on
    plot(mean(a,2),'k','linewidth',2)
    axis tight
    title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
    cnt=cnt+4;
end

cnt=3;
for i=3:8
    subplot(6,4,cnt)
    a=flash(c&onOff>0,:,i)';
    b=repmat(mean(a(1:400,:)),9000,1);
    a=a-b;
%     b=repmat(max(a),9000,1);
%     a=a./abs(b);
    a=[a(4501:end,:); a(1:4500,:)];
    plot(a)
        hold on
    plot(mean(a,2),'k','linewidth',2)
    axis tight
    title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
    cnt=cnt+4;
end


cnt=2;
for i=3:8
    subplot(6,4,cnt)
    a=moments(c&onOff<0,:,i)';
    b=repmat(mean(a(1:400,:)),5000,1);
    a=a-b;
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
    plot(a)
        hold on
    plot(mean(a,2),'k','linewidth',2)
    axis tight
    title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
    cnt=cnt+4;
end

cnt=4;
for i=3:8
    subplot(6,4,cnt)
    a=moments(c&onOff>0,:,i)';
    b=repmat(mean(a(1:400,:)),5000,1);
    a=a-b;    
%     b=repmat(max(a),5000,1);
%     a=a./abs(b);
    a=[a(2501:end,:); a(1:2500,:)];
    plot(a)
        hold on
    plot(mean(a,2),'k','linewidth',2)
    axis tight
    title(['ON cell ND',int2str(10-i),' n=',int2str(sum(c&onOff>0))])
    cnt=cnt+4;
end






cnt=1;
figure
for i=1:9
    subplot(3,3,i)
    a=flash(c&onOff<0,:,i)';
    b=repmat(mean(a(1:400,:)),9000,1);
    a=a-b;
    maxFlashNeg=max(a(2500:2700,:));
    maxFlashPos=max(a(5000:5200,:));
    semilogy(maxFlashPos./maxFlashNeg);
%     difOFF(i,1:15)=(maxFlash-maxMoment);
    
    hold on
    a=flash(c&onOff>0,:,i)';
%     a(:,7)=[];
    b=repmat(mean(a(1:400,:)),9000,1);
    a=a-b;
    maxFlashNeg=max(a(7000:7200,:));
    maxFlashPos=max(a(500:700,:));
    semilogy(maxFlashPos./maxFlashNeg,'r');
%     difON(i,1:15)=(maxFlash-maxMoment);
%     plot(maxFlash-maxMoment,'r');
    
    axis([0 16 0.1 10])
    line([0 16],[1 1],'color','k')
    title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end

errorbar([(0.87:1:9)' (1.12:1:9.5)'],[mean(difOFF,2) mean(difON,2)],[std(difOFF,0,2) std(difON,0,2)],'*')
hold on
bar([mean(difOFF,2) mean(difON,2)])
set(gca,'xtick',1:9,'xticklabel',('987654321')')
title('Flash peak - Moment peak for positive stimulus. Blue - OFF, red - ON cells')






cnt=1;
for i=1:9
    subplot(3,3,i)
    a=moments(c&onOff<0,:,i)';
    maxMomentPos=max(a(500:1000,:));
    maxMomentNeg=max(a(3000:3500,:));
    plot(maxMomentPos./maxMomentNeg);
    
    hold on
    a=moments(c&onOff>0,:,i)';
    maxMomentPos=max(a(3000:3500,:));
    maxMomentNeg=max(a(500:3100,:));
    plot(maxMomentPos./maxMomentNeg,'r');
    
    axis tight
    title(['ND',int2str(10-i),' n=',int2str(sum(c&onOff<0))])
end

errorbar([(0.87:1:9)' (1.12:1:9.5)'],[mean(difOFF,2) mean(difON,2)],[std(difOFF,0,2) std(difON,0,2)],'*')
hold on
bar([mean(difOFF,2) mean(difON,2)])
set(gca,'xtick',1:9,'xticklabel',('987654321')')
title('Flash peak - Moment peak for positive stimulus. Blue - OFF, red - ON cells')







a=moments(c&onOff<0,:,3)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(3315,:)<-1.5)



a=flash(c&onOff>0,:,3)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
find(a(5157,:)<-1.5)



a=flash(c&onOff>0,:,6)';
b=repmat(mean(a(1:400,:)),9000,1);
a=a-b;
b=repmat(max(a),9000,1);
a=a./abs(b);
figure
plot(a)
find(a(5300,:)<-1.5)




a=moments(c&onOff>0,:,5)';
b=repmat(mean(a(1:400,:)),5000,1);
a=a-b;
b=repmat(max(a),5000,1);
a=a./abs(b);
figure
plot(a)
find(a(3151,:)<-0.96)






figure
for i=3:8
    subplot(2,3,i-2)
    a=flash(c&onOff>0,:,i)';
    b=repmat(mean(a(1:400,:)),9000,1);
    a=a-b;
    % b=repmat(max(a),9000,1);
    % a=a./abs(b);
    a1=max(a(500:1000,:));
    a2=max(a(7000:7500,:));
    plot(sort(a1./a2))    
    hold on
    a=flash(c&onOff<0,:,i)';
    b=repmat(mean(a(1:400,:)),9000,1);
    a=a-b;
    % b=repmat(max(a),9000,1);
    % a=a./abs(b);
    a1=max(a(2480:2980,:));
    a2=max(a(4930:5430,:));
    plot(sort(a2./a1),'r')
    title(['ND',int2str(10-i)])
end




figure
plot(a)
