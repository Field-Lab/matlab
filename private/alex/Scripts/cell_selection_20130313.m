figure
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


figure

i=4
a=flash(c&onOff<0,:,i)';
reb=mean(a(1000:1250,:))./mean(a(1:400,:))>1&~isinf(mean(a(1000:1250,:))./mean(a(1:400,:)));
k=a(:,reb);
b=repmat(mean(k(1:400,:)),9000,1);
k=k-b;
% k=k(:,k(700,:)<0);
subplot(2,1,1)
plot(k)
hold on
plot(mean(k,2),'k','linewidth',2)
axis tight
title(['ND',int2str(10-i),' n=',int2str(size(k,2))])

k=a(:,~reb);
subplot(2,1,2)
b=repmat(mean(k(1:400,:)),9000,1);
k=k-b;
% k=k(:,k(700,:)<0);
plot(k)
hold on
plot(mean(k,2),'k','linewidth',2)
axis tight
title(['ND',int2str(10-i),' n=',int2str(size(k,2))])




figure

i=4
a=flash(c&onOff<0,:,i)';
a=a(:,32);
plot(a)

f=find(diff(timing_flash));
figure
cnt=1;
for i=3:8    
    
    for cellN=[5 38 17 23]
        subplot(6,4,cnt)
        a=flash(c&onOff<0,:,i)';
        a=a(:,cellN);
        
        %     b=repmat(max(a),9000,1);
        %     a=a./abs(b);
        %     a(:,8)=[];
        plot(a)
        axis tight
        title(['ND',int2str(10-i),' cell ID =',int2str(cellN)])
        for k=1:4
            line([f(k) f(k)],[0,max(a)],'color','k')
        end
        cnt=cnt+1;
    end
end

figure
cnt=1;
for i=3:8    
    
    for cellN=[3 4 10 11]
        subplot(6,4,cnt)
        a=flash(c&onOff>0,:,i)';
        a=a(:,cellN);
        plot(a)
        axis tight
        title(['ND',int2str(10-i),' cell ID =',int2str(cellN)])
        for k=1:4
            line([f(k) f(k)],[0,max(a)],'color','k')
        end
        cnt=cnt+1;
    end
end

