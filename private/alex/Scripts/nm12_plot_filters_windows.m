clear
load(['S:\data\alexandra\MEA_data\sine_analysis\HL10_accumParameters'])

path2save='S:\data\alexandra\MEA_data\sine_analysis\HL10_pictures\';

figure
set(gcf,'position',[1  31 1680 946]);
nds='87654321';
for i=1:size(LinearFilter,4)
    cnt=1;
    for j=1:2:16
        subplot(2,4,cnt)  

        a=reshape(mean(LinearFilter(:,1:2:end,j:j+1,i),2),500,2);
        plot(a)
        hold on        
        a=reshape(mean(LinearFilter(:,2:2:end,j:j+1,i),2),500,2);
        plot(a)
        title(['ND',nds(cnt)],'FontWeight','bold')
        cnt=cnt+1;
        line([0 500],[0,0],'color','k')
        hold off
    end
    subplot('Position',[0.5 0.96 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    name=names{i}(1:end-4);
    if i>=33&&i<67
        name=[name(1:9),'_1',name(10:end)];
    elseif i>=67
        name=[name(1:9),'_2',name(10:end)];
    end
    title(name,'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,name,'.emf']);
    
end

for i=1:103
    name=names{i}(1:end-4);
    if i>=33&&i<67
        name=[name(1:9),'_1',name(10:end)];
    elseif i>=67
        name=[name(1:9),'_2',name(10:end)];
    end
    disp(name)
end