clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','onOff','bothContr')
clear a b c
for i=1:8
    subplot(2,4,i)
    b=hc_gain(bothContr&onOff>0,i);  
    b1=lc_gain(bothContr&onOff>0,i);
    plot(b,b1,'or')
    hold on
    
    b=hc_gain(bothContr&onOff<0,i);
    b1=lc_gain(bothContr&onOff<0,i);
    plot(b,b1,'o')
    axis([-0. 0.15 -0. 0.07])
end

figure
for i=1:8
    subplot(2,4,i)
    b=hc_gain(bothContr&onOff>0,i);  
    onHC{i}=b(b>0&b<0.15);   
    b1=lc_gain(bothContr&onOff>0,i);
    onLC{i}=b1(b1>0&b1<0.15);
    ons=b./b1;
    ons=ons(ons>0&ons<15);
    onRat{i}=ons;
    plot(ons,'or')
    hold on
    
    b=hc_gain(bothContr&onOff<0,i);
    offHC{i}=b(b>0&b<0.15);
    b1=lc_gain(bothContr&onOff<0,i);
    offLC{i}=b1(b1>0&b1<0.15);
    offs=b./b1;
    offs=offs(offs>0&offs<15);
    offRat{i}=offs;
    plot(offs,'o')
    axis([-0. 150 0. 15])
    line([0 150],[5,5],'color','k')
    line([0 150],[1,1],'color','k')
    p(i) = ranksum(offs,ons);
end

for i=1:7
    pon(i)=ranksum(onHC{i},onHC{i+1});
    pon1(i)=ranksum(onLC{i},onLC{i+1});
    poff(i)=ranksum(offHC{i},offHC{i+1});
    poff1(i)=ranksum(offLC{i},offLC{i+1});
end



clear p_on p_off
for i=1:7
    for j=i+1:8
        p_on(i,j-1) = ranksum(onRat{i},onRat{j});
        p_off(i,j-1) = ranksum(offRat{i},offRat{j});
    end
end


figure
subplot(1,2,1)
a=p_on;
a(a>=0.05)=2;
a(a>0&a<0.05)=1;
imagesc(a)
set(gca,'xtick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
set(gca,'ytick',1:7,'yticklabel',{'8','7','6','5','4','3','2'})











