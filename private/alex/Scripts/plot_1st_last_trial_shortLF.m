cd('S:\user\alexandra\scripts')
clear
date='20121023';
filter_length=500;
nd='87654321';
filterspath=['F:\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);



path2save=['F:\',date,'\FFFlicker_LF_plots_forPaper\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end

close(figure(1))
figure(6)
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(length(nd));
% simple filter plot
for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    name=str2num(fits(unit).name(end-18:end-15));
    if name==32
        break
    end
end
for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    name=fits(unit).name(1:end-15);
    cellType=common_res_fit(1,:);
    load([filterspath,name,'_FFFlicker_linear_filter'])
    HighFilters=HighFilters./abs(repmat(cellType,500,1));
    cnt=1;
    for i=1:12:12*length(nd)
        k=subplot(rows,cols,cnt);
        tmp=HighFilters(:,[i,i+11]);
        plot(tmp)
        line([0 350],[0 0],'color',[1 1 1]*0.5)
        hleg=legend(int2str(HighSpikeCount([i, i+11])'));
        legend('boxoff')
        set(hleg,'Fontsize',8)
        title(['ND',nd(cnt)])
        set(gca,'XGrid','on')
        set(gca,'xtick',50:50:350)
        axis([0 350 -2 2])
        cnt=cnt+1;
    end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([name,'   High CONTRAST 1st and last trial in ND'],'FontSize',15,'FontWeight','bold','Interpreter','none')
    
    saveas(gcf,[path2save,name,'.emf'])

end
