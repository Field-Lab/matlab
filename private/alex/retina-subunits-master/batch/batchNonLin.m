cellTypes = {'off midget'};
clf
set(gcf,'Position',[342 291 1079 515]);
for iType = 1:length(cellTypes)
    switch cellTypes{iType}
        case 'off midget'
            plotID = 1;
        case 'on midget'
            plotID = 2;
        case 'sbc'
            plotID = 3;
    end
    subplot(2,3,plotID);
    ids = getTypeIds(dat.cellTypes,cellTypes{iType});
    r2_SUB = zeros(1,length(ids));
    r2_LN = zeros(1,length(ids));
    for i=1:length(ids)
        try
            fileName = sprintf('%g-subunit-greedy-0.33-fit.mat',ids(i));
            if isfile(fullfile(dat.analPath,'subunit',fileName));
                load(fullfile(dat.analPath,'subunit',fileName));
            else
                fileName = sprintf('%g-subunit-singletons-fit.mat',ids(i));
                load(fullfile(dat.analPath,'subunit',fileName));
            end
            r2_SUB(i) = res.out_SUB.r2;
            r2_LN(i) = res.out_LN.r2;
            if r2_SUB(i) > 0.2;
                plotNonLin(res.fit_SUB.g);
            end
        catch
            fprintf('(batchNonLin) cannot load cell %g\r',i);
        end
        
    end
    axis equal;
    xlim([-2.5 2.5]);
    ylim([-1 2.5]);
    set(gca,'XTick',[-2 -1 0 1 2]);
    set(gca,'YTick',[-1 0 1 2]);
    
    subplot(2,3,3+plotID);
    plot(r2_LN,r2_SUB,'.','MarkerSize',25,'Color','k');
    axis square equal
    box off
    set(gca,'TickDir','out');
    xlim([0 0.4]);
    ylim([0 0.4]);
    hold on
    plot([0 1],[0 1],'Color','k');
end
%%
% plot the data set in the lower left corner
axes('Position',[0.025 0.025 0.2 0.2]);
axis off
h = text(0,0,lower(sprintf('%s',dat.dataSet)));
set(h,'FontSize',14);
set(h,'Color','k');
