clf
for i_celltype = 1:2
    x = aggregated_scores.celltype{i_celltype}.scores_WN;
    y = aggregated_scores.celltype{i_celltype}.scores_NSEM;
    
    x(find(x<0)) = 0;
    y(find(y<0)) = 0;
    
    MS = 20;
    subplot(1,2,i_celltype); hold on
    plot(linspace(0,1,100),linspace(0,1,100),'k');
    plot(x,y, 'g.', 'markersize',MS);
    xlim([0 1]);
    ylim([0 1]);
    
end