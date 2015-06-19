low_cut = .5
high_cut = 1.5
colors = {'r','g','b','c'};
exp_strings = {'expA','expB','expC'};
i_exp = 1
BD   = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection)); 

MS_OUT = 30;
MS_IN = 12;
exp_string = 'expA'
for i_celltype = 1:2

    if i_celltype == 1, celltype = 'ONP';  cellgroup = allcells{i_exp}.ONP;  end
    if i_celltype == 2, celltype = 'OFFP'; cellgroup = allcells{i_exp}.OFFP; end

    ps_WN   = raster_scores.celltype{i_celltype}.scores_WN.ps_filter;
    ps_NSEM = raster_scores.celltype{i_celltype}.scores_NSEM.ps_filter;

    plotstring_surround = sprintf('%s.', 'r');
    if i_celltype == 1, plotstring_center = 'w.'; titlestring = 'ONP'; end
    if i_celltype == 2, plotstring_center = 'k.'; titlestring = 'OFFP';end

    vals_WN = NaN(1, length(ps_WN));
    for i_cell = 1:length(ps_WN)
        vals_WN(i_cell) = geomean(exp(ps_WN{i_cell}));
    end
    vals_NSEM = NaN(1, length(ps_NSEM));
    for i_cell = 1:length(ps_NSEM)
        vals_NSEM(i_cell) = geomean(exp(ps_NSEM{i_cell}));
    end
    subplot(1,2,i_celltype); hold on; set(gca,'fontsize',8)
    plot(linspace(low_cut,high_cut,100),linspace(low_cut,high_cut,100),'k')

    vals_WN(find(vals_WN<=low_cut))      = low_cut;
    vals_NSEM(find(vals_NSEM<=low_cut))  = low_cut;
    vals_WN(find(vals_WN>=high_cut))     = high_cut;
    vals_NSEM(find(vals_NSEM>=high_cut)) = high_cut;

    plot(vals_WN, vals_NSEM, plotstring_surround,  'markersize', MS_OUT);
    plot(vals_WN, vals_NSEM, plotstring_center,'markersize', MS_IN);

    xlim([low_cut high_cut]);
    ylim([low_cut high_cut]);

    set(gca,'xtick', [.5 .75 1 1.5]);
    set(gca,'ytick', [.5 .75 1 1.5]);

    xlabel('White Noise');
    ylabel('Natural Scenes');
    title(sprintf('Mean Gain of Raster Fitted PS Filter: %s %s', exp_string, titlestring))
end