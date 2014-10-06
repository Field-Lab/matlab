% load pieces_analysis.mat before running this script



% identify number of time points
tp = 1; for ss = 1:length(combos); tp = max([tp combos{ss}]); end

% make rigure
figure(1);clf;

% plot timecourses
ha=make_loop_slider_list(1,1,length(spatial_sensitivity));

% display the spatial sensitivity image for each subset of the chunks
while 1
    % get the slider position
    ss = round(get(ha,'Value'));
    cla;
    
    % plot spatial sensitivity
    subplot('Position',[.05 .1 .9 .8]);
    imagesc(spatial_sensitivity{ss});
    axis image
    set(gca,'xtick',[],'ytick',[])
    
    % plot which subset of the chunks this is
    srplot = subplot('Position',[.05 .9 .9 .05]);
    sr = zeros(1,tp);
    sr(combos{ss}) =1;
    imagesc(sr)
    set(gca,'xtick',[-0.5:tp+0.5],'ytick',[])
    
    uiwait
end

% figure;tt=[];for ss=1:length(spatial_sensitivity);tt=[tt sum(sum(spatial_sensitivity{ss}))];end;plot(tt,'.-')
