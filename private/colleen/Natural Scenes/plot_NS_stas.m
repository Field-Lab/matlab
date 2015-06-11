% make NSEM figure comparing NSEM STAs and vision WN STAs
cell_ids = [997 2228 3362 5691 7307];
frame = 27;
date = '2013-08-19-6';
sta_datarun = 'data000';
filepath = ['/Users/colleen/Desktop/NSEM_blocked_spikes/',date, '/NS STAs/'];


for i = 1:length(cell_ids)
    cell_id = cell_ids(i);
    load([filepath, num2str(cell_id), '_NS_STA.mat']);
    NS_STA = sta;
    
    datarun.names.rrs_sta_path = ['/Volumes/Analysis/',date, '/', sta_datarun, '/', sta_datarun, '.sta'];
    
    %%% ACTUALLY LOAD UP DATA
    opt = struct('verbose',0,'load_sta',1);
    datarun = load_data(datarun,opt);
   
   
    cell_ind = find(datarun.cell_ids == cell_id);
    if ~isempty(cell_ind)
    sta = datarun.stas.stas{cell_ind};
    sta = norm_image(sta);
    
    fig=figure('PaperSize',[12 8], 'Position',[1, 1, 800, 400]);
    
    
    h = subplot('position', [0.03 0.25, 0.45,0.4]);
    imagesc(NS_STA(:,:,27))
    colormap gray
    axis off
    axis equal
    axis tight
    title({'Natural Scenes', ['Frame ', num2str(frame), '/30']})
    
    g = subplot('position', [0.52 0.25, 0.45,0.4]);
    imagesc(squeeze(sta(:,:,:, 27)))
    axis off
    axis equal
    axis tight
    title({'White Noise', ['Frame ', num2str(frame), '/30']})
    uicontrol('style','text','string',[date, ' Cell ', num2str(cell_id)],'position',[275 300 250 25], 'FontSize', 22, 'FontWeight', 'bold','FontName','Helvetica', 'BackgroundColor', [1 1 1]);
    set(gcf, 'PaperPositionMode', 'auto')
    print(fig,'-dpdf',sprintf('%s%s.pdf',filepath,['Cell_',num2str(cell_id)]));
    else
        disp(['ID ', num2str(cell_id), ' is missing']);
    end
end
