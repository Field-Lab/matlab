% make NSEM figure comparing NSEM STAs and vision WN STAs
clear
cell_ids = [77 123 140 169 170 202 317 349 351 406 496 502 542 631 661 694 861 1036 1292 1475 1580 1640 2252 2284 2299 2358 2390 2508 2672 2746 2751 2779 2823 3337 3428 3590 3679 3706 3725 3741 3812 3992 4023 4039 4141 4236 4237 4504 4550 4643 4730 4801 4953 4987 5026 5074 5089 5120 5162 5163 5182 5211 5212 5299 5599 5600 5747 5792 5823 5852 5869 5881 5912 6122 6216 6256 6272 6336 6376 6383 6588 6707 6768 6980 6992 6994 7174 7370 7520];
frame = 27;
date = '2013-10-10-0';
sta_datarun = 'data000';
filepath = ['/Users/colleen/Desktop/NSEM_blocked_spikes/',date, '/NS STAs/'];

    datarun.names.rrs_sta_path = ['/Volumes/Analysis/',date, '/', sta_datarun, '/', sta_datarun, '.sta'];
    
    %%% ACTUALLY LOAD UP DATA
    opt = struct('verbose',0,'load_sta',1);
    datarun = load_data(datarun,opt);
for i = 1:length(cell_ids)
    cell_id = cell_ids(i);
    load([filepath, num2str(cell_id), '_NS_STA.mat']);
    NS_STA = sta;
    

   
   
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
