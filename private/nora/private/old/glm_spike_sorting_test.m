piece = '2012-08-09-3';
cells = [841, 1276, 1426, 1772, 2101];

analysis_path = ['/Volumes/Analysis/' piece];
exp_info = experimentinfoNSEM(piece);

%% original

WN_datarun = load_data([piece '/' exp_info.dr.slvWN '-from-' exp_info.dr.mas]);
WN_datarun = load_neurons(WN_datarun);
prepped_data_OG = interleaved_data_prep(WN_datarun, [3600 1200], 59, 'cell_spec', cells, 'visual_check', 0);

for i = 1:length(cells)
   IDP_plot_raster(prepped_data, i) 
   title(['On Par ' num2str(cells(i)) ' Original Repeats'])
   set(gcf, 'Position', [ 500 500 800 200])
end

%% mvision data
% ref_data = [analysis_path '/' exp_info.dr.slvWN '-from-' exp_info.dr.mas];
% new_data = [analysis_path '/' exp_info.dr.slvWN '-mVision/' exp_info.dr.slvWN];
% ref_data = [analysis_path '/' exp_info.dr.mas '-mVision/' exp_info.dr.mas];
% mvision_map = crossIdentifyNeuronIDs(ref_data, new_data, [842,1276,1231,1878,1908]);

new_cells = [766,1457,1231,1775,1908];
WN_datarun = load_data([piece '/' exp_info.dr.slvWN '-mVision/' exp_info.dr.slvWN '/' exp_info.dr.slvWN]);
WN_datarun = load_neurons(WN_datarun);
prepped_data_new = interleaved_data_prep(WN_datarun, [3600 1200], 59, 'cell_spec', new_cells, 'visual_check', 1);

%% check STAs
tic; [blockedmoviecell,~,~] = loadmoviematfile(piece, 'WN', '8pix_Identity_8pix','fitmovie'); toc
for i = 1:59
   fitmovie{i} = blockedmoviecell{i}.matrix; 
end

for i = 2:length(cells)
    STA = STA_from_blocks(prepped_data_OG.fitspikes(:,i), fitmovie);
    figure(1); imagesc(sum(STA(:,:,24:27), 3)); title(cells(i))
    STA = STA_from_blocks(prepped_data_new.fitspikes(:,i), fitmovie);
    figure(2); imagesc(sum(STA(:,:,24:27), 3)); title('Spectra')
    pause()
    close all
end

%%
for i = 1:length(cells)
   figure;
   IDP_plot_PSTH(prepped_data_OG, i); hold on
   IDP_plot_PSTH(prepped_data_new, i);
   xlim([0 2])
   title(['On Par ' num2str(cells(i))])
   set(gcf, 'Position', [ 500 500 800 200])
end

%% original

NSEM_datarun = load_data([piece '/' exp_info.dr.slvNSEM '-from-' exp_info.dr.mas]);
NSEM_datarun = load_neurons(NSEM_datarun);
prepped_data_OG = interleaved_data_prep(NSEM_datarun, [7200 3600], 59, 'cell_spec', cells, 'visual_check', 1);

% for i = 1:length(cells)
%    IDP_plot_raster(prepped_data, i) 
%    title(['On Par ' num2str(cells(i)) ' Original Repeats'])
%    set(gcf, 'Position', [ 500 500 800 200])
% end

%% mvision data
ref_data = [analysis_path '/' exp_info.dr.slvNSEM '-from-' exp_info.dr.mas];
new_data = [analysis_path '/' exp_info.dr.slvNSEM '-mVision/' exp_info.dr.slvNSEM];
mvision_map = crossIdentifyNeuronIDs(ref_data, new_data, cells);

%%
new_cells = [767,1278,1232,1833,1893];
NSEM_datarun = load_data([piece '/' exp_info.dr.slvNSEM '-mVision/' exp_info.dr.slvNSEM '/' exp_info.dr.slvNSEM]);
NSEM_datarun = load_neurons(NSEM_datarun);
prepped_data_new = interleaved_data_prep(NSEM_datarun, [7200 3600], 59, 'cell_spec', new_cells, 'visual_check', 1);

%% check STAs
tic; [blockedmoviecell,~,~] = loadmoviematfile(piece, 'NSEM', '8pix_Identity_8pix','fitmovie'); toc
for i = 1:59
   fitmovie{i} = blockedmoviecell{i}.matrix; 
end

for i = 2:length(cells)
    STA = STA_from_blocks(prepped_data_OG.fitspikes(:,i), fitmovie);
    figure(1); imagesc(sum(STA(:,:,24:27), 3)); title(cells(i))
    STA = STA_from_blocks(prepped_data_new.fitspikes(:,i), fitmovie);
    figure(2); imagesc(sum(STA(:,:,24:27), 3)); title('Spectra')
    pause()
    close all
end

%%
for i = 1:length(cells)
   figure;
   IDP_plot_PSTH(prepped_data_OG, i); hold on
   IDP_plot_PSTH(prepped_data_new, i);
   xlim([0 2])
   title(['On Par ' num2str(cells(i))])
   set(gcf, 'Position', [ 500 500 800 200])
end