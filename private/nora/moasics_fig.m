pieces = {'2013-08-19-6', '2012-09-27-3'};
cell_Types = {'On Parasol', 'Off Parasol'};

for i_piece = pieces
   info = experimentinfoNSEM(i_piece{1});
   datarun = load_data([i_piece{1} '/' info.dr.mas]);
   datarun = load_params(datarun);
   for i_cell = cell_Types
       figure; plot_rf_fit(datarun, i_cell{1}, 'labels', false)
       axis equal
       xlim([0 64])
       ylim([0 32])
       title([i_piece{1} ' ' i_cell{1}])
       box on
       set(gca, 'XTick', [])
       set(gca, 'YTick', [])
   end    
end