clear;

datarunpath = '/Volumes/Analysis/2014-11-24-3/mv/data009/data009';
datarun = load_data(datarunpath);
datarun = load_neurons(datarun);
stastempfolder = split(datarun.names.rrs_prefix, filesep);
stastempfolder = fullfile(join(stastempfolder(1:(end-1)), filesep), 'stastemp');

ncells = length(datarun.cell_ids);
cell_has_lightresponse = zeros(ncells, 1);
for k = 6:6
    cellid = datarun.cell_ids(k);
    load(fullfile(stastempfolder, sprintf('sta_%s.mat', num2str(cellid))))
    
    % Possibly try to take out some offset in each sta
    sta = remove_offset_from_sta(sta, sta{1});
    
    % Plot
    plot_sta(sta, 0.1, 1)
    title(sprintf('Cell ID: %s, index: %d', num2str(cellid), num2str(k)))
    
%     % Prompt user: was there a light response?
%     choice = questdlg('Was there a light response?', ...
%         'STA Prompt', ...
%         'Yes','No','No');
%     
%     % Handle response
%     switch choice
%         case 'Yes'
%             cell_has_lightresponse(k) = 1;
%         case 'No'
%             cell_has_lightresponse(k) = 0;
%     end
end

fprintf('Fraction of STAs with a light response: %s.\n', nnz(cell_has_lightresponse)/numel(cell_has_lightresponse));