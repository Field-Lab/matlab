% compute cone weights for subsets of the full run
% gauthier 2010-10

% load datarun with full STAs
datarun = load_data('2008-08-27-5/data003/data003/data003');
datarun = load_neurons(datarun);
datarun = load_params(datarun);

% load dataruns for partials
parts = 1:2;
for dd=1:length(parts)
    dataruns{dd} = load_data(sprintf('2008-08-27-5/data003/data003/pieces/part%d/part%d',parts(dd),parts(dd)));
    dataruns{dd} = load_neurons(dataruns{dd});
    dataruns{dd} = load_sta(dataruns{dd},'load_sta',[],'save_sta',true);
end

% choose cells
cell_spec = {1,2,3,4,5};

cell_ids = get_cell_ids(datarun,cell_spec);

figure(3);clf;hold on

% compute RFs for each subset
datarun_new = datarun;
for cc=1:length(cell_ids)
    fprintf('.')
    cell_id = cell_ids(cc);
    cell_index = get_cell_indices(datarun,cell_id);
    
    % only compute combined STA if the cell is in both
    if ~isempty(intersect(get_cell_ids(dataruns{1},cell_id),cell_id)) && ...
        ~isempty(intersect(get_cell_ids(dataruns{2},cell_id),cell_id))
       
        n_spikes_1 = length(dataruns{1}.spikes{get_cell_indices(dataruns{1},cell_id)});
        n_spikes_2 = length(dataruns{2}.spikes{get_cell_indices(dataruns{2},cell_id)});
    
        new_sta = get_sta(dataruns{1},cell_id)*n_spikes_1;
        new_sta = new_sta + get_sta(dataruns{2},cell_id)*n_spikes_1;
        new_sta = new_sta / (n_spikes_1 + n_spikes_2);
        
        datarun_new.stas.stas{cell_index} = new_sta;
        
        figure(3);plot(n_spikes_1,n_spikes_2,'.');
    end
end
fprintf('\n')

% compute RFs
datarun_new = get_sta_summaries(datarun_new,cell_ids,'verbose',1,'keep_stas',0,'keep_rfs',1,'fig_or_axes',1,...
    'marks_params',struct('strength','vector length','filter',fspecial('gauss',15,0.7),'thresh',5));
