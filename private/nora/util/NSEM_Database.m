%% Easy to get in matlab
function NSEM_Database(datarun_class_name)

%% pick cell types
%datarun_class = load_data('2011-12-04-0/data000');
datarun_class = load_data(datarun_class_name);
datarun_class = load_params(datarun_class);
cell_types = [];
for i = 1:length(datarun_class.cell_types)
    cell_type = datarun_class.cell_types{i}.name;
   
    % automatically sort some types
    if any(strcmpi(cell_type, {'On Parasol', 'Off Parasol', 'On Midget', 'Off Midget', 'SBC', 'On Large', 'Off Amacrine'}))
         % automatically include
        cell_types = [cell_types i];
        manual = 0;
    elseif any(strcmpi(cell_type, {'junk', 'crap', 'unclassified'}))
        % automatically don't include
        manual = 0;
    elseif length(datarun_class.cell_types{i}.cell_ids)==0
        manual = 0;
    else
        % make the person choose
        manual = 1;
    end
    
    % person chooses
    if manual
        f = figure;
        plot_rf_fit(datarun_class, cell_type);
        title(cell_type)
        x = waitforbuttonpress;
        if ~x
            cell_types = [cell_types i];
        end
        close(f)
    end
end

%%
for i = cell_types
    figure;
    plot_rf_fit(datarun_class, datarun_class.cell_types{i}.name);
    title(datarun_class.cell_types{i}.name)
end

% Mosaics and overall quality
% Spike Rate
% generator signal
% parasol RF size
% stability

%% from notebook
% Stimuli
% Piece details: culture, RPE, temperature, eccentricity

%% spike sorti
end