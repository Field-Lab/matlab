% assume the image_catalog.m for the piece has already been loaded, and a datarun specified



% LOAD NEURONS (for cell ids) AND EI (for array boundaries, electrode positions)
datarun = load_neurons(datarun,'sort_cell_ids',true);
datarun = load_ei(datarun,datarun.cell_ids(1));
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);


% LOAD AXON PATHS

if ~isempty(axon_spreadsheet_path)

    % read in spreadsheet
    [junk, junk, axons_xls] = xlsread(axon_spreadsheet_path);

    % load coordinates of the path for each axon
    num_axons = size(axons_xls,2)/2;
    axons = cell(num_axons,1);
    for aa=1:num_axons
        start_col = aa*2-1;
        % determine if contains numeric data
        if isnumeric(axons_xls{2,start_col})
            % load entries
            axons{aa}(:,1) = cell2mat({axons_xls{2:end,start_col}});
            axons{aa}(:,2) = cell2mat({axons_xls{2:end,start_col+1}});
            % remove the trailing NaN entries, which correspond to empty fields
            axons{aa} = axons{aa}(~isnan(axons{aa}(:,1)),:);
        end
    end

    % transform to array space
    axons_orig = axons;
    % make transformation, based on which image was used to draw the axons
    switch axon_source
        case 'F1'
            TA1F1_reverse = cp2tform(bp_f1,ip_f1,'lwm'); % these should be from the TA1F1 transformation
            T=maketform('composite',TA1F1_reverse,fliptform(TA1));
        case 'A1'
            T=fliptform(TA1);
    end
    % apply to each axon
    for aa=1:num_axons
        if ~isempty(axons{aa})
            axons{aa} = tforminv(T,axons_orig{aa});
        end
    end
end


