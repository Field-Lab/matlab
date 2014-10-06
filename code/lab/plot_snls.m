function datarun = plot_snls(datarun, cell_spec, varargin)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  datarun = my_function(datarun, cell_spec, <params>)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false               show output
% figure            []                  figure or axes to plot in. if 0, make new figure.
%
%
%
% 2009-09  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('figure', 0);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% prepare figure
if params.figure==0
    fig_num = figure;
else
    fig_num = params.figure;
end
figure(fig_num);clf

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% parameters of the loop
start_index=1; index_min=1; index_max=length(cell_indices);

% set up figure with scroll bar
slider = make_loop_slider_list(start_index,index_min,index_max);

% for each panel
while 1
    
    % identify which panel this is
    k=round(get(slider,'Value'));
    
    % clear axes
    cla

    % get cell index and id
    cell_index = cell_indices(k);
    cell_id = datarun.cell_ids(cell_index);
    
    % get snl
    snl = datarun.stas.snls{cell_index};

    % note cell type
    ct = find_cell_types(datarun,cell_id);
    if ct ~= 0;
        ct = datarun.cell_types{ct}.name;
    else
        ct = 'unclassified';
    end
    
    % if it doesn't exist, point this out
    if isempty(snl)
        title(sprintf('cell ID %d (%s) -- no SNL computed',cell_id,ct))

    else

        % make title
        fig_title = sprintf('cell ID %d (%s)',cell_id,ct);

        % plot SNL
        plot_snl_(snl.gen_signal,snl.spikes,'foa',gca,'fit',snl.fit_params,'fig_title',fig_title)

    end
    
    uiwait;
end


fprintf('here')

