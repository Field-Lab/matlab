function [all_corr,cell_ids] = find_matches(datarun,axons,axon_id,varargin)
% find_matches     overlay EIs on an axon, scroll through different EIs
%
% usage:  find_matches(datarun,axons,axon_id,varargin)
%
% arguments:     datarun - datarun struct
%                  axons - axons cell array, each element an Nx2 matrix of axon path
%                 xon_id - which axon to plot
%               varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% ei_cutoff         -1              see plot_ei_
% ei_scale          1               see plot_ei_
% cell_ids          []             	which cell IDs to plot.  if empty, choose all cell ids from electrodes near the axon start
% radius           	1         		radius of electrodes in which to get cell ids
% sort              true            plot EIs by order of correlation (highest first)
%
%
% 2010-01  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius', 1);
p.addParamValue('ei_cutoff', -1);
p.addParamValue('ei_scale', 1);
p.addParamValue('cell_ids', []);
p.addParamValue('sort', true);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% for a specified axon and electrode, will generate plots of the axon overlaid on all nearby EIs


% PARAMETERS
fig = 1;



% IDENTIFY POTENTIALLY MATCHING EIs, IF NEEDED

if ~isempty(params.cell_ids)
    cell_ids = params.cell_ids;

    % show user what's going on
else

    % get axon start
    axon_start = axons{axon_id}(1,:);

    % identify nearest electrode
    dists = sum([datarun.ei.position(:,1) - axon_start(1)   datarun.ei.position(:,2) - axon_start(2)].^2,2);
    [junk,electrode]=min(dists);

    % identify nearby electrodes
    array_info = load_array_info(datarun,2);
    search_electrodes = array_info.electrodeMap.getAdjacentsTo(electrode, params.radius);

    % get list of all possible cell IDs
    cell_ids = [];
    for ee = search_electrodes'
        cell_ids = [cell_ids; cell_ids_from_electrode(ee)];
    end
    % pare down to actual cell IDs
    cell_ids = intersect(cell_ids,datarun.cell_ids);

    % show user what's going on
    fprintf('axon %d.  closest electrode %d.  also searching ',axon_id,electrode)
    for ee = setdiff(search_electrodes,electrode)
        fprintf('%d ',ee')
    end
    fprintf('\n')
    
end


% compute correlations
eicorr = qd_compare_ei_fit(datarun,axons,axon_id,cell_ids);
corr = eicorr;

if params.sort
    % sort cell_ids
    [~,i] = sort(corr(:), 1, 'descend');
    cell_ids = cell_ids(i);
    corr     = corr(i);
end







% PLOT EACH EI WTIH AXON

% generate predicted EI
predicted_ei = generate_predicted_ei(axons{axon_id},datarun.ei.position);

% prep figure
figure(fig);clf

% create slider control
ha = make_loop_slider_list(1, 1, length(cell_ids), {@slider_plot, datarun, axons, axon_id, cell_ids, corr, params});

% plot once before any clicks
slider_plot(ha, [], datarun, axons, axon_id, cell_ids, corr, params);


function slider_plot(handle, ~, datarun, axons, axon_id, cell_ids, all_corr, params)
% display one frame of the STA

% get the slider position
cc = round(get(handle,'Value'));
cell_id = cell_ids(cc);
rho = all_corr(cc);
cla;

% load the ei
datarun = load_ei(datarun,cell_id);

% compute correlation with predicted ei
%ei = datarun.ei.eis{get_cell_indices(datarun,cell_id)};
%rho = corr(predicted_ei,max(abs(ei),[],2));

% plot spatial sensitivity
qd_plot_aligned_rgc(datarun,axons,cell_id,axon_id,'foa',gca,'ei_scale',params.ei_scale,...
    'ei_cutoff',params.ei_cutoff,'title',sprintf('corr = %0.3f, ',rho));

