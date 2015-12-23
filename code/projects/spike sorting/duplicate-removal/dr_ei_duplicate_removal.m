function [result] = dr_ei_duplicate_removal(datarun, varargin)
% EI_DUPLICATE_REMOVAL     Remove duplicates using electrophysiological images
%
% NOTE: this function is one of two cell comparison methods that can be
% used by running dr_remove_duplicates_using_dendrogram
%
% usage:  result = ei_duplicate_removal(dataset, <params>)
%
% arguments: datarun - datarun structure to remove duplicates from
%            varargin - struct or list of optional parameters (see below)
%
% outputs:   datarun - updated datarun structure
%
%
% optional params, their default values, and what they specify:
%
% verbose                  false           show output
% fig_or_axes              []              figure or axes to plot in. if 0, make new figure. if empty, don't plot
% electrode_threshold      5               electrodes with signals below electrode_threshold are set to 0
% significant_electrodes   10              there must be at least this number of electrodes with a significant signal
% corr_threshold           0.95            threshold for determining if a pair of cells are duplicates one another
% space_only               true            use only the spatial information from ei
%
%
%
% 11/08 tamachado, based on EI mapping code by MG and GDF
%

% SET UP ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addRequired('datarun',@isstruct);
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('significant_electrodes', 10, @isnumeric);
p.addParamValue('electrode_threshold', 5, @isnumeric);
p.addParamValue('corr_threshold', 0.8, @isnumeric);
p.addParamValue('space_only', true);

% resolve user input and default values
p.parse(datarun, varargin{:});

% get params struct
datarun = p.Results.datarun;
fig_or_axes = p.Results.fig_or_axes;
significant_electrodes = p.Results.significant_electrodes;
electrode_threshold = p.Results.electrode_threshold;
corr_threshold = p.Results.corr_threshold;
space_only = p.Results.space_only;
verbose = p.Results.verbose;

disp(sprintf('Threshold set to %0.2f',corr_threshold))

% set up plot axes
plot_axes = set_up_fig_or_axes(fig_or_axes);

% get cell numbers
indices = get_cell_indices(datarun, datarun.cell_ids);

% use spatial and temporal information
if ~space_only
    if verbose
        fprintf('EI Duplicate Finding (Spatiotemporal) ...');
        start_time = clock;
    end

    nElectrodes = size(datarun.ei.eis{indices(1)},1);
    eiDepth = size(datarun.ei.eis{indices(1)},2);
    
    sumVector = zeros(nElectrodes*eiDepth,1);
    for i=1:length(indices)
        t = reshape(datarun.ei.eis{indices(i)},nElectrodes*eiDepth,1);
        sumVector = sumVector + t;
    end

    [sorted, ind] = sort(sumVector,'descend');

    mat = zeros(length(indices), nElectrodes * eiDepth);
    for i=1:length(indices)
        t = reshape(datarun.ei.eis{indices(i)},nElectrodes*eiDepth,1);
        if isempty(t), error('EI is empty'); end
        mat(i,:) = t;
    end
    
    % get rid of bad cells that don't have eis
    mat = mat';
    mat(:,sum(isnan(mat),1) > 0) = [];

    %correlate
    corr=corrcoef(mat);
end


%use only spatial info   
if space_only 
    if verbose
        fprintf('EI Duplicate Finding (Spatial) ...');
        start_time = clock;
    end
    
    exclude_list = [];
    mat = zeros(size(datarun.ei.eis{indices(1)},1),length(indices));
    
    %figure out which cells are below threshold
    for i=1:length(indices)
        st = max(datarun.ei.eis{indices(i)},[],2);
        st(st < electrode_threshold)=0;
        if length(find(st > electrode_threshold)) < significant_electrodes
            exclude_list = [exclude_list indices(i)]; %#ok<AGROW>
        end
        mat(:,i) = st;
    end
    
    % get rid of bad cells that don't have eis
    mat(:,sum(isnan(mat),1) > 0) = [];

    %correlate
    corr = corrcoef(mat);
end


result = corr;
% 
% %make a plot visualizing each ei in space relative to the others
% if ~isempty(plot_axes)
%     axes(plot_axes)
%     % plot the thing
% end
% 
% killed = 0; edges = 0;
% dg = zeros(length(indices),length(indices)); 
% dg2 = dg;
% for ii = 1:length(indices)
%     
%     if ii > length(indices)-killed, break; end
%     
%     %find all cells above threshold
%     [sorted_corr] = sort(result(ii, result(ii,:) >= corr_threshold));
%     
%     corr_indices = zeros(length(sorted_corr),1);
%     
%     for jj = 1:length(sorted_corr)
%         corr_indices(jj) = find(result(ii,:) == sorted_corr(jj));
%     end
%     
%     corr_indices = wrev(corr_indices);
%     
%     dg2(ii,corr_indices) = result(ii,corr_indices);
%     dg2(ii,find(dg2(ii,:) == 0)) = -1;
%     dg(ii,corr_indices) = 1;
%     for jj = 1:length(corr_indices)
%         edges = edges + 1;
%         edge(edges,1:2) = [ii corr_indices(jj)];
%     end
% end


% display how long it took
if verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

