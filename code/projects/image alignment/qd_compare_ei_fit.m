function eicorrs = qd_compare_ei_fit(datarun,axons,axon_id,cell_ids, varargin)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  result = my_function(arg1, <params>)
%
% arguments:     arg1 - argument 1
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
%
%
% date author
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('matchopts', {});

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


cell_indices = get_cell_indices(datarun,cell_ids);
num_cells = length(cell_indices);

position = datarun.ei.position;
datarun = load_ei(datarun,cell_ids);
eis = cell(num_cells, 1);
for cc=1:length(cell_indices)
    cell_index = cell_indices(cc);
    eis{cc} = datarun.ei.eis{cell_index};
end

eicorrs = zeros(num_cells,1);
parfor cc=1:num_cells
    match = compute_ei_axon_match(axons{axon_id}, eis{cc}, position, params.matchopts{:});
    eicorrs(cc) = match.corr;
end

if nargout < 1,
    figure();
    hist(eicorrs,40);
    clear eicorrs;
end
