function SE_TC = standard_error_time_course(datarun, cell_spec, varargin)
% standard_error_time_course     compute standard erro of the time course for a collection of cells
%
% usage:  SE_TC = standard_error_time_course(datarun, cell_spec, varargin)
%
%  NOTE: currently this function only works for BW time courses
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    SE_TC - standard error of the time course
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           true          	show output
% norm_one          true            make the norm of each time course equal
%                                   to one before averaging
%
%
% 2014-01  GDF
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', true, @islogical);
p.addParamValue('norm_one', true, @islogical);
p.addParamValue('peak_norm', false, @islogical);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;





% BODY OF THE FUNCTION

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% initialize
TC_ensemble = [];
num_tcs = 0;
wrong_size = 0;

% loop through cells
for cc = 1:length(cell_indices)

    % get cell index and id
    cell_index = cell_indices(cc);

    % get time course of this cell
    tc = datarun.stas.time_courses{cell_index};
    
    if params.norm_one; tc = tc ./ norm(tc); end
    if params.peak_norm; tc = tc ./ abs(ext(tc)); end


    if ~isempty(tc)
        % incorporate into time courses to be averaged
        num_tcs = num_tcs + 1;

        % if this is the first one, use it to initialize
        if isempty(TC_ensemble)
            TC_ensemble(cc,:) = tc;
            
        else % otherwise, just add it in
            
            % be sure it is the correct size
            if length(tc) == length(TC_ensemble(1,:))
                TC_ensemble(cc,:) = tc;
            else
                % if not, note this
                wrong_size = wrong_size + 1;
                num_tcs = num_tcs - 1;
            end
        end

    end
end

% compute average
SE_TC = std(TC_ensemble, 0, 1) / sqrt(num_tcs-1);


% display how long it took
if params.verbose
    if wrong_size > 0
        extra_text = sprintf('\n    (%d cells were not included because their time course had a different number of colors)',wrong_size);
    else
        extra_text = '';
    end
    fprintf('Time course found for %d of %d cells%s.\n',num_tcs,length(cell_indices),extra_text)
end

