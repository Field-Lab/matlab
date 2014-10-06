function avg_time_course = average_time_course(datarun, cell_spec, varargin)
% average_time_course     compute average time course for a collection of cells
%
% usage:  avg_time_course = average_time_course(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    avg_time_course - a vector containing the average time course
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           true          	show output
% norm_one          true            make the norm of each time course 1
%                                   before averaging
%
%
% 2009-06  gauthier
% 2014-01, added norm_one keyword
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
avg_time_course = [];
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
        if isempty(avg_time_course)
            avg_time_course = tc;
            
        else % otherwise, just add it in
            
            % be sure it is the correct size
            if size(tc) == size(avg_time_course)
                avg_time_course = avg_time_course + tc;
            else
                % if not, note this
                wrong_size = wrong_size + 1;
                num_tcs = num_tcs - 1;
            end
        end

    end
end

% compute average
avg_time_course = avg_time_course / num_tcs;


% display how long it took
if params.verbose
    if wrong_size > 0
        extra_text = sprintf('\n    (%d cells were not included because their time course had a different number of colors)',wrong_size);
    else
        extra_text = '';
    end
    fprintf('Time course found for %d of %d cells%s.\n',num_tcs,length(cell_indices),extra_text)
end

