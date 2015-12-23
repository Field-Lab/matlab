function [time_courses, noTc] = get_time_courses_matrix(datarun, cell_ids, varargin)

% get_time_courses     get time course vectors of several cells and puts
% them into a matrix
%
% usage:  get_time_courses(datarun, cell_ids)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_ids - which cells (see get_cell_ids for options)
%                
% output:   time_courses: a TxC matrix with T time points and C cells
%
% 2012-08  sneha 

p = inputParser;
p.addParamValue('norm_one', false, @islogical);

p.parse(varargin{:});


b = 1;
noTc = [];
for a = 1:length(cell_ids)
        cell_index = get_cell_indices(datarun,cell_ids(a));
        if(isempty(datarun.stas.time_courses{cell_index,1}))
            noTc(1,b) = cell_ids(1,a);
            b = b+1;
        else
        time_courses(:, a) = datarun.stas.time_courses{cell_index,1};
        
        if p.Results.norm_one
            time_courses(:,a) = time_courses(:,a) ./ norm(time_courses(:,a));
        end
        
        end
end
time_courses = time_courses(:,any(time_courses));

