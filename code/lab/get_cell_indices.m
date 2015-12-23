function [cell_numbers, cell_type, cell_type_number] = get_cell_indices(datarun, cell_specification)
% get_cell_indices     get cell numbers matching to the cell specification
%
% usage:  [cell_numbers, cell_type, cell_type_number ] = get_cell_indices( datarun, cell_specification )
%
% arguments:      datarun -	datarun struct containing the following fields
%                               datarun.cell_ids - list of cell ids
%                               datarun.cell_types - cell types listed in standard format
%      cell_specification -	one of the following:
%                               vector of cell ids
%                               cell array of cell type names and/or numbers
%                               'all'
%
% outputs:   cell_numbers - row vector of cell numbers
%               cell_type - name of the cell type(s)
%        cell_type_number - ordinal number of the cell type(s)
%
%
% NOTE: if more than one cell type was specified, then the following applies:
%           if there are one or zero output arguments, all cell_ids will be returned as a list
%           if there are more than one output arguments, all outputs will be returned as cell arrays, grouped by cell type
%
%
% examples:
%
%       get_cell_indices( datarun, 1153 )
%           ans = [ <cell number of cell ID 1153> ]
%
%       get_cell_indices( datarun, [1153 3005] )
%           ans = [ <cell numbers of cell IDs 1153 and 3005> ]
%
%       get_cell_indices( datarun, 'ON parasol' )
%           ans = [ <cell numbers of ON parasol cells> ]
%
%       get_cell_indices( datarun, {'ON parasol','OFF parasol'} )
%           ans = [ <cell numbers of ON and OFF parasol cells> ]
%
%       get_cell_indices( datarun, {'ON parasol',2} )
%           ans = [ <cell numbers of ON and OFF parasol cells> ]
%
%       cell_nums = get_cell_indices( datarun, {'ON parasol','OFF parasol'} )
%           cell_nums = [ <cell numbers of ON and OFF parasol cells> ]
%
%       [cell_numbers, cell_type, cell_type_number ] = get_cell_indices( datarun, {'ON parasol','OFF parasol'} )
%           cell_numbers = { [ <cell numbers of ON parasol cells> ], [ <cell numbers of OFF parasol cells> ] }
%           cell_type = {'ON parasol','OFF parasol'}
%           cell_type_number = {1,2}
%
%
% gauthier 2008-03-02
% 


% identify which kind of specification it is
switch class(cell_specification)
    
    case {'int16','uint16','int32','uint32','int64','uint64','single','double'}  % a list of cell numbers
        cell_numbers = nums_from_ids(datarun, cell_specification);
        cell_type = [];
        cell_type_number = [];
    
    case 'char'  % a cell type name
        if strcmpi(cell_specification,'all')
            cell_numbers = 1:length(datarun.cell_ids);
            cell_type = [];
            cell_type_number = [];
        else
            [cell_ids, cell_type, cell_type_number] = get_cell_ids_by_type( datarun, cell_specification);
            cell_numbers = nums_from_ids(datarun, cell_ids);
        end

    case 'cell'  % a list of names or cell type numbers
        % handle one element of the list at a time
        for cc = 1:length(cell_specification)
            [ cell_ids{cc}, cell_type{cc}, cell_type_number{cc} ] = ...
                get_cell_ids_by_type( datarun, cell_specification{cc} );
            % convert
            cell_numbers{cc} = nums_from_ids(datarun, cell_ids{cc});
        end

        % if only one cell type was specified, don't return a cell array
        if length(cell_specification) == 1
            cell_numbers = cell_numbers{1};
            cell_type = cell_type{1};
            cell_type_number = cell_type_number{1};
        else % if more than one type was specified, check number of output arguments
            switch nargout
                case {0,1} %if nothing more than cell numbers will be returned 
                    % return the union of cell numbers
                    cell_numbers = union_cell_numbers(cell_numbers);
            end
        end

    otherwise
        error('Cell specification not recognized.  Type ''help get_cell_indices''.')
end


function [cell_ids, cell_type, cell_type_number] = get_cell_ids_by_type(datarun, cell_spec)
% return information about the type specified in cell_specification
[cell_type_number, cell_type] = get_cell_type_nums(datarun, cell_spec);
cell_ids = datarun.cell_types{cell_type_number}.cell_ids;


function cell_numbers = union_cell_numbers(cell_array)
% take the union of vectors in a cell array and sort them

cell_numbers = cell_array{1};

for cc = 2:length(cell_array)
    cell_numbers = union(cell_numbers,cell_array{cc});
end

cell_numbers = sort(cell_numbers);

