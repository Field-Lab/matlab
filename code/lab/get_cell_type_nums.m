function [cell_type_nums, cell_type_names] = get_cell_type_nums(datarun, cell_specs)
% GET_CELL_TYPE_NUMS    Parse cell spec to get cell type number
% usage: cell_type_nums = get_cell_type_nums(datarun, cell_specs)
% arguments: datarun
%            cell_spec - Can be just the number, in which case we basically
%                        do nothing.  Can be a string, in which
%                        case we search for the matching cell type in
%                        DATARUN and return its index.  Or can be an array
%                        of numbers or a cell array of numbers and strings.
%
% 2010-01 phli, abstracted out of GET_CELL_INDICES, now handles cell array
%

if iscell(cell_specs)
    cell_type_nums = zeros(size(cell_specs));
    cell_type_names = cell(size(cell_specs));
    for i = 1:numel(cell_specs)
        [cell_type_nums(i), cell_type_names{i}] = parse_cell_type(datarun, cell_specs{i});
    end
    
    if numel(cell_specs) == 1
        cell_type_names = cell_type_names{1};
    end
else
    [cell_type_nums, cell_type_names] = parse_cell_type(datarun, cell_specs);
end



function [cell_type_num, cell_type_name] = parse_cell_type(datarun, cell_type)

switch class(cell_type)
    case 'double' % cell type is specified as the cell type number
        cell_type_num = cell_type;
        cell_type_name = datarun.cell_types{cell_type_num}.name;
    
    case 'char' % cell type is specified by its name
        % look through list of cell types for a name match (case insensitive)
        for i = 1:length(datarun.cell_types)
            if strcmpi(datarun.cell_types{i}.name, cell_type)
                cell_type_num = i;
                cell_type_name = datarun.cell_types{cell_type_num}.name;
                return
            end
        end
        
        % if no match was found, give an error
        error('Cell type ''%s'' not recognized.', cell_type)
        
    otherwise
        error('Cell specification not recognized.  Type ''help get_cell_type_number''.')
end

