function params_out = add_forked_params(params,p_temp)
% return list of forked parameters
%
%
%
% gauthier 2008-10
%



% make struct of parameters which
%   1) the user specified, i.e. their value is not 'default value'
%   2) are relevant to this fork, i.e. are found in p_temp.Parameters
params_in = make_struct_to_pass(params,cell_double(p_temp.Parameters));

% resolve user input and default valuess
p_temp.parse(params_in);

% add them to the params struct
params_out = struct_union(params,p_temp.Results);





function dbl_cell = cell_double(cell_in)
% double entries in a cell
% e.g. {'a','b'}  -->  {'a','a','b','b'}

for cc=1:length(cell_in)
    dbl_cell{cc*2 - 1} = cell_in{cc};
    dbl_cell{cc*2} = cell_in{cc};
end



function struct_out = struct_union(struct_in_1,struct_in_2)
% take the union of two structs, giving priority to the second if field names are duplicated

% include all fields from the first struct
struct_out = struct_in_1;

% get list of fields in second struct
add_fields = fieldnames(struct_in_2);

% incorporate them
for ff = 1:length(add_fields)
    struct_out.(add_fields{ff}) = struct_in_2.(add_fields{ff});
end
