function struct_out = make_struct_to_pass(struct_in,priority_list)
% this function makes a struct of parameters called 'struct_out'.
% it copies values from struct_in, renaming them according to instructions encoded in priority_list.
%    if priority_list is empty, use all fields from struct_in, and keep the names identical
%    in this case, make_struct_to_pass simply serves to remove fields with the value 'default value'

if ~isempty(priority_list)
    % parse priority_list to create the variable field_lookup
    p = inputParser;
    p.KeepUnmatched = true;
    p.parse(priority_list{:});
    field_lookup = p.Unmatched;
else
    field_lookup = field_names(struct_in);
end

% field_lookup is a struct.  Each field is a potential field in struct_out, and the value
% indicates which fields of struct_in should be copied to struct_out
%
% the value of each field in field_lookup is either a string, or a cell array of strings.
%   if a string, then the string is the name of the field in struct_in whose value should be copied to struct_out
%   if a cell array, then the cell array is the ordered list of names of fields in struct_in which should be copied to struct_out.
% note that a field from struct_in is NOT copied to struct_out if the value of the field is 'default value'


% get the names of fields to potentially put in struct_out
new_fields = fieldnames(field_lookup);


% initialize struct_out 
struct_out = struct;

% fill in each potential field of struct_out, if there's anything to fill in
for ff = 1:length(new_fields)

    % get the name of the new field
    new_field_name = new_fields{ff};

    % get the name(s) of fields in struct_in whose values might be used
    source_fields = field_lookup.(new_field_name);

    % if there's only one name, put it into a cell array
    if ~iscell(source_fields)
        source_fields = {source_fields};
    end

    % look at all candidate source fields, until finding one whose value isn't 'default value'
    for ss = 1:length(source_fields)

        % if the field doesn't exist, there must be a bug in the code!
        if ~isfield(struct_in,source_fields{ss})
            error('This function must take ''%s'' as an optional argument name, but it does not!  Is something misspelled?',source_fields{ss})
        end

        % if this field has something other than 'default value'...
        if ~isequal(struct_in.(source_fields{ss}),'default value')
            % use that value in struct_out
            struct_out.(new_field_name) = struct_in.(source_fields{ss});
            % and stop searching
            break
        end
        % otherwise, the loop continues.  if no field has something other than 'default value', then the field is not present in struct_out
    end

end

