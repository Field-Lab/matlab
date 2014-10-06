function parsed = parse_rrs_prefix(datarun)

parsed = struct([]);
rrs_fields = {'rrs_prefix', 'rrs_neurons_path','rrs_params_path','rrs_ei_path','rrs_sta_path','rrs_globals_path'};

if ~isfield(datarun, 'names') || isempty(datarun.names)
    return
end

for i = 1:length(rrs_fields)
    fieldname = rrs_fields{i};
    if isfield(datarun.names, fieldname) && ~isempty(datarun.names.(fieldname))
        parsed = parse_rrs_prefix(datarun.names.(fieldname));
        return
    end
end