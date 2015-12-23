function pth = get_rrs_neurons_path(datarun)
% Get rss_neurons_path from one of several locations in DATARUN
% 2010-01 phli

if isfield(datarun.names, 'rrs_neurons_path') && ~isempty(datarun.names.rrs_neurons_path)
    pth = datarun.names.rrs_neurons_path;
    return;
end

if all({'cell_file', 'experiment'}, @(field) isfield(datarun.names, field) && ~isempty(datarun.names.(field)))
    splitfname = split(datarun.names.cell_file, '.');
    ext = splitfname{end};
    if strcmp(ext, 'neurons')
        pth = [server_path datarun.names.experiment '/' datarun.names.cell_file];
        return;
    end
end

pth = [];