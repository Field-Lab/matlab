function datarun = get_lazyhex_ei_neighbors(datarun, varargin)

opts = inputParser();
opts.addParamValue('pitch', []);
opts.parse(varargin{:});
opts = opts.Results;


if isempty(opts.pitch)
    ai = load_array_info(datarun.ei.java_ei.getArrayID());
    opts.pitch = ai.spacing;
end


if ~isfield(datarun.ei, 'neighbors')
    datarun = get_ei_neighbors(datarun);
end


for i = 1:length(datarun.ei.neighbors)
    datarun.ei.neighbor_struct(i) = get_lazyhex_ei_neighbors(i, datarun.ei.neighbors{i}, datarun.ei.position, 'pitch', opts.pitch);
end