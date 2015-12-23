function cell_type_str = load_params_cell_type(datarun, cell_ids, varargin)

opts = inputParser();
opts.addParamValue('java_params', []);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.java_params)
 opts.java_params = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);
end

cell_type_str = cell(size(cell_ids));
for i = 1:numel(cell_ids)
    cell_type_str{i} = opts.java_params.getClassIDs.get(uint32(cell_ids(i)));
end