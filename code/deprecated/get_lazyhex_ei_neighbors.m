function neighbor_struct = get_lazyhex_ei_neighbors(electrode, neighbors, positions, varargin)

opts = inputParser();
opts.addParamValue('pitch', []);
opts.parse(varargin{:});
opts = opts.Results;


if isempty(opts.pitch)
    min_distance = ipdm(positions, 'Subset', 'SmallestFew', 'Limit', 1);
    opts.pitch = full(min_distance(min_distance > 0));
end


% For lazy-hexagonal, some neighbors will be separated by PITCH, some by sqrt(5) * PITCH
distances = ipdm(positions([electrode; neighbors],:));
inline = distances(2:end,1) < (opts.pitch*1.05);
neighbor_struct.inline = neighbors(inline);
neighbor_struct.acrosslines = neighbors(~inline);