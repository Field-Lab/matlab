function nyc = map2manhattan(map, varargin)
% MAP2MANHATTAN     Convert an indexed map to pixel edge lines (Manhattan grid style)

opts = inputParser();
opts.addParamValue('indices', sort(setdiff(unique(map), 0)));
opts.addParamValue('quiet', false);
opts.parse(varargin{:});
opts = opts.Results;

if ischar(map), map = load_stimpixmap(map); end

if ~opts.quiet
    for i = 1:length(opts.indices)
        fprintf('=');
    end
    fprintf('\n');
end

nyc = cell(length(opts.indices),1);
for i = 1:length(opts.indices)
    if ~opts.quiet
        fprintf('*');
    end
    
    [x,y] = mask2manhattan(map == opts.indices(i));
    nyc{i}.x = x;
    nyc{i}.y = y;
end
if ~opts.quiet, fprintf('\n'); end