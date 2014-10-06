function map = load_stimpixmap(varargin)

if nargin == 1
    filename = varargin{1};
else
    piece = varargin{1};
    mapname = varargin{2};
    mapnum = 0;
    if nargin > 2, mapnum = varargin{3}; end
    
    filename = sprintf('%s%s/%s/map-%04d.txt', server_data_path, piece, mapname, mapnum);
end

map = dlmread(filename);