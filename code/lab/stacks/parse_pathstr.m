function [path, page] = parse_pathstr(pathstr)

splitpath = split(pathstr, ':');

if length(splitpath) == 1
    path = pathstr;
    page = 1;
else
    path = join(splitpath(1:end-1), ':');
    page = str2double(splitpath{end});
end
