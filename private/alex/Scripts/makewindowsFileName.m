function ouName = makewindowsFileName(inPath,driveletter)

if nargin == 1
    driveletter = 'S:';
end

ouName = strrep(strrep(inPath,'/mnt/muench_data',driveletter),'/','\');