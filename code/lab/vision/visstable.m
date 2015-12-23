function visstable()
% VISSTABLE     Swap out volatile or local Vision path for stable path

unload_vislocal();
javaswappath(vision_path_stable(), vision_path_volatile());