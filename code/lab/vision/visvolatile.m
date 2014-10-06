function visvolatile()
% VISVOLATILE     Swap out stable or local Vision path for volatile path

unload_vislocal();
javaswappath(vision_path_volatile, vision_path_stable);