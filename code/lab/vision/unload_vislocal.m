function unload_vislocal()

if exist('vision_path_local', 'file')
    javarmpath_quiet(vision_path_local);
end