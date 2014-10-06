function is_root_path = is_windows_root_path(pth)
is_root_path = ~isempty(regexpi(pth, '^[a-z]+:[/\\]'));