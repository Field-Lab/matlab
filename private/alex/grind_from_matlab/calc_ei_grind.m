function calc_ei_grind(date, path2analysis)

tmp = dir([path2analysis, 'data*']);


for i = 1:length(tmp)
    dat_dir = fullfile('/Volumes/Data/', date, tmp(i).name);
    if ~isdir(dat_dir)
        dat_dir = fullfile('/Volumes/Archive/', date, tmp(i).name);
    end
    an_dir = fullfile(path2analysis,  tmp(i).name);
    if ~exist(fullfile(an_dir, [tmp(i).name,'.ei']), 'file')        
        my_command = ['/Volumes/Lab/Development/scripts/grind -o ', dat_dir, ' ', an_dir];
        system(my_command);
    end
end