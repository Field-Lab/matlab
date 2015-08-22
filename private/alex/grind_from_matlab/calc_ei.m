function calc_ei(date, path2analysis)

tmp = dir([path2analysis, 'data*']);


for i = 1:length(tmp)
    dat_dir = fullfile('/Volumes/Data/', date, tmp(i).name);
    an_dir = fullfile(analysis_path,  tmp(i).name);
    my_command = ['/Volumes/Lab/Development/scripts/grind -o ', dat_dir, ' ', an_dir];
    system(my_command);
end