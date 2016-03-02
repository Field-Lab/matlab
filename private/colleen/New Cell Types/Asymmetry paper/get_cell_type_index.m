% get cell_type index

function [match] = get_cell_type_index(j_values)

cell_types_of_interest = {'ON parasol', 'OFF parasol', 'ON midget', 'OFF midget', 'ON large 1', 'ON large 2', 'OFF large 1', 'OFF large 2'};
[~, txt,array_size] = xlsread('/Volumes/Lab/Users/crhoades/Large Cell Data ARVO.xlsx');
array_size = cell2mat(array_size(2:end,4));

for j= 2:j_values%:size(txt,1)
    run_opts.date=strtrim(txt{j,1}); % one slash at the end
    temp1 = strtrim(txt{j,2});
    temp2 =  strtrim(txt{j,3});
    run_opts.concatname=temp1; % Name (or modified name) of run, no slashes\
    run_opts.dataname = temp2;
    
    % Sometimes the data has two versions of the concate name
    run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',run_opts.dataname, '/',run_opts.dataname];
    % run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
    datarun.names.rrs_params_path=['/Volumes/Analysis/', run_opts.file_name, '.params'];

    opt=struct('verbose',1,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',false, 'load_sta', 0, 'load_sta_params', 0, 'load_all',0);
    opt.load_sta_params.save_rf = 1; % Necessary to get the sta to load correctly
    %     try
    datarun=load_data(datarun,opt);
    for i = 1:size(cell_types_of_interest,2)
        for z = 1:size(datarun.cell_types,2)
            if strcmpi(cell_types_of_interest(i), datarun.cell_types{z}.name)
                match(j,i) = z;
            end
        end
        
    end
    
    clear datarun
    
end
