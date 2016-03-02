clear
close all




[~, txt,array_size] = xlsread('/Volumes/Lab/Users/crhoades/Large Cell Data ARVO.xlsx');
array_size = cell2mat(array_size(2:end,4));
% txt = {'
% };
for j= 25:size(txt,1)
    clear dataparam fitparam
    % piece = txt(j,1:3);
    dataparam.date=[strtrim(txt{j,1}), '/']; % one slash at the end
    temp1 =  strtrim(txt{j,3});
    stimulus =  strtrim(txt{j,6});
    dataparam.concatname=temp1; % Name (or modified name) of run, no slashes\
    
    % Sometimes the data has two versions of the concate name
    % run_opts.file_name = [run_opts.date, '/', run_opts.concatname, '/',  run_opts.concatname];
    
    
    % Wrong Movie Information
    dataparam.file_name_wrong = [dataparam.date, '/', dataparam.concatname, '/', 'wrongMovie/wrongMovie'];
    
    wrong_stimulus = [stimulus(1:end-5), num2str(str2num(stimulus(end-4:end))+11111)];
    dataparam.mdf_file_wrong=['/Volumes/Analysis/stimuli/white-noise-xml/', wrong_stimulus, '.xml'];
    
    % Right Movie Information
%     dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
    %
        dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname, '/',dataparam.concatname];

    dataparam.mdf_file_right=['/Volumes/Analysis/stimuli/white-noise-xml/', stimulus, '.xml'];
    
    
    % fit parameters
    fitparam.fit_n_one_filters = true;
    fitparam.fit_n_two_filters = true;
    fitparam.fit_surround_sd_scale = false;
    fitparam.fit_surround =false;
    fitparam.fit_surround_amp_scale =false;
    fitparam.initial_n_one_filters =12;
    fitparam.initial_n_two_filters = 12;
    
    % fitparam.interpolate = false;
    % fit the channels separately (useful for blue/green or SBCs)
    fitparam.independent_fit = 0;
    fitparam.num_frames = 30; % both have to be run with the name number of frames
    % how many stixels per STA on average can be actually noise
    fitparam.false_stixels =0.25;
    
    
    dataparam.cell_type = {'ON large 1'};
    
    % list specific cell (1), or run for a whole cell type (0)
    dataparam.select_cells = 0;
    dataparam.cell_specification = [6512]; %can be anything if select_cells =0
    
    
    dataparam.filepath=['/Volumes/Lab/Users/crhoades/Fitting/',dataparam.date,'/',dataparam.concatname,'/'];
    
    fit_large_cell_sta_function(dataparam, fitparam)
end
