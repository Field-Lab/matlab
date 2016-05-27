clear
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
dsave = '/Volumes/Lab/Users/Nora/Ella_Data/2016-02-17-1/';

for cell_type = {'On','Off'}
    cells = get_cell_ids(datarun_class, [cell_type{1} ' Parasol']); % cell ids to fit
    datarun_class = load_sta(datarun_class, 'load_sta', cells);
    
    %% WN
    test_data = 'data027';
    fit_data = 'data026';
    test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
    TestData = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1', 'seed', 22222);
    fit_datarun = load_data([Analysis_Path '/' fit_data '/' fit_data], struct('load_neurons', 1, 'load_params', 1));
    monitor_refresh = 100/median(diff(fit_datarun.triggers));
    
    for i = 1:length(cells)
        disp(i)
        glm_cellinfo.cid           = cells(i);
        cell_savename = [cell_type{1} 'Par_' num2str(cells(i))];
        master_idx         = find(fit_datarun.cell_ids == cells(i));
        fitspikes = align_spikes_triggers(fit_datarun.spikes{master_idx}, fit_datarun.triggers, 100, monitor_refresh);
        figure;
        %[STA,center] = STA_Test(fitspikes, WNStimData.FitMovie, 1, 1/monitor_refresh);
        eval(['WNCellData.' cell_savename '.FitSpikes = fitspikes;'])
        eval(['WNCellData.' cell_savename '.STA = datarun_class.stas.stas{master_idx};'])
        eval(['WNCellData.' cell_savename '.TestSpikes = TestData.testspikes(:,i);'])
    end
end

% Load stim
xml_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-8-1-0.48-11111.xml';
WNStimData.FitMovie = get_WN_movie(xml_file, fitframes);
WNStimData.testmovie = TestData.testmovie;
WNStimData.monitor_refresh = monitor_refresh;

save([dsave 'WNStimData_201602171.mat'],'WNStimData');
save([dsave 'WNCellData_201602171.mat'],'WNCellData');
%}
