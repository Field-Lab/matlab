
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BAYESIAN CONE FINDING STEP 1 OF 5    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load a datarun
datarun_to_load = 'plantain';


% choose datarun
if ~exist('datarun','var')

    start_time = 0;
    clear old_cone_finding_path

    switch datarun_to_load
        case 'blueberry'
            load([single_cone_path 'saved/blueberry.mat'])
            start_time = 6400;
            
        case 'peach'
            load([single_cone_path 'saved/peach.mat'])
            
        case 'plantain'
            load([single_cone_path 'saved/plantain.mat'])
            
        case 'plantain-1'
            load([single_cone_path 'saved/plantain-1.mat'])
            old_cone_finding_path = 'plantain-old';
            
        case 'apricot'
            load([single_cone_path 'saved/apricot.mat'])
            start_time = 3600;
            
        case 'kiwi'
            load([single_cone_path 'saved/kiwi.mat'])
            
        case 'grapes'
            load([single_cone_path 'saved/grapes.mat'])
            
        case 'butterfly'
            load([single_cone_path 'saved/butterfly.mat'])
            
        case 'plum'
            load([single_cone_path 'saved/plum.mat'])
            
        case 'cherry'
            load([single_cone_path 'saved/cherry.mat'])
            
        case 'mango'
            load([single_cone_path 'saved/mango.mat'])
            
        case 'apple'
            load([single_cone_path 'saved/apple-13.mat'])
            
        case 'plantain_last_hour'
            load([single_cone_path 'saved/plantain_last_hour.mat'])
            old_cone_finding_path = 'plantain';
            start_time = 3601;
            
        case 'plantain_15'
            load([single_cone_path 'saved/plantain_15.mat'])
            old_cone_finding_path = 'plantain';
            
        case 'plantain_30'
            load([single_cone_path 'saved/plantain_30.mat'])
            old_cone_finding_path = 'plantain';
            
            
    end
    datarun = load_index(datarun);
    datarun = load_neurons(datarun);
    datarun = load_params(datarun);
    datarun.stas.java_sta = load_java_sta(datarun);
    datarun.stimulus = guess_stimulus(datarun.stimulus);

    % only needed if SNLs will be computed (they are usually stored in datarun already)
    if ~isfield(datarun.stas,'snls') || isempty(datarun.stas.snls)
        datarun = load_java_movie(datarun);
    end

    % load cones from old cone finding
    if ~exist('old_cone_finding_path','var')
        old_cone_finding_path = [datarun_to_load '-old'];
    end
    cone_data = import_single_cone_data([],old_cone_finding_path);
    datarun.cones.centers = cone_data.cone_centers;
    datarun.cones.types = cone_data.cone_types;
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
    %load([single_cone_path datarun.names.nickname '/Wc.mat'])
end


% to load datarun initially
if 0
   
    data_piece = '2011-06-30-0';
    data_condition = 'rf-3';
    
    datarun = load_data(data_piece, data_condition);
    datarun = load_index(datarun);
    datarun = load_params(datarun,'verbose',1);
    datarun = load_neurons(datarun);
    datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
    datarun.stimulus = guess_stimulus(datarun.stimulus);
    datarun = get_sta_summaries(datarun,'all','verbose',1,'keep_stas',0,'keep_rfs',1,'fig_or_axes',1,...
        'marks_params',struct('strength','vector length','filter',fspecial('gauss',15,0.7),'thresh',5));
    datarun = load_java_movie(datarun);
    % should use thousands of stimuli, or all
    datarun = get_snls(datarun, {1,2,3,4},'frames',-2:0,'start_time',0,'stimuli',10000,'new',false);
    save([single_cone_path 'saved/' datarun.names.nickname],'datarun','-v7.3')
end
