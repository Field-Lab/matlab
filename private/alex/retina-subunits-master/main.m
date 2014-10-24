
    
%% cone preprocessing

datarun = load_data(fullfile(server_path(), '2012-09-24-1/data003/data003'));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones(datarun,'15'); % load_cones(datarun, 'Analysis');
datarun = make_mosaic_struct(datarun);

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'bayes');
    
%% cone preprocessing

datarun = load_data(fullfile(server_path(), '2011-06-30-0/streamed/data003/data003'));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones(datarun,'10'); % load_cones(datarun, 'Analysis');
datarun = make_mosaic_struct(datarun);

datarun.names.rrs_movie_path='/Volumes/Analysis/2011-06-30-0/data003/data003.movie';%BW-2-4-0.48-11111-300x300-60.35.xml';

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'bayes');

    
%%
dat = loadData('stanford','2011-06-30-0/data003');

batchAnal('stanford','2011-06-30-0/data003','off midget','stc',1,0,0.33);

batchAnal('stanford','2011-06-30-0/data003','off midget','subunit',1,0,0.33)

batchAnal('stanford','2011-06-30-0/data003','off midget','subunit',0,1,0.33)


%%
dat = loadData('stanford','2011-10-25-5/data001-0');

batchAnal('stanford','2011-10-25-5/data001-0','off midget','stc',1,0,0.33);

batchAnal('stanford','2011-10-25-5/data001-0','off midget','subunit',1,0,0.33)

batchAnal('stanford','2011-10-25-5/data001-0','off midget','subunit',0,1,0.33)

%%


dat = loadData('stanford','2012-09-24-1/data003');

batchAnal('stanford','2012-09-24-1/data003','off midget','stc',1,0,0.33);

batchAnal('stanford','2012-09-24-1/data003','off midget','subunit',1,0,0.33)

batchAnal('stanford','2012-09-24-1/data003','off midget','subunit',0,1,0.33)





%% extra stuff
% load_sta
% fill in missing stimulus parameters, if desired
if params.guess_stimulus
    frame_rate=datarun.stimulus.interval*(1000/datarun.stimulus.refresh_period);
    if frame_rate<70
        display_type='oled';
    else
        display_type='crt';
    end
    datarun.stimulus = guess_stimulus(datarun.stimulus,'display_type',display_type);
end

% load_java_movie
    % identify the stimulus in the java movie, and guess about parameters that were not specified
    
    frame_rate=datarun.stimulus.interval*(1000/datarun.stimulus.refresh_period);
    if frame_rate<70
        display_type='oled';
    else
        display_type='crt';
    end
    
    the_stimulus = guess_stimulus(stimulus_from_java_movie(movie),'display_type',display_type);
    stimsource = sprintf('XML movie %s',movie_xml_path);


