%% An index file needs to be made in order for this script to work.
% you can look at the index file for plantain (2008-08-27-5) for a template

%% data files:

% apple files
datarun = load_data('2010-03-05-2', 'rf-12-apple-gf');
datarun = load_data('2010-03-05-2', 'rf-13-apple-gf');
datarun = load_data('2010-03-05-2', 'rf-15-apple-gf');

% plantain file
datarun = load_data('2008-08-27-5', 'rf-1-gf');

%blueberry
datarun = load_data('2008-08-26-2', 'rf-1-blueberry');


%% load info from neurons params and sta files
% get sta summaries
tic
datarun = load_index(datarun);
%datarun = load_data(datarun_spec);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun.stimulus = guess_stimulus(datarun.stimulus);
datarun = get_sta_summaries(datarun,'all','verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct('strength','vector length','filter',fspecial('gauss',15,0.7),'thresh',5));


% load java movie
datarun = load_java_movie(datarun);


% should use thousands of stimuli, or all
start_time = 0;
datarun = get_snls(datarun, {1,2,3,4,5},'frames',:,'start_time',start_time,'stimuli',10000,'new',false);
%save([single_cone_path 'saved/' datarun.names.nickname],'datarun','-v7.3')
toc

%% If needed: fake a RGB STA from BW

for cc = 1:length(datarun.stas.rfs)
    temp_rf = datarun.stas.rfs{cc};
    temp_rf = repmat(temp_rf, [1 1 3]);
    datarun.stas.rfs{cc} = temp_rf;
end

%%
tic
bcf = bayesian_cone_finding_loop(datarun,bcf_params);
toc

tic
bcf = bayesian_cone_finding_multicore(datarun, bcf_params);
toc

% make some false cone loations and types if needed
datarun.cones.centers = [];
datarun.cones.types = [];

choose_magic_number(datarun,bcf,bcf_params);

tic
save_bayesian_cones(datarun,bcf,bcf_params, 50, '50-test',true,[]);
toc

