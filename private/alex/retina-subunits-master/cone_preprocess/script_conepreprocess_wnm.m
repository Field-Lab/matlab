% ssh snl-e@wharf.dhcp.snl.salk.edu
% /Applications/MATLAB7.13.app/bin/matlab -nodesktop
% addpath(genpath('/snle/lab/matlab-standard/'));
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
% cd /snle/home/snl-e/matlab-standard/private/freeman

data_path = '2012-09-21-2/streamed/data011/data011';
movie_xml_path = '/snle/acquisition/movie-xml2/BW-2-5-0.48-33333-300x300-60.35.xml';
cone_finding = '_snle_acquisition_2012-09-21-2_data011_data011-localmax-t11';

datarun_wn = load_data(data_path);
datarun_wn = load_params(datarun_wn);
datarun_wn = load_neurons(datarun_wn);
datarun_wn.names.movie_xml_path = movie_xml_path;
datarun_wn = load_java_movie(datarun_wn);  
datarun_wn = load_cones(datarun_wn,cone_finding);
datarun_wn.stimulus.refresh_period = datarun_wn.stimulus.java_movie.getRefreshTime; 
datarun_wn = conepreprocess_wnm(datarun_wn);
datarun_wn = get_sta_summaries(datarun_wn, {1,2,3,4,5}, 'keep_stas', false, 'marks_params', struct('robust_std_method', 6));
conepreprocess_save(datarun_wn);

% cd /share/users-lcv/
% scp snl-e@wharf.dhcp.snl.salk.edu:/snle/acquisition/freeman/_date_/_run_ .