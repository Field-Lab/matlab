function calc_sta_grind(path2stimfile, path2analysis, frequency)

[wn_movie_name, stix_size] = get_wn_movie_names(path2stimfile, frequency);
tmp = dir([path2analysis, 'data*']);

for i=1:length(tmp)

    datapath = fullfile(path2analysis,  tmp(i).name);
    movie_id = str2num(tmp(i).name(end-2:end))+1;
    
    if ~isempty(wn_movie_name{movie_id}) % calculate sta        
        if stix_size(movie_id)<3
            config_file = 'primate-1cone_ath.xml';
            my_command = ['/Volumes/Lab/Development/scripts/grind -p -c /Volumes/Lab/Development/vision-xml/current/', config_file,...
            ' ', datapath, ' ', wn_movie_name{movie_id}];
        else
            config_file = 'primate_ath.xml';
            my_command = ['/Volumes/Lab/Development/scripts/grind -l -c /Volumes/Lab/Development/vision-xml/current/', config_file,...
                ' ', datapath, ' ', wn_movie_name{movie_id}];
        end
    else % make params file
        my_command = ['/Volumes/Lab/Development/scripts/vision-calc-grind ',...
            '''Make Parameters File'' ''MainFilePath::/', datapath, ''' ',...
            '''nThreads::10'' ''STAFitCalculator::false'' ''TimeCourseCalculator::false'' ''AutoCalculator::false'''];
    end
    system(my_command);
end
