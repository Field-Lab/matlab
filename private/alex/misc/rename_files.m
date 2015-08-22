path2data='/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';

name_extension='-from-data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015_data016_data017_data018_data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029';

add_extension = '';
% add_extension='-from-d04_20';

range = [1:29];
% rename files
for i=range
   tmp=['00', int2str(i)];
   myName=['data', tmp(end-2:end), name_extension];
   myNewName=['data', tmp(end-2:end), add_extension];
%    mydir=[path2data, myNewName];
   mydir=[path2data, myName];
   a=dir(mydir);
   for j=1:length(a)
       [b,c]=regexp(a(j).name, myName);
       if ~isempty(b)
           movefile([mydir,'/',a(j).name], [mydir,'/',myNewName, a(j).name(c+1:end)])
       elseif ~isempty(regexp(a(j).name, '.ei', 'once'))
           movefile([mydir,'/',a(j).name], [mydir,'/',myNewName, '.ei'])
       end
   end
end

%rename folders

for i=range
   tmp=['00', int2str(i)];
   myName=['data', tmp(end-2:end), name_extension];
   myNewName=['data', tmp(end-2:end), add_extension];
   movefile([path2data, myName], [path2data, myNewName])
end

date = '2015-08-17-1';
analysis_path = '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';
calc_ei('2015-08-17-1', '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/')

calc_sta('/Volumes/Data/2015-08-17-1/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/')

filepath = '/Volumes/Data/2015-08-17-1/Visual/stimuli.lisp';
[wn_movie_name, stix_size] = get_wn_movie_names(filepath);

addpath(genpath('/home/vision/alex/matlab/private/alex'))


for i=37
    datapath = ['/Volumes/Analysis/2015-03-09-2/d37-40-norefit-test/data0', int2str(i)];
    
    if ~isempty(wn_movie_name{i}) % calculate sta
        if stix_size(i)<3
            config_file = 'primate-1cone_ath.xml';
            my_command = ['/Volumes/Lab/Development/scripts/grind -p -c /Volumes/Lab/Development/vision-xml/current/', config_file,...
            ' ', datapath, ' ', wn_movie_name{i}];
        else
            config_file = 'primate_ath.xml';
            my_command = ['/Volumes/Lab/Development/scripts/grind -l -c /Volumes/Lab/Development/vision-xml/current/', config_file,...
                ' ', datapath, ' ', wn_movie_name{i}];
        end
    else % make params file
        my_command = ['/Volumes/Lab/Development/scripts/vision-calc-grind ',...
            '''Make Parameters File'' ''MainFilePath::/', datapath, ''' ',...
            '''nThreads::10'' ''STAFitCalculator::false'' ''TimeCourseCalculator::false'' ''AutoCalculator::false'''];
    end
    system(my_command);
end

% for ei calculation, use:
date = '2015-08-17-1';
analysis_path = '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';
for i = 1:10
    if i<10
        j = ['0', int2str(i)];
    else
        j = int2str(i);
    end
    dat_dir = ['/Volumes/Data/', date, '/data0', j];
    an_dir = [analysis_path, '/data0', j];
    my_command = ['/Volumes/Lab/Development/scripts/grind -o ', dat_dir, ' ', an_dir];
    system(my_command);
end
