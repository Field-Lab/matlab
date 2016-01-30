path2data='/Volumes/Analysis/2011-10-25-9/d02-10-norefit/';

name_extension='-from-data002_data003_data004_data005_data006_data007_data008_data009_data010';

add_extension = '';
% add_extension='-from-d04_20';

range = [2:10];
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
%        elseif ~isempty(regexp(a(j).name, '.ei', 'once'))
%            movefile([mydir,'/',a(j).name], [mydir,'/',myNewName, '.ei'])
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


addpath(genpath('/home/vision/alex/matlab/private/alex'))
addpath(genpath('/Users/alexth/test4/matlab/private/alex'))

%  calc_ei_grind('2015-03-09-2', '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/')

% calc_ei_grind('2015-08-17-1', '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/')
% calc_ei_grind('2015-08-17-1', '/Volumes/Analysis/2015-08-17-1/d18-40-norefit/')

% calc_ei_grind('2015-08-17-5', '/Volumes/Analysis/2015-08-17-5/d01-29-norefit/')
% calc_ei_grind('2015-08-17-5', '/Volumes/Analysis/2015-08-17-5/d24-48-norefit/')

% calc_sta_grind('/Volumes/Data/2015-08-17-1/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-1/d18-40-norefit/')
% calc_sta_grind('/Volumes/Data/2015-08-17-5/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-5/d01-29-norefit/')
% calc_sta_grind('/Volumes/Data/2015-08-17-5/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-5/d24-48-norefit/')
% calc_sta_grind('/Volumes/Archive/2007-02-06-4/Visual/stimuli.lisp', '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/')
calc_sta_grind('/Volumes/Data/2011-10-25-9/Visual/stimuli.lisp', '/Volumes/Analysis/2011-10-25-9/d02-10-norefit/', 60.35)
% 
calc_sta_grind('/Volumes/Archive/2007-02-06-4/Visual/stimuli.lisp', '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/')
get_wn_movie_names('/Volumes/Data/2011-10-25-9/Visual/stimuli.lisp', 60.35);
calc_ei_grind('2007-02-06-4', '/Volumes/Analysis/2007-02-06-4/d00-015-norefit/')
