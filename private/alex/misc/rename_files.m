path2data='/Volumes/Analysis/2015-08-17-5/d24-48-norefit/';

name_extension='-from-data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039_data040_data041_data042_data043_data044_data045_data046_data047_data048';

add_extension = '';
% add_extension='-from-d04_20';

range = [24:48];
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


addpath(genpath('/home/vision/alex/matlab/private/alex'))
addpath(genpath('/Users/alexth/test4/matlab/private/alex'))

% calc_ei_grind('2015-08-17-1', '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/')
calc_ei_grind('2015-08-17-1', '/Volumes/Analysis/2015-08-17-1/d18-40-norefit/')

calc_ei_grind('2015-08-17-5', '/Volumes/Analysis/2015-08-17-5/d01-29-norefit/')
calc_ei_grind('2015-08-17-5', '/Volumes/Analysis/2015-08-17-5/d24-48-norefit/')

% calc_sta_grind('/Volumes/Data/2015-08-17-1/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-1/d18-40-norefit/')
% calc_sta_grind('/Volumes/Data/2015-08-17-5/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-5/d01-29-norefit/')
calc_sta_grind('/Volumes/Data/2015-08-17-5/Visual/stimuli.lisp', '/Volumes/Analysis/2015-08-17-5/d24-48-norefit/')



