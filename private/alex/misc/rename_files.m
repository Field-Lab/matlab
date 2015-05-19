path2data='/Volumes/Analysis/2010-09-24-1/d05-36-norefit/';

name_extension='-from-data005_data006_data007_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036';

add_extension='-from-d05_36';

range = [5:7 21:36];
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





