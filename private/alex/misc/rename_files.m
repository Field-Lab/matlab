path2data='/Volumes/Analysis/2011-10-25-0/d01-06-norefit/';

name_extension='-from-data001_data002_data003_data004_data005_data006';

add_extension='-from-d01_06';

range = 1:6;
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





