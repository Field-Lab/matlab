disp('Running startup...');

matlabPath = '/snle/home/snl-e/matlab-standard/code';
addpath([matlabPath '/lab/utilities'])
addpath(genpath2(matlabPath, {'.svn'}));
addpath('/snle/home/snl-e/matlab-standard/specific_path_code/matlab_xunit_2.0.1/matlab_xunit/xunit');
addpath(genpath2('/snle/home/snl-e/matlab-standard/private/phli', {'.svn' '.git'}));
addpath(genpath2('/snle/home/snl-e/matlab-standard/private/freeman', {'.svn' '.git'}));

phlMatlabPath = '/Users/peterli/Dropbox/matlab';
addpath(genpath2(phlMatlabPath, {'.svn' '.git'}));
addpath('/Users/peterli/MATLAB', '-begin');

addpath(genpath('/Users/alexth/test4/matlab/private/alex'))

%Local Java builds: Vision, Cell-Finder
vislocal();

%For CAC Cornell/Purdue Matlab Teragrid
java.lang.System.setProperty('sun.security.ssl.allowUnsafeRenegotiation','true');

%For loading TrakEM2 project objects
trakem2_jinit();

clear matlabPath phlMatlabPath;
