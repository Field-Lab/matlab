function mexFilePath = cacMex(allFiles, compileFiles)
%CACMEX(allFiles, compileFiles) - compiles c or c++ files to mex files on
% TUC and returns resulting mexfile.  
% allFiles - a char matrix where each row corresponds to a file needed for
%  the compilation.  This should include all headers and source files needed
%  by the compilation.
% compileFiles - a char matrix where each row corresponds to a source file
%  that would be placed on the mex commandline.  The first row should
%  contain the actual mex interface file (includes mex.h).
%
% EXAMPLES:
%   fList = {'mexfile.cpp','class1.cpp'};
%   compileFiles = strvcat(fList);
%   fList ={'mexfile.cpp','class1.cpp','class1.h','class2.h'}
%   allFiles = strvcat(fList);
%   %Move all files to TUC
%   sendFileToCAC('mexFile.cpp');
%   sendFileToCAC('class1.cpp');
%   sendFileToCAC('class1.h');
%   sendFileToCAC('class2.h'');
%   %Now create a job that submits the function
%   j = createJob(sched);
%   createTask(j,@cacMex,1,{allFiles,compileFiles});
%   Grab commandwindow output incase compile files
%   set(j,'PathDependencies',{'\\storage01\matlab\naw47'});
%   submit(j);
%   waitForState(j);
%   %Output contains the compiled mex file, use gridFTP to retrieve it
%   a = getAllOutputArguments(j);
%   ftp = gridFTP();
%   [path,name,ext] = fileparts(a{1});
%   fn = [name ext];
%   ftp.get(fn,a{1});

%files is a cell array of the files needed for compilation, with the
%mex-interface file the first one on the list.  The files need to be on the
%path, but should just contain the names of the files.

homedir=pwd;

%Locate mexopts, uses Visual Studio 9 64 bit compiler.
mopts = fullfile(matlabroot,'toolbox','local','mexopts.bat');

%Create a temporary directory
tDir = tempname();
[success,message,mid] = mkdir(tDir);
if success == 0
    error('Creating temporary directory failed: %s', message);
end

%Unpack files into a cell array
allFiles = cellstr(allFiles);
compileFiles = cellstr(compileFiles);
%Move all files to the directory
for i=1:length(allFiles)
    fullpath = which(allFiles{i});
    [success,message,mid] = copyfile(fullpath,tDir);
    if success == 0
        error('Copying file to temporary directory failed: %s', message);
    end
end
cd(tDir);

%Issues the mex command, add any arguments that you need.
mexStr = ['mex -f ' mopts];
%Add all of the files on to the end of the command
fileList = ' ';
for i=1:length(compileFiles)
    fileList = [fileList ' ' compileFiles{i}];
end
mexStr = [mexStr fileList];
fprintf('mexStr %s\n',mexStr);
eval(mexStr);

%verify a mexw64 file was created, you could beef this up by checking for
%the expected file name fairly easily if you would like to.
f = ls('*.mexw64');
if length(f) == 0
    error('No compiled file created.  Check command output for errors');
end

%Copy the mex files back to the job home
jobLocation = getenv( 'MDCE_JOB_LOCATION' );
storageLocation = getenv( 'MDCE_STORAGE_LOCATION' );
jobHome = fullfile(storageLocation, jobLocation);
copyfile(f,jobHome);

%Clean the temporary directory
[success, message, messageid] = rmdir(tDir,'s');
if success == 0
    warning('Failure cleaning up temporary directory, check output:%s', message);
end

mexFilePath = fullfile(jobLocation,f);
