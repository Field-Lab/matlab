function [j,filePath] = cacMexSubmitter(sched)

sendFileToCAC('cacMex.m');

%Move all the files to be compiled over, 
%see also gridFTP
sendFileToCAC('files\block.h');
sendFileToCAC('files\energy.h');
sendFileToCAC('files\graph.cpp');
sendFileToCAC('files\graph.h');
sendFileToCAC('files\maxflow.cpp');
sendFileToCAC('files\graphcut_denoising_new.cpp');
sendFileToCAC('files\graphcut_denoising.h');

%Assemble the compile file list.  These are the source files that will need to be on
%the mex line.
fileList = {'graphcut_denoising_new.cpp','graph.cpp','maxflow.cpp'};
compileFiles = strvcat(fileList);
%Assemble the full file list, this includes all additional header files and include files.
fileList = {'graphcut_denoising_new.cpp','graph.cpp','maxflow.cpp',...
    'block.h','energy.h','graph.h','graphcut_denoising.h'};
allFiles = strvcat(fileList);

%Create our single task job.
j = createJob(sched);
createTask(j,@cacMex,1,{allFiles,compileFiles});

alltasks = get(j, 'Tasks');
set(alltasks, 'CaptureCommandWindowOutput', true);
a = get(sched,'ParallelSubmitFcn');
set(j,'PathDependencies',{a{3}});

submit(j);
waitForState(j);
a = getAllOutputArguments(j);
%a now contains the fullpath to the file, gridFTP that back to pwd
ftp = gridFTP();
[path,name,ext] = fileparts(a{1});
fn = [name ext];
ftp.get(fn,a{1});
filePath = fullfile(pwd,fn);

