function cacrun_bigdata(remoteDir)
% the basic flow of operation:
% 1. create temporary directory on node
% 2. copy work file to temporary directory on node
% 3. save output to node temp directory
% 4. copy output back to remote home directory
% 5. clean-up temp directory on node

% change to temp directory
% unique name to avoid clobbering
UUID = char(java.util.UUID.randomUUID);
TMPDIR = ['t:\' UUID];
mkdir(TMPDIR);
orig_dir = pwd;

cd(TMPDIR);

% copy data file to node for reading

copyfile(fullfile(remoteDir,'bigdata.mat'),TMPDIR);
% load data into memory
load('bigdata.mat','bigdata');


% silly operation 
% in a real example this would be some sort of computation
reallybigdata = bigdata * 2;

% save output data to node
save('reallybigdata.mat','reallybigdata'); 


% copy back to shared file-system 
copyfile('reallybigdata.mat',fullfile(remoteDir,'reallybigdata.mat'));
% clean-up temp files
delete('bigdata.mat');
delete('reallybigdata.mat');
% leave directory before deleting
cd(orig_dir);
rmdir(TMPDIR);
