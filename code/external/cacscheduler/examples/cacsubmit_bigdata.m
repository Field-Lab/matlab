function status = cacsubmit_bigdata(scheduler)
script_home = fileparts(which('cacsubmit_bigdata'));
run_script = fullfile(script_home,'cacrun_bigdata.m');
remoteDir = scheduler.SubmitFcn{3};
fprintf('creating large data array\n')

% this is a dummy array, in practice you would want
% to upload your data file in advance to your remote
% directory 

bigdata = 1:1000;
save('bigdata.mat','bigdata'); 

ftp = gridFTP();
fprintf('uploading large data array\n')
ftp.put(fullfile(pwd,'bigdata.mat'),fullfile(remoteDir,'bigdata.mat'));

% in practice the job submission script would already be on
% remote file system

fprintf('uploading job script\n')
ftp.put(run_script,fullfile(remoteDir,'cacrun_bigdata.m'));
ftp.close();
j = createJob(scheduler);
j.PathDependencies = {remoteDir};
t = createTask(j,@cacrun_bigdata,0,{remoteDir});
fprintf('job is being submitted (do not cancel)\n')
j.submit();
fprintf('waiting for %s to finish (ok to cancel)\n',j.name)
wait(j);
% the job script copied the data from the matlab node to the
% remote home directory
ftp.get(fullfile(pwd,'reallybigdata.mat'),fullfile(remoteDir,'reallybigdata.mat'));
fprintf('loading data file\n');
load('reallybigdata.mat','reallybigdata');

% verify that we have the expected result

if bigdata * 2 == reallybigdata
    fprintf('success\n');
    status = 'success';
else
    fprintf('failure\n');
    status = 'failure';
end
end
