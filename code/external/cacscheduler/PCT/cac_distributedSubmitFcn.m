function cac_distributedSubmitFcn(scheduler, job, props, ...
                                        clusterHost, remoteDataLocation,language)
% CACNONSHAREDSIMPLESUBMITFCN() - PCT interface function to submit a
% distributed job to the TUC cluster at CAC.  The end user should not call this
% function directly, it is called automatically by creating a cac scheduler
% object and caling createJob/submit.
% 
% EXAMPLES
%  cacsched - create sched object
%  j = createJob(sched);
%  t =  createTask(j, @rand, 1, {42});
%  t2 = createTask(j,@rand,1,{25});
%  j.submit(); 
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

                                   
% if ~ischar(clusterHost)
%     clusterHost='tuc0001.cac.cornell.edu';
%     %error('distcomp:genericscheduler:SubmitFcnError', ...
%     %    'Hostname must be a string');
% end
% if ~ischar(remoteDataLocation)
%     remoteDataLocation='\\storage01\matlab\xxx',
%     %error('distcomp:genericscheduler:SubmitFcnError', ...
%     %    'Remote Data Location must be a string');
% end


if nargin==5
    language='jsdlccs';
end


%scheduler.UserData = { clusterHost ; remoteDataLocation };
localDataLocation = scheduler.DataLocation;

% Set the name of the decode function which will be executed by
% the worker. The decode function must be on the path of the MATLAB
% worker when it starts up. This is typically done by placing the decode
% function in MATLABROOT/toolbox/local on the cluster nodes, or by
% prefixing commandToRun (created below) with a command to cd to the
% directory where the decode function's M file exists.
decodeFcn = 'cac_distributedDecodeFcn';

% Read the number of tasks which are to be created. This property
% cannot be changed.
numberOfTasks = props.NumberOfTasks;

% A unique file and directory name for the job. This is used to create
% files and a directory under scheduler.DataLocation
jobLocation = props.JobLocation;

% A cell array of unique file names for tasks. These are used to create
% files under jobLocation
taskLocations = props.TaskLocations;

% Since PBS jobs will be submitted from a UNIX host on the cluster, 
% a single quoted string will protect the MATLAB command.
quote = '''';

% The MATLAB command to be run on a cluster node to execute a task.
matlabexe = props.MatlabExecutable;
matlabargs = props.MatlabArguments;
%commandToRun = [ props.MatlabExecutable ' ' props.MatlabArguments ];
%get(props);

remPath = scheduler.UserData.RemotePath;
userHome = scheduler.UserData.UserHome;

%Find out if we are doing a zip transfer
fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
transferProps = fm.getTransferProperties();
%ALL proxies MUST contain this key
if transferProps.containsKey('mode')
    mode = transferProps.get('mode');
    if strcmp(char(mode),'file')
        gangTransfer=false;
    else
        gangTransfer=true;
    end
else
    error('FileMover appears to have a bad TransferProxy.  Add a valid TransferProxy');
end

% Copy the matlab_metadata.mat file to the remote host.
%FilesToZip = {};
localMetaDataFile = [ localDataLocation '/matlab_metadata.mat' ];
%FilesToZip{end+1} = localMetaDataFile;
copyDataToCAC(localMetaDataFile, remPath, clusterHost);

% Copy the local job directory to the remote host.
localJobDirectory = [ localDataLocation '/' jobLocation ];
%Can zip handle a directory?
%FilesToZip{end+1} = localJobDirectory;
copyDataToCAC(localJobDirectory, remPath, clusterHost);

% Copy the local job files to the remote host.
localJobFiles = [ localDataLocation '/' jobLocation '.*' ];
%FilesToZip{end+1} = localJobFiles;
copyDataToCAC(localJobFiles, remPath, clusterHost);
%All file transfers have been requested, so if gang, signal
if gangTransfer
    fm.setTransferProperty('jobID',get(job,'Name'));
    fm.isTransferDone();
    %pause(5);
end
%zip it up, but DON't prefix with jobLocation or wildcards will see it!
%zipName = ['LJZIP_' jobLocation '.zip'];
%zipFileName = fullfile(localDataLocation, zipName);
%zip(zipFileName,FilesToZip);
%Create the remote job directory so that ccs errors can go somewhere
%ftp = gridFTP();
%ftp.createDir(jobLocation);
%remoteZipLocation = [remoteDataLocation '\' zipName];
%ftp.put(zipFileName,zipName);

% Submit tasks which the scheduler will execute by starting MATLAB workers.
jobinfoarray= cell(1,numberOfTasks);
for i = 1:numberOfTasks
    taskLocation = taskLocations{i};
    remoteJobDir = [ remoteDataLocation '\' jobLocation ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We pull some properties from the static ClusterInfo class to
    % configure the job.
    cacprops = getClusterInfoProperties();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This is where CAC diverges from (eg) PBS.
    %There is no batch script, instead we have to provide the environment
    %settings, commandToRun(in parts), and stdout/stderr locations
    
    %taskLocation appears to come in like this: JobNN/TaskN, which is a
    %backwards slash, just ensure everything looks OK:
    taskLocation = strrep(taskLocation,'/','\');
    %Create environment variables, we'll need 6:
    %MDCE_DECODE_FUNCTION,MDCE_STORAGE_LOCATION,MDCE_STORAGE_CONSTRUCTOR
    %MDCE_JOB_LOCATION,MDCE_TASK_LOCATION,MDCE_DEBUG
    %we kind of have to hack this as there isn't a great way to do a map
    %+2 for unzipping
    env = cell(6,1);
    env{1} = sprintf('MDCE_DECODE_FUNCTION:%s',decodeFcn);
    env{2} = sprintf('MDCE_STORAGE_LOCATION:%s',remoteDataLocation);
    env{3} = sprintf('MDCE_STORAGE_CONSTRUCTOR:%s',props.StorageConstructor);
    env{4} = sprintf('MDCE_JOB_LOCATION:%s',props.JobLocation);
    env{5} = sprintf('MDCE_TASK_LOCATION:%s',taskLocation);
    env{6} = sprintf('MDCE_DEBUG:%s','TRUE');
    %If we are performing a zip operation, then we need to add to
    %properties to signal the worker
    %zipFile MUST have a value
    zipFile = '';
    if gangTransfer
        env{7} = sprintf('CACSCHEDULER_ZIP_JOB:%s','TRUE');
        %Retrieve from props
        transferProps = fm.getTransferProperties(); 
        %HACK
        %The transfer proxy creates all paths in the zip file appropriatly
        %but the creation of the remoteZipLocation includes those paths
        %as well.  So we hack the remoteDataLocation to jettison the path.
        zipFile = transferProps.get('zipFile');

        remoteZipLocation = [userHome '\' zipFile];

        env{8} = sprintf('CACSCHEDULER_ZIP_LOCATION:%s',remoteZipLocation);
    end
    
    
    %Construct stdout/stderr streams
    stdoutFileLocation = [remoteDataLocation '\' taskLocation '.ou' ];
    stderrFileLocation = [remoteDataLocation '\' taskLocation '.er' ];
    
    %arrayify arguments
    args = cell(1,1);
    args{1} = matlabargs;
    
    % verify upload completed
    if gangTransfer
        ftp = gridFTP();
        %fprintf('verify location\n');
        zfile = ftp.isFile(remoteZipLocation);
        if ~zfile
            error('CACSCHEDULER:JobSubmission Failure','Unable to upload zip file');
        end
    end
    
    %Start calling java
    hpc = edu.cornell.cac.tuc.matlab.JSDLMediator(1);
    try
        jobinfo = hpc.jobStart(matlabexe, args, stdoutFileLocation, ...
            stderrFileLocation,env,cacprops,1);
    catch
        fprintf('Submission of Task%d failed, retrying in 5 seconds...\n',i);
        pause(5);
        try
            jobinfo = hpc.jobStart(matlabexe, args, stdoutFileLocation, ...
            stderrFileLocation,env,cacprops,1);
        catch
            error('CACSCHEDULER:JobSubmission Failure','Unable to submit Task%d\n',i);
        end
    end
    jobinfoarray{1,i} = jobinfo;
end
%Package the data to be persisted
s.CARObjects = jobinfoarray;
s.ZipLocation = zipFile;
s.RemoteDataLocation = remoteDataLocation;
s.RemotePath = remPath;
setJobSchedulerData(scheduler, job, s);
%set(job,'UserData',jobinfoarray);
