function downloadJob(sched,job,verbose,forceNoZip)
%Download the job (assumed to be completed) specified to the local dataDir

if nargin < 3
    verbose = 0;
    forceNoZip=0;
end

%Retrieve properties from sched and Job
userData = sched.UserData;
clusterHost = userData.ClusterHost;
remoteDataLocation = userData.RemoteDataLocation;
schedData = getJobSchedulerData(sched,job);

localDataLocation = sched.DataLocation;
jobLocation = job.pGetEntityLocation;
remoteDataLocation = userData.RemotePath;

if verbose
    fprintf('Attempting to download %s from %s\n', get(job,'Name'),remoteDataLocation);
end

try
    %Copy the zipped file to localDataLocation
    %Find out if we are doing a zip transfer
    fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
    transferProps = fm.getTransferProperties();
    %ALL proxies MUST contain this key
    if transferProps.containsKey('mode')
        mode = transferProps.get('mode');
        %forceNoZip must be False in order to try zip download.
        if ~strcmp(char(mode),'file') && ~forceNoZip
            remoteFile = schedData.ZipLocation;
            %remoteFile = transferProps.get('zipFile');
            %Create a local filename for this bad boy
            lFile = java.io.File(localDataLocation);
            if verbose
                fprintf('Downloading main zip file...\n');
            end
            fm.get(lFile,remoteFile);
            jn = get(job,'Name');
            if verbose
                fprintf('Downloading Task1.out.mat...\n');
            end
            fm.get(java.io.File(lFile,jn), [remoteDataLocation '/' jn '/Task1.out.mat']);
        else
            if verbose
                fprintf('Copying State file...\n');
            end
            iCopyJobStateFile(localDataLocation, remoteDataLocation, ...
                jobLocation, clusterHost);
            if verbose
                fprintf('Copying Job Directory...\n');
            end
            iCopyJobDirectory(localDataLocation, remoteDataLocation, ...
                jobLocation, clusterHost);
            if verbose
                fprintf('Copying Job Files...\n');
            end
            iCopyJobFiles(localDataLocation, remoteDataLocation, ...
                jobLocation, clusterHost);
        end
    else
        error('FileMover appears to have a bad TransferProxy.  Add a valid TransferProxy');
    end
    % restore schedData
    setJobSchedulerData(sched,job,schedData);
catch err
    error('distcomp:genericscheduler:GetJobOutput', ...
        ['Unable to access files in directory %s on host %s', ...
        'because of error\n%s\n', ...
        'You will need to manually copy files from %s\n', ...
        'on host %s\n', ...
        'to the local directory %s\n', ...
        'To stop seeing this message, cancel Job %s.'
        ], remoteDataLocation, ...
        clusterHost, err.message, remoteDataLocation, clusterHost, ...
        localDataLocation, num2str(job.ID));
end

function iCopyJobStateFile(localDataLocation, remoteDataLocation, ...
    jobLocation, clusterHost)
remoteJobStateFile = [ remoteDataLocation '/' jobLocation '.state.mat' ];
copyDataFromCAC(localDataLocation, remoteJobStateFile, clusterHost);


function iCopyJobDirectory(localDataLocation, remoteDataLocation, ...
    jobLocation, clusterHost)
remoteJobDirectory = [ remoteDataLocation '/' jobLocation ];
copyDataFromCAC(localDataLocation, remoteJobDirectory, clusterHost);


function iCopyJobFiles(localDataLocation, remoteDataLocation, ...
    jobLocation, clusterHost)
remoteJobFiles = [ remoteDataLocation '/' jobLocation '.*' ];
copyDataFromCAC(localDataLocation, remoteJobFiles, clusterHost);