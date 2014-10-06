function state = cac_getJobStateFcn(scheduler, job, state)
%CACGetJobState() - Retrieves the state of a job running on TUC.  The end user
%should not call this function directly, but instead access these functions
%through the cac scheduler object and PCT interface functions.
%
% EXAMPLES
%  cacsched - create sched object
%  j = createJob(sched);
%  t =  createTask(j, @rand, 1, {42});
%  t2 = createTask(j,@rand,1,{25});
%  j.submit(); 
%  waitForState(j);
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing
%
%

userData = getJobSchedulerData(scheduler,job);

schedData = userData; % old reference

clusterHost = scheduler.UserData.ClusterHost;
remoteDataLocation = userData.RemoteDataLocation;
remPath = userData.RemotePath;
localDataLocation = scheduler.DataLocation;
jobLocation = job.pGetEntityLocation;
remoteDataLocation = remPath;

% Find the state of the job without copying files!
local_state = iGetLocalJobState(job);
%fprintf('local_state %s\n',local_state);

% Shortcut if the job state is already finished or failed
jobInTerminalState = strcmp(local_state, 'finished') || strcmp(local_state, ...
                                                  'failed');
% and we have already done the last mirror
if jobInTerminalState
    state = local_state;
    return;
end

hpc = edu.cornell.cac.tuc.matlab.JSDLMediator(1);

cars = userData.CARObjects;
st = {};

isrunning = 0;
ispending = 0;
isfailed = 0;

for i = 1:length(cars)
    st{end+1} = lower(char(hpc.jobStatus(cars{i})));
    %Debug print here:
    %fprintf('Job status of %s is %s.\n', char(cars{i}.getJobID()), st{end});
    if strcmp(st{end},'running')
        isrunning = 1;
    end
    if strcmp(st{end},'pending')
        ispending = 1;
    end
    if strcmp(st{end},'failed')
        isfailed = 1;
    end
    % specify to HPC job manager stopping a queued job
    if strcmp(st{end},'cancelled')
        isfailed = 1;
    end
    %     if ~( strcmp(st{end},'finished') || strcmp(st{end},'failed') )
    %         state = st{end};
    %     end
end

if all(strcmp('finished',st))
    state = 'finished';
end
if all(strcmp('failed',st))
    state = 'failed';
end
if all(strcmp('cancelled',st))
    state = 'cancelled';
end
if all(strcmp('running',st))
    state = 'running';
end
%fprintf('isrunning %d ispending %d\n',isrunning,ispending);
if isrunning == 0 && ispending == 0
    %  fprintf('CHANGE STATE TO FINISHED\n');
    if strcmp('running',state) || strcmp('queued',state) || strcmp('cancelled',state)
        %fprintf('CHANGE STATE TO FINISHED from: %s\n',state);
        if isfailed == 1
            % any type of failure means job failed
            state = 'failed';
        else
            state = 'finished';
        end
    end
end
if isrunning == 1 && ~strcmp('running',state)
    %  fprintf('CHANGE STATE TO RUNNING\n');
    state = 'running';
end
if ~strcmp(state,local_state)
    %fprintf('job state not the same hpc: %s local: %s updating state\n',state,local_state);
    job.pSetState( state );
end

%fprintf('Current state of %s is %s.\n', get(job,'name'),state);

% If the job has just finished, copy back the data from the cluster.
if strcmp(state, 'finished') || strcmp(state,'failed') || strcmp(state,'cancelled')
    fprintf('Downloading completed job: %s.\n', get(job,'name'));
    % remote Job files don't contain the schedData object
    % that was added right before the job was submitted
    % retrieve the local info and restore after download
    schedData = getJobSchedulerData(scheduler,job);
    try
        %Copy the zipped file to localDataLocation
        %Find out if we are doing a zip transfer
        fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
        transferProps = fm.getTransferProperties();
        %ALL proxies MUST contain this key
        if transferProps.containsKey('mode')
            mode = transferProps.get('mode');
            if ~strcmp(char(mode),'file')
                dAndZipStart = tic;
                remoteFile = schedData.ZipLocation;
                %remoteFile = transferProps.get('zipFile');
                %Create a local filename for this bad boy
                lFile = java.io.File(localDataLocation);
                fm.get(lFile,remoteFile);
                dAndZipStop = toc(dAndZipStart);
                downloadStart = tic;
                jn = get(job,'Name');
                fm.get(java.io.File(lFile,jn), [remoteDataLocation '/' jn '/Task1.out.mat']);
                downloadStop = toc(downloadStart);
            else
                iCopyJobDirectory(localDataLocation, remoteDataLocation, ...
                    jobLocation, clusterHost);
                iCopyJobFiles(localDataLocation, remoteDataLocation, ...
                    jobLocation, clusterHost);
            end
        else
            error('FileMover appears to have a bad TransferProxy.  Add a valid TransferProxy');
        end
        % handle aborted jobs
        if strcmp(state,'failed') || strcmp(state,'cancelled') || isfailed == 1
            %Debug print here:
            %fprintf('Process aborted, updating local state\n');
            % check local state of downloaded file
            local_state = iGetLocalJobState(job);
            if ~strcmp( local_state, 'finished' )
                % for whatever reason the job terminated strangely
                % this will either be finished or failed or cancelled from
                % HPC status
                %fprintf('updating job state to: %s\n',state);
                job.pSetState( state );
            end
            % update finish time
            serializer = job.pReturnSerializer;
            serializer.putFields(job, ...
                    {'finishtime' 'state'}, ...
                    {char(java.util.Date) 'finished'});

            tasks = job.Tasks;
            for i = 1:numel(tasks)
                if ~strcmp( tasks(i).State, 'finished' )
                    % if the task status isn't finshed set to failed
                    % this can be due to the HPC job failing before
                    % matlab could update the state
                    %fprintf('update task(%d) state to failed\n',i);
                    tasks(i).pSetState( 'failed' );
                    % set error 
                    serializer.putFields(tasks(i), ...
                    {'erroridentifier', 'errormessage', 'finishtime' 'state'}, ...
                    {'distcomp:task:Cancelled', 'job ended unexpectedly contact help@cac for further assistance', char(java.util.Date) 'finished'});
                end
            end
        end
        % restore schedData
        % sync remote files with local copies
        setJobSchedulerData(scheduler, job, schedData);
        %ftp = gridFTP();
        %zipFile = fullfile(localDataLocation, [get(job,'Name') '_out.zip']);
        %remoteFile = [remoteDataLocation '/LJZIP_' get(job,'Name') '.zip'];
        %ftp.get(zipFile,remoteFile);
        %unzip(zipFile,localDataLocation);
        %jn = get(job,'Name');
        %Should probably execute a delete on that zip as destroy won't get
        %it!
        %ftp.delete(remoteFile);
        %         %Set Task1 to finished
        %         task1State = fullfile(localDataLocation, jn, 'Task1.state.mat');
        %         fid = fopen(task1State,'w');
        %         fwrite(fid,'finished');
        %         fclose(fid);
        %         %Set Job to finished
        %         jobState = fullfile(localDataLocation, [jn '.state.mat']);
        %         fid = fopen(jobState,'w');
        %         fwrite(fid,'finished');
        %         fclose(fid);
        %retrieve task1.out.mat
        %task1Out = fullfile(localDataLocation, jn,'Task1.out.mat');
        
        %ftp.get(task1Out, [remoteDataLocation '/' jn '/Task1.out.mat']);
        
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
end
%fprintf('Download and unzip took %4.4fs\n',dAndZipStop);
%fprintf('Download only took %4.4fs\n',downloadStop);

function state = iGetLocalJobState(job)
serializer = job.pReturnSerializer;
if ~isempty(serializer)
    state = char(serializer.getField(job, 'state'));
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

