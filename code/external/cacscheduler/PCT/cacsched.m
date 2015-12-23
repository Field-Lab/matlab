function sched= cacsched(userName);
%CACSCHED() - construct a scheduler object to interface with TUC located at
%the CAC in Ithaca, NY.  This script creates an object called 'sched' that
%is passed to createJob functions in order to initiate jobs on TUC.
%
% EXAMPLES:
%  sched = cacsched('CACUSERNAME');
%  dj = createJob(sched);
%  dj = createParallelJob(sched);
%  dj = createMatlabPoolJob(sched);
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

%% Check arguments
error(nargchk(1,1,nargin,'struct'));
0/1;
%% check java
try
    cm = edu.cornell.cac.tuc.cacscheduler.globus.CertManager.getInstance();
catch exc
    if findstr('Undefined variable "edu"',exc.message)
        error('Error loading cacscheduler.  Please verify that classpath.txt is configured.  Contact help@cac.cornell.edu for further assistance.');
    end
    error(['Unknown error:' exc.message]);
end
%Learn a bit about the system
verName = ['R' version('-release')];
ClusterMLRoot = strcat('C:\PROGRA~1\MATLAB\',verName);
CACHome = fileparts(which('runtests'));
clusterHost = 'tuc000.cac.cornell.edu';


userPath = verName; % optional sub-directory, i.e. "laptop",
                    % defaults to MATLAB version
userHome = ['\\storage01\matlab\' userName ];
RemoteDataLocation = [userHome '\' userPath];



%Extract any home relative path from remoteDataLocation
n = regexp(RemoteDataLocation,['storage01\\matlab\\(?<userid>\w+)\\' ...
                    '(?<path>.+)'],'names');
if ~isempty(n)
        remPath = ['~/' n.path];
    userHome = ['\\storage01\matlab\' n.userid];
else
        remPath = '~';
    userHome = [];
end
remPath = strrep(remPath,'\','/');


sched = findResource('scheduler','type','generic');
set(sched,'SubmitFcn',{@cac_distributedSubmitFcn,clusterHost,RemoteDataLocation});
set(sched,'ParallelSubmitFcn',{@cac_parallelSubmitFcn,clusterHost,RemoteDataLocation});
set(sched,'GetJobStateFcn',@cac_getJobStateFcn);
set(sched,'CancelJobFcn',@cac_cancelJobFcn);
set(sched,'CancelTaskFcn',@cac_cancelTaskFcn);
set(sched,'DestroyJobFcn',@cac_destroyJobFcn);
set(sched,'ClusterOsType','pc');
set(sched,'ClusterMatlabRoot',ClusterMLRoot);
set(sched,'HasSharedFilesystem',false);
%Check to see if the jobs directory exists yet
if ~exist(fullfile(CACHome, verName),'file')
    try
        mkdir(fullfile(CACHome, verName));
    catch
        s = sprintf('Unable to create DataLocation directory for storage of job files. Please manually create %s and rerun cacsched.',...
            fullfile(CACHome, verName));
    end
end
%Else, we're already all set
%Check to see if the remote jobs directory exists yet

% authenticate to myproxy server

cm.pingMarian(verName);
if ~com.mathworks.jmi.Support.useSwing()
    fprintf('You are using a terminal connection.\n');
    fprintf('This is a BETA feature to support terminal users\n');
    fprintf('If this doesn''t work, sorry!\n');
    host = input('Enter proxy server [myproxy.cac.cornell.edu]: ', ...
                 's');
    if isempty(host)
        %fprintf('setting default');
        host = 'myproxy.cac.cornell.edu';
    end
    %fprintf('host: %s\n',host);
    proxyUsername = input(sprintf('Enter  username [%s]: ', ...
                                  userName),'s');
    if isempty(proxyUsername)
        proxyUsername = userName;
    end
    %fprintf('host: %s\n',host);
    %fprintf('username: %s\n',proxyUsername);

  fprintf('Please enter password: ');
    command = '/bin/sh -c ''read -s password;echo $password;''';
    newline = true;
[status,output] = system(command);

if newline
    fprintf('\n');
end

if status == 0
    rawPassword = regexprep(output, '\n', '');
    password = iProcessBackspaces(rawPassword);
else
    message = 'Reading password failed.';
    error(message);
end
  
  cm.setMyProxySettings(proxyUsername,password,host);
  %  error('can not run\n');
end

ftp = gridFTP();
if ~ftp.isDir(RemoteDataLocation)
    %    fprintf('creating directory\n');
    ftp.createDir(RemoteDataLocation);
end



set(sched,'DataLocation',fullfile(CACHome, verName));

sched.UserData.UserName = userName;
sched.UserData.UserHome = userHome;
sched.UserData.RemotePath = remPath;
sched.UserData.ClusterHost = clusterHost;
sched.UserData.RemoteDataLocation = RemoteDataLocation;

%Set the filetransferproxy
%fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
%zipper = edu.cornell.cac.tuc.cacscheduler.globus.ftp.CollectAndZipProxy();
%fm.setProxy(zipper);

% Removes the backspaces and backspaced characters from a character
% array.
% from MATLAB 2010b
function noBackspaces = iProcessBackspaces(maybeBackspaces)

%short circuit if there are no backspaces
backspaceChar = char(8);

if all (maybeBackspaces ~= backspaceChar)
    noBackspaces = maybeBackspaces;
    return
end

noBackspaces = blanks(length(maybeBackspaces));

%next index to insert at
index = 1;
for rawIndex = 1:length(maybeBackspaces)

    %if this is a backspace
    if strcmp(backspaceChar,maybeBackspaces(rawIndex))
        %if there's something to backspace over
        if index > 1
            %move back one character
            index = index - 1;
        end
        %else ignore the backspace because there's nothing to the left
    else
        %put the character in the next spot and advance the index
        noBackspaces(index) = maybeBackspaces(rawIndex);
        index = index +1;
    end
end
if index > 1
    noBackspaces = noBackspaces(1:index-1);
else
    noBackspaces = '';
end


