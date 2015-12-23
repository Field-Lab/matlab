function s = qpeek(sched,job)

userData = sched.UserData;
clusterHost = userData.ClusterHost;
remoteDataLocation = userData.RemoteDataLocation;

jName = get(job,'Name');
logFile = [remoteDataLocation '\' jName '\' jName '.cacscheduler.ou'];
%Try to download file
tLogFile = tempname();
ftp = gridFTP();
if ftp.isFile(logFile)
    ftp.get(tLogFile,logFile);
else
    %warning('No log file has been created yet');
    s ='';
    return;
end
ftp.close();


%Now read the file
s = {};
fid = fopen(tLogFile);
while 1
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    s{end+1} = line;
end
fclose(fid);

s = strvcat(s);


