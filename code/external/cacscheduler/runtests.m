function runtests(sched,level)
%runTests(sched,level) - test function to test cacscheduler Install.
%   level == 1 - run basic tests that don't have external communication
%   level == 2 - run component tests that require external resources
%   level == 3 - run integration tests that simulate actual use
%   level == 0 - run all levels
%
% EXAMPLES:
%  cd CACSCHEDULERHOME;
%  cacsched;
%  runtests(sched,3);

%We know runTests must be on the path, so we use that.
CACHome = fileparts(which('runtests'));
addpath(fullfile(CACHome, 'examples'));

if level == 1 || level == 0
    fails = 0;
    fprintf('**Running basic install tests to ensure install is OK***\n');
    %Test to see if the userhome from cacsched looks normal
    runUserHomeTest(sched);
    %Test that cacscheduler jars are in the javaclasspath correctly.
    runJavaClassPathTest(sched);
    %test that DataLocation exists and is writable.
    runDataLocationTest(sched);
    fprintf('**Basic install tests complete.  Resolve any errors before continuing.******\n');
    fprintf('\n');
end

if level == 2 || level==0
    fprintf('**Running functional tests of cacscheduler components.***\n');
    %Test gridFTP is working
    rungridFTPTest(sched);
    %See if we can submit a simple 1 core job
    runJobSubmit(sched);
    fprintf('***Functional component tests passed***\n');
end

if level == 3 || level ==0
    fprintf('Running integration tests that simulate actual use of TUC by a user.\n');
    getFreeCores();

    distFails = runDistributedJobTest(sched);
    parFails = runParallelJobTest(sched);
    parCanFails = runParallelCancelJobTest(sched);
    distCanFails = runDistributedCancelJobTest(sched);
    tempDirFails = runTempDirTest(sched);
    
    %Report
    if distFails > 0
        fprintf('Distributed Job test failed, see above for reasons and contact CAC for help.\n');
    elseif parFails > 0
        fprintf('Parallel Job test failed, see above for reasons and contact CAC for help.\n');
    elseif parCanFails > 0
        fprintf('Parallel cancel test failed, see above for reasons and contact CAC for help.\n');
    elseif distCanFails > 0
        fprintf('Distributed cancel test failed, see above for reasons and contact CAC for help.\n');
    elseif tempDirFails >0
        fprintf('Temporary directory test failed, see above for reasons and contact CAC for help.\n');
    else
        fprintf('Tests passed, things appear to be working for you!\n');
    end
    rmpath(fullfile(CACHome, 'examples'));
end
%End function
end


    function runUserHomeTest(sched)
        ud1 = get(sched,'SubmitFcn');
        ud2 = get(sched,'ParallelSubmitFcn');
        %Check that both are the same
        if ~strcmp(ud1{3},ud2{3})
            fprintf('A different userhome has been specified for distributed and parallel jobs.  This is likely a bad configuration.\n');
            fprintf('SubmitFcn Home: %s\n', ud1{3});
            fprintf('SubmitFcn Home: %s\n', ud2{3});
            fprintf('Edit cacsched to correct this.\n');
        else
            %Make sure that the username is valid.
            %'\\storage01\matlab\aaa11
            %'\\storage01\matlab\aaa11\maybeotherstuff'
            userHome = ud1{3};
            n = regexp(userHome,'\\storage01\\matlab\\(?<userName>\w+)','names');
            %Shouldn't be USERNAME
            if strcmpi(n.userName,'username')
                warning('cacsched reports that the username is %s, you should edit cacsched.m to replace %s with your cac username.\n'...
                    , n.userName, n.userName);
            elseif ~(isnan(str2double(n.userName(1))) && ~isnan(str2double(n.userName(end))))
                warning('Username from cacsched (%s) doesn''t look valid, this should be your cac username.\n', n.userName);
                fprintf('If that isn''t your cac username (the one you use to log on to MyProxy, correct in cacsched\n');
            end
        end
        %end function
    end

    function out = runJavaClassPathTest(sched)
        jars = javaclasspath('-static');
        cp = which('classpath.txt');
        CACHome = fileparts(which('runTests'));
        %Ensure we have a cacscheduler.jar and cxf-2.2.7.jar in there somewhere and
        %that they exist.
        cxfOK = 0;
        ljOK = 0;
        for i = 1:length(jars)
            if ~isempty(findstr('cxf-2.2.7.jar',jars{i}))
                %See if the path is real, this is easiest done with a java File.
                jf =  java.io.File(jars{i});
                if ~jf.exists()
                    %See if I can find the file at CACHOME/lib/cxf-2.2.7.jar
                    jf_real = java.io.File(fullfile(CACHome,'lib','cxf-2.2.7.jar'));
                    warning('Found cxf-2.2.7.jar in the javaclasspath, but it doesn''t seem to exist on your file system.\n');
                    fprintf('MATLAB has loaded %s, which says cxf-2.2.7 is located %s but I can''t find it there.\n', cp);
                    if jf_real.exists()
                        fprintf('I found cxf-2.2.7.jar here:%s, please change classpath.txt to reflect that path.',char(jf_real));
                    end
                else
                    %fprintf('Found %s, classpath looks OK\n', char(jf.getAbsolutePath()));
                    cxfOK = 1;
                end
            end
            %Also check cacscheduler.jar
            if ~isempty(findstr('cacscheduler.jar',jars{i}))
                jf = java.io.File(jars{i});
                if ~jf.exists()
                    jf_real = java.io.File(fullfile(CACHome,'lib','cacscheduler.jar'));
                    warning('Found littejohn.jar in the javaclasspath, but it doesn''t seem to exist on your file system.\n');
                    fprintf('MATLAB has loaded %s, which says litteljohn.jar is located %s but I can''t find it there.\n', cp);
                    if jf_real.exists()
                        fprintf('I found litteljohn.jar here:%s, please change classpath.txt to reflect that path.',char(jf_real));
                    end
                else
                    %fprintf('Found %s, classpath looks OK\n', char(jf.getAbsolutePath()));
                    ljOK = 1;
                end
            end
        end
        if cxfOK && ljOK
            fprintf('Classpath entries check out, Java looks good, trying to load cacscheduler classes\n');
            try
                CM = edu.cornell.cac.tuc.cacscheduler.globus.CertManager.getInstance();
                cmVer = CM.getClass().getPackage().getImplementationVersion();
                fm = edu.cornell.cac.tuc.cacscheduler.globus.ftp.FileMover.getInstance();
                fmVer = fm.getClass().getPackage().getImplementationVersion();
                if strcmp(cmVer, fmVer)
                    fprintf('Successfully loaded cacscheduler version %s.  Java is OK.\n', char(fmVer));
                else
                    fprintf('Couldn''t load cacscheduler version info.  Java looks bad.\n');
                end
            catch exc
                warning('Failed loading cacscheduler classes, Java looks bad. error below.');
                disp(exc.message);
            end
        end
        %end function
    end

    function     runDataLocationTest(sched)
        %Make sure that DataLocation is OK.
        pth = get(sched,'DataLocation');
        %Make sure it exists
        jf = java.io.File(pth);
        if jf.exists()
            %ensure writable
            try
                fn = fullfile(pth,'test.out');
                fid = fopen(fn,'w');
                fwrite(fid,[1,2,3]);
                fclose(fid);
                fprintf('DataLocation is writable and ready to go\n');
            catch exc
                warning('Unable to write to DataLocation:' + pth);
                fprintf('Open cacsched.m and change the ''DataLocation'' property to a writeable directory\n');
            end
        else
            warning('DataLocation '+pth+'doesn''t exist');
            fprintf('Please create %s or edit cacsched.m to point to an existing directory\n', pth);
        end
        %end function
    end

    function rungridFTPTest(sched)
        %Let's see if gridFTP is working.
        CACHome = fileparts(which('runTests'));
        ftp = gridFTP();
        ftp.list('',1);
        localPath = fullfile(CACHome,'runtests.m');
        %Check to see if the remote file already exists
        if ftp.isFile('runtests.m')
            ftp.delete('runtests.m');
        end
        ftp.put(localPath,'runtests.m');
        files = ftp.list('*.m',1);
        foundruntests = 0;
        for i = 1:length(files)
            if strcmp(files{i}.name,'runtests.m')
                foundruntests = 1;
                break;
            end
        end
        if foundruntests == 1
            fprintf('gridFTP appears to be working OK!\n');
        else
            fprintf('gridFTP failed!\n');
        end
        ftp.close();
    end

    function runJobSubmit(sched)
        fprintf('Attempting to submit a simple test job to TUC.\n');
        %Do a simple job submission to see if we can submit a non-matlab job
        try
            hpc = edu.cornell.cac.tuc.matlab.JSDLMediator(1);
            A(1) = java.lang.String('/all');
            %jobinfo = hpc.jobStart(matlabexe, args, stdoutFileLocation, stderrFileLocation,env,1);
            %props = java.util.HashMap();
            % grab cacprops to be consistent
            cacprops = getClusterInfoProperties();
            jobinfo = hpc.jobStart('ipconfig', A, '', '',A,cacprops,1);
            try
                str2double(jobinfo.getJobID());
                fprintf('Service communication looks OK\n.');
            catch exc
                warning('Encountered an error communicating with the scheduler.\n');
                disp(exc.message);
            end
        catch exc
            warning('Encountered an error communicating with the scheduler.\n');
            if ismac() && ~isempty(strfind(exc.message,'Keystore was tampered with, or password was incorrect'))
            	%this is the special case of the broken cacerts from the Apple javaupdate
            	fprintf('You''ve likely just bumped into a Mac bug related to a 12/09 Java update.\n');
            	fprintf('You will need sudo powers on this machine to modify the default password on the cacerts file.\n');
            	fprintf('Open a termina and enter the following line to correct this (no KNOWN side effects)\n.');
            	fprintf('sudo keytool --storepasswd --new changeit --keystore /System/Library/Frameworks/JavaVM.framework/Resources/Deploy.bundle/Contents/Home/lib/security/cacerts --storepasswd changeme\n');
            	fprintf('Once complete, rerun runtests(sched,2)\n');
            else
            	disp(exc.message);
            end
        end
    end

    function getFreeCores()
        %Get tucUsage, which may require a updateContrib, which may require a
        %path addition
        CACHome = fileparts(which('runTests'));
        if ~exist('updateContrib.m','file')
            fprintf('Contrib folder not on path, temporarily adding it.\n');
            addpath(fullfile(CACHome,'contrib'));
        end
        %OK, we've done our best to get updateContrib available.  Try getStatus
        if ~exist('getTUCStatus.m','file')
            %is it because contrib still isn't on the path?
            if ~exist('updateContrib.m','file')
                fprintf('Odd, can''t seem to find updateContrib after path add, might not be able to give you TUC current usage.\n');
                rempath(fullfile(CACHome,'contrib'));
                freeCores = -1;
                return;
            else
                %it exists, call it
                updateContrib();
            end
        end
        %OK, contrib is on the path and we've updated contrib.
        if ~exist('getTUCStatus.m','file')
            fprintf('Can''t seem to retrieve TUC status.  Hopefully we''re not busy!\n');
        else
            getTUCStatus();
        end
    end
    
    function fails = runDistributedJobTest(sched)
    fails = 0;
    fprintf('Running distributed job test...\n');
    [j,result] = cacsubmit(sched);
    %Make sure we got back two matrices
    if (length(result) ~= 2)
        warning('cacscheduler:DistributedJobTest', 'Distributed Job test should have returned 2 results, but returned %d', length(result));
        fails = fails+1;
    end
    %Make sure both matrices are 3x3
    for i =1:2
        [r,c] = size(result{i});
        if (r ~= 3 && c ~= 3)
            warning('cacscheduler:DistributedJobTest','Distributed Job test should have returned 3x3 matrices, instead matrix %d is %dx%d',i,r,c);
            fails = fails+1;
        end
    end
    %Make sure matrices are non-zero
    for i=1:2
        if (sum(sum(result{i})) == 0)
            warning('LitleJohn:DistributedJobTest','Distributed Job test should have returned non-zeros matrices');
            fails = fails+1;
        end
    end
    end
    function fails = runDistributedCancelJobTest(sched)
    % test that we can cancel a job
    fails = 0;
    fprintf('Running distributed cancel job test...\n');
    j = caccancelsubmit(sched);
    if ~strcmp(j.state,'finished')
        warning('LittleJoh:DistributedJobCancelTest','Distributed Cancel test should have returned "finished" state');
        fails = fails+1;
    end
    if ~findstr('cancelled',j.tasks(1).ErrorMessage)
        warning('LittleJoh:DistributedJobCancelTest','Distributed Job cancel test should have returned reason for finished state');
        fails = fails+1;
    end
    end
    
    function fails = runParallelJobTest(sched)
    fprintf('Running parallel job test...\n');
    fails = 0;
    [j,result] = cacparsubmit(sched);
    if (length(result) ~= 4)
        warning('cacscheduler:ParallelJobTest','Parallel Job test should have returned 4 results, but returned %d', length(result));
        fails = fails+1;
    end
    %make sure that all are equals
    for i=1:length(result);
        if result{i} ~= 136
            warning('cacscheduler:ParallelJobTest','Parallel Job should have returned all values being equal to 136, but instead I process %d returned %4.2f',i,result{i});
            fails = fails+1;
        end
    end
    end
    
    function fails = runParallelCancelJobTest(sched)
    % test that we can cancel a job
    fails = 0;
    fprintf('Running parallel cancel job test...\n');
    j = cacparcancelsubmit(sched);
    if ~strcmp(j.state,'finished')
        warning('LittleJoh:ParallelJobCancelTest','Parallel Job cancel test should have returned "finished" state');
        fails = fails+1;
    end
    if ~findstr('cancelled',j.tasks(1).ErrorMessage)
        warning('LittleJoh:ParallelJobCancelTest','Parallel Job cancel test should have returned reason for finished state');
        fails = fails+1;
    end
    end
    
    function fails = runDirJobTest(sched)
    fails = 0;
    fprintf('Running unique directory distributed test....\n');
    [j,result] = cacextraDirsubmit(sched);
    if (length(result) ~= 4)
        warning('cacscheduler:DistributedJobTest','Distributed Job with unique dir should have returned 4 results but returned %d', length(result));
        fails = fails + 1;
    end
    
    for i =1:length(result)
        [r,c] = size(result{i});
        if (r ~= 25 && c ~= 75)
            warning('cacscheduler:DistributedJobTest','Distributed Job test should have returned 25x75 matrices, instead matrix %d is %dx%d',i,r,c);
            fails = fails+1;
        end
    end
    end
    
    function fails = runTempDirTest(sched)
    fails = 0;
    fprintf('Running temp directory test...\n');
    status = cacsubmit_bigdata(sched);
    if ~strcmp(status,'success')
        fails = 1;
    end
    end
