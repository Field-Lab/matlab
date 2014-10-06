function out = cacCheckClassPath()
CACHome = fileparts(which('runTests'));
jars = javaclasspath('-static');
jarfiles = {'cacscheduler.jar','bcprov-jdk15-1.43.jar', 'bcprov-jdk16-143.jar','cog-jglobus-1.7.0.jar','commons-logging-1.1.1.jar','cryptix-asn1.jar','cryptix.jar','cryptix32.jar','cxf-2.2.10.jar','log4j-1.2.15.jar','not-yet-commons-ssl-0.3.11.jar','puretls.jar'};
foundjarfiles = {};
prevLine = 'aquaDecorations.jar';
cp = which('classpath.txt');
if isempty(CACHome)
    fprintf('could not locate cacscheduler directory\n');
else
    fprintf('cacscheduler home: %s\n',CACHome);
end
%Ensure we have a cacscheduler.jar and cxf-2.2.10.jar in there somewhere and
%that they exist.
isOK = 1;
ibsOK = 1;
wsOK = 1;
os = computer;
for i = 1:length(jars)
    % on unix/mac trailing whitespace is a problem
    % does not appear to be a problem on windows
    m=regexp(char(jars{i}),'\s$','match');
    if ~isequal(m,{})
        if ~strcmp(os,'PCWIN')     
            fprintf('ERROR: path contains whitespace:\nSTART:%s:END\n\n',jars{i});
            wsOK = 0;
            % whitespace problem, skip other checks
            continue;
        else
            jars{i} = regexprep(jars{i},'\s$','');
            %fprintf('jars{i}: %s\n',jars{i});
        end
        
    end
    if findstr('ib6https.jar',jars{i})
        fprintf('ERROR comment out line:\n$matlabroot/java/jarext/ice/ib6https.jar\n\n## $matlabroot/java/jarext/ice/ib6https.jar\n');
        ibsOK = 0;
    end
    for j=1:length(jarfiles)
        %fprintf('j: %s\n',jarfiles{j});
        m=regexp(char(jars{i}),[jarfiles{j} '$'],'match');
        if strcmp(jarfiles{j},m)
            %if strcmp('littlejohn.jar',m)
            %    if ~strcmp(prevLine,jars{(i-1)})
            %        fprintf('wrong location %s\n',jars{i-1});
            %    end
            %end
            %fprintf('found %s\n',jars{i});
            jf = java.io.File(jars{i});
            jarOK = 1;
            if ~jf.exists()
                jf_real = java.io.File(fullfile(LJHome,'lib',jarfiles{j}));
                jarOK = 0;
                fprintf('MATLAB has loaded %s,\n which says %s is located at:\n  %s\n\n but I can''t find it there.\n\n', jarfiles{j},jarfiles{j},jars{i});
                if jf_real.exists()
                    jarOK = 0;
                    fprintf('\nI found %s here:\n  %s\n\n please change classpath.txt to reflect that path.\n\n',jarfiles{j},char(jf_real));
                end
            end
            if jarOK
                foundjarfiles{end+1} = jarfiles{j};
            end
        end
    end
end
missingjars = setdiff(jarfiles,foundjarfiles);
if ~isempty(missingjars)    
    fprintf('missing the following jar files\n');
    fprintf('please add to classpath.txt:\n\n### START ###\n\n');
end
for j=1:length(missingjars)
    %fprintf('missing %s\n',missingjars{j});
    %fprintf('please add to classpath.txt:\n');
    fprintf('%s\n',fullfile(LJHome,'lib',missingjars{j}));
    isOK = 0;
end
if ~isempty(missingjars)    
    fprintf('\n### END ###\n');
end
if isOK & wsOK & ibsOK
    fprintf('Classpath entries check out, Java looks good\n\n');
    try
        CM = edu.cornell.cac.tuc.cacscheduler.globus.CertManager.getInstance();
        cmVer = CM.getClass().getPackage().getImplementationVersion();
        %fm = edu.cornell.cac.tuc.littlejohn.globus.ftp.FileMover.getInstance();
        %fmVer = fm.getClass().getPackage().getImplementationVersion();
        if strcmp(cmVer, cmVer)
            fprintf('Successfully loaded cacscheduler version %s.  Java is OK.\n', char(cmVer));
        else
            fprintf('Couldn''t load cacscheduler version info.  Java looks bad.\n');
        end
    catch exc
        warning('Failed loading cacscheduler classes, Java looks bad. error below.');
        disp(exc.message);
    end
else
    fprintf('\n\nclasspath does not appear to be valid\n');
    if ~wsOK
        fprintf('you have extra blank spaces at the end of lines in your classpath,\nplease correct.\n');
    end
        
end
%end function
end
