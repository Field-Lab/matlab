function cacMakeClassPath(override)
startline = 'maci64=$matlabroot/java/jarext/aquaDecorations.jar';
skipline = '$matlabroot/java/jarext/ice/ib6https.jar';
jarfiles = {'cacscheduler.jar','bcprov-jdk15-1.43.jar', 'bcprov-jdk16-143.jar','cog-jglobus-1.7.0.jar','commons-logging-1.1.1.jar','cryptix-asn1.jar','cryptix.jar','cryptix32.jar','cxf-2.2.10.jar','log4j-1.2.15.jar','not-yet-commons-ssl-0.3.11.jar','puretls.jar'};
LJHome = fileparts(which('runtests'));
LJLib = fullfile(LJHome,'lib');
classpath = which('classpath.txt');
fprintf('found classpath: %s\n',classpath);
mroot = matlabroot;
fs = findstr(classpath,mroot);
if nargin==1 | fs == 1
    fprintf('Creating a custom classpath\n');
    fid = fopen(classpath);
    % assume that current directory is userpath
    tmpfile = fullfile(pwd,'classpath.txt.cac');
    realfile = fullfile(pwd,'classpath.txt');

    fout = fopen(tmpfile,'w');
    while ~feof(fid)
        tline = fgets(fid);
        
        tline_l = lower(tline);
        % skip lines that already include CACSCHEDULER lib
        if findstr(LJLib,tline_l) == 1
            fprintf('skipping: %s',tline);
            continue
        end
        % scan file jarfiles to skip
        skipLine = 0;
        for i = 1:length(jarfiles)
            if findstr(tline,jarfiles{i})
                fprintf('skipping: %s',tline);
                skipLine = 1;
                continue;
            end
        end
        if skipLine
            continue;
        end
        if findstr(tline,'littlejohn.jar')
            fprintf('skipping: %s',tline);
            continue;
        end
        if findstr(tline,'cxf-2.2.7.jar')
            fprintf('skipping: %s',tline);
            continue;
        end
        if findstr(tline,skipline) == 1
            fprintf(fout,'## %s',tline);
            fprintf('commenting out: %s',tline);
            continue
        end
        fwrite(fout,tline);
        % look for startline and add custom jars after it
        if findstr(tline,startline) == 1           
            for i = 1:length(jarfiles)
                s = fullfile(LJLib,jarfiles{i});
                fprintf(fout,'%s\n',s);
                fprintf('adding %s\n',s);
            end
        end
    end

    fclose(fid);
    fclose(fout);
    % rename any existing file
    cf = java.io.File(realfile);
    if cf.exists()
        movefile(realfile,fullfile(pwd,'classpath.txt.cac.orig'),'f');
    end

    % all done, now rename file
    movefile(tmpfile,realfile,'f');
    
else
    fprintf('looks like you already have a custom classpath\n\nTo overwrite the existing classpath \nrun cacMakeClassPath(1)\n\n');
end
