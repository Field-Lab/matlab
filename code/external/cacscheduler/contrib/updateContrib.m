function bool = updateContrib()
%updateContrib

%We use a hashmap to keep track of files and their versions
import java.util.*;
contents = HashMap;

contribListing = 'http://cacscheduler.cac.cornell.edu/caccontrib/listing.cac';
try
    contribContents = urlread(contribListing);
catch ex
    error('Unable to retrieve contrib listing!');
end
%Local Contrib mirror
contribDir = fileparts(which('updateContrib'));
%Local Contrib contents
mirror = fullfile(contribDir,'mirror.cac');
if exist(mirror,'file')
    %First time we've been run or mirror.lj has been deleted, so pull
    %everything
    fid = fopen(fullfile(contribDir,'mirror.cac'),'r');
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        [fname,ver,update] = strread(tline,'%s%d%s','delimiter','\t');
        contents.put(java.lang.String(fname),ver);
    end
    fclose(fid);
end

[fname,ver,update] = strread(contribContents,'%s%d%s','delimiter','\t');
updated = 0;
notupdated = 0;
for i=1:length(fname)
    %See if the current file is up to date locally.
    if contents.containsKey(java.lang.String(fname{i}))
        %Compare versions;
        localVer = contents.get(java.lang.String(fname{i}));
        if ver(i) > localVer
            fetchFile(fname{i});
            contents.put(java.lang.String(fname{i}),ver(i));
            updated = updated+1;
        else
            fprintf('%s already at version %d, skipping...\n',fname{i},localVer);
            notupdated = notupdated+1;
        end
    else
        %No local version, add it
        fetchFile(fname{i});
        contents.put(java.lang.String(fname{i}),ver(i));
        updated = updated+1;
    end
end

%Update the contents
fid = fopen(mirror,'w');
try
    keys = contents.keySet();
    itr = keys.iterator();
    while itr.hasNext()
        fname = itr.next();
        val = contents.get(fname);
        %Write to file
        fprintf(fid,'%s\t%d\t%s\n',char(fname),val,date);
    end
    fclose(fid);
catch
    %Just try and force the closing of the file.
    fclose(fid);
end
        
fprintf('Update complete.  %d of %d files updated\n',updated, updated+notupdated);

function fetchFile(fname)
%Create urlreadpath
contribBase = 'http://cacscheduler.cac.cornell.edu/caccontrib/';
localBase = fileparts(which('updateContrib'));
try
    fileContents = urlread([contribBase fname]);
catch ex
    error('Unable to retrieve contrib File: %s!',fname);
end
fprintf('Updating %s...\n',fname);
%Write file
localFile = fullfile(localBase,fname);
%fopen 'w' will delete anything that current exists
fid = fopen(localFile,'w');
fwrite(fid,fileContents);
fclose(fid);
