function bool = clearContrib()
%clearContrib clears out all contrib functions.  This is useful to rescue
%from an update that has failed and left contrib in an unknown state.

contribDir = fileparts(which('updateContrib'));
%Local Contrib contents
mirror = fullfile(contribDir,'mirror.lj');
if exist(mirror,'file')
    fid = fopen(fullfile(contribDir,'mirror.lj'),'r');
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        [fname,ver,update] = strread(tline,'%s%d%s','delimiter','\t');
        %Try to delete the file, but don't worry if it doesn't exist
        if exist(fname{1},'file')
            delete(fullfile(contribDir, fname{1}));
        end 
    end
    fclose(fid);
end

%Delete the mirror
try
    delete(mirror);
catch
    fprintf('Unable to delete %s, delete manually if it exists.\n',mirror);
end