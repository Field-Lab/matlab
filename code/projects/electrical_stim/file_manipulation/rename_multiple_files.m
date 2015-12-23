
% once you are in the directory with the files to be renamed
fileNames = dir;
for i = 1:length(fileNames)
    oldName = fileNames(i).name;
    
    nameChanged = 0;
    for j = 100:300
        if ~isempty(strfind(oldName, ['m' num2str(j)]))
            newName = strrep(oldName, ['m' num2str(j)], ['m' num2str(j+77)]);
            nameChanged = 1;
        end
    end
    
    for j = 10:99
        if (~nameChanged) && ~isempty(strfind(oldName, ['m' num2str(j)]))
            newName = strrep(oldName, ['m' num2str(j)], ['m' num2str(j+77)]);
            nameChanged = 1;
        end
    end
    
    for j = 1:9
        if (~nameChanged) && ~isempty(strfind(oldName, ['m' num2str(j)]))
            newName = strrep(oldName, ['m' num2str(j)], ['m' num2str(j+77)]);
            nameChanged = 1;
        end
    end
    
    %newName = strrep(oldName, 'B2-R', 'A2-R'); % replaces any instances of the first string with
                                               % the second string (if there are no instances, does nothing)
    %if ~strcmpi(oldName, newName)
    if nameChanged
        movefile(oldName, ['../data015_new/' newName])
    else
        warndlg(['name = ' oldName ' wasn''t changed for some reason'])
    end
end