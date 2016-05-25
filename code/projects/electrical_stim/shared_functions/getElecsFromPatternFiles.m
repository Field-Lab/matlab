function elecsUsed = getElecsFromPatternFiles(filePath)
%%% Get electrodes associated with pattern numbers 
% input: filePath points to 'pattern_files' folder
% e.g., filePath = '/Volumes/Analysis/2012-09-24-3/data008/pattern_files/';
% LG 3/2014

prevPattern = 0; 
files = dir(filePath);
for i = 3:size(files,1)
    fileName = files(i).name;
    pattern = str2double(fileName(8:find(fileName == '_',1)-1));
    if pattern ~= prevPattern
        temp = load(fullfile(filePath, fileName));
        if isfield(temp,'pattern')
            for ii = 1:size(temp.pattern,2)
                allPatterns(pattern,ii) = temp.pattern(ii).channel;
            end
        elseif isfield(temp,'Pattern')
            for ii = 1:size(temp.Pattern,2)
                allPatterns(pattern,ii) = temp.Pattern(ii).channel;
            end
        end
    end
    prevPattern = pattern;
end
elecsUsed = sparse(allPatterns);
end
