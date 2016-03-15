function elecsUsed = getElecsFromPatternFile(filePath,p)
%%% Get electrodes associated with specific pattern (p) 
% input: filePath points to 'pattern_files' folder and pattern number p
% e.g., filePath = '/Volumes/Analysis/2014-11-05-05/data005/pattern_files/';
% p=1371;
% Gonzalo Mena 5/2015


prevPattern = 0; 
files = dir(filePath);
indp=strmatch(['pattern' num2str(p) '_'], char(files.name));

for i = indp
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