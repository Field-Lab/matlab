
%cd to desired directory

patternNos = [1:8 10:24 26:56 58:64];

%patternNos = 1:124;

nPatterns = length(patternNos);
for i = 1:nPatterns
    mkdir(['p' num2str(patternNos(i))])
end

files = dir;
for i = 1:length(files)
    disp(i)
    for j = 1:nPatterns
        if ~isempty(strfind(files(i).name, ['p' num2str(patternNos(j)) '_m']))
            copyfile(files(i).name, ['p' num2str(patternNos(j)) filesep files(i).name])
            delete(files(i).name)
        end
    end
end