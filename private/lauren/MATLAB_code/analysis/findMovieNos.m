function movieNos = findMovieNos(pathToAnalysisData,patternNo)
% get movie numbers associated with a particular patterned stimulation
%  inputs: pathToAnalysisData        string to Analysis directory
%          patternNo                 number corresponding to pattnern
% outputs: movieNos                  vector containing movie numbers 
% % LG 8/2014
movieNos = [];
patternNoString = ['p' num2str(patternNo)];
files = dir([pathToAnalysisData filesep patternNoString]);

for i = 1:length(files)
    if strfind(files(i).name, patternNoString) == 1
        mIndices = strfind(files(i).name, 'm');
        movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
    end
end
movieNos = sort(movieNos);
end