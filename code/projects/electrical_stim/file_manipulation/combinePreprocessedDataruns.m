% Code to join together two preprocessed dataruns that are 

% Select the two paths to be joined
% path1 = '/Volumes/Analysis/2015-09-23-8/data006';
% path2 = '/Volumes/Analysis/2015-09-23-8/data007';
% % *** Need to update the code to handle 2 electrode stimulation ***

path1 = '/Volumes/Analysis/2015-10-06-3/data001';
path2 = '/Volumes/Analysis/2015-10-06-3/data002';

% Create new directory to save the combined analysis
% newDirName = '/Volumes/Analysis/2015-09-23-8/data006-data007/';
newDirName = '/Volumes/Analysis/2015-10-06-3/data001-data002/';
if ~exist(newDirName,'dir')
    mkdir(newDirName);
end

cd(newDirName); 

%%
allPatterns = 1:512; % Update this to automatically tell what's in here
% Find highest movie number in the first data directory
cutoffMovieNo = 177; 
for p = 1:length(allPatterns);
    patternNo = allPatterns(p); 
    movieNos_dir1 = findMovieNos(path1,patternNo);
    [finalAmp_dir1, ~, ~, ~, ~, ~] = getStimAmps(path1,patternNo, movieNos_dir1(end));
  
    movieNos_dir2 = findMovieNos(path2,patternNo);
    [firstAmp_dir2, ~, ~, ~, ~, ~] = getStimAmps(path2, patternNo, movieNos_dir2(1));
          
    % Check if there was overlap. 
    if finalAmp_dir1 == firstAmp_dir2
        
%         movieNos1(end); %movie 177
        % Save movies2 in new data directory, increment movie file by
        % cutoffMovieNo - 1
        
    else
        copyfile(fullfile([path1,['p' num2str(p)],)
    end
end
%%
for i = 1:length(files)
    disp(i)
    for j = 1:nPatterns
        if ~isempty(strfind(files(i).name, ['p' num2str(patternNos(j)) '_m']))
            copyfile(files(i).name, ['p' num2str(patternNos(j)) filesep files(i).name])
            delete(files(i).name)
        end
    end
end