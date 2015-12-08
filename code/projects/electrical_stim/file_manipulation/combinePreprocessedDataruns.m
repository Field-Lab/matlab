% Code to join together two preprocessed dataruns that are 

% % newDirName = '/Volumes/Analysis/2015-09-23-8/data006-data007/';
% newDirName = '/Volumes/Analysis/2015-11-09-3/data001-data002/';
path1 = uigetdir('/Volumes/Analysis', 'Select e-stim datarun part1'); 
[exp_name,run1] = fileparts(path1);
path2 = uigetdir(exp_name, 'Select e-stim datarun part1');
[~,run2] = fileparts(path2); 

% Create new directory to save the combined analysis
% newDirName = '/Volumes/Analysis/2015-09-23-8/data006-data007/';
newDirName = fullfile(exp_name,[run1 '-' run2]);
if ~exist(newDirName,'dir')
    mkdir(newDirName);
end

cd(newDirName); 

%% For each pattern, move files into a new folder.
patternfolders = dir(fullfile(path1,'p*'));
allPatterns = 1:size(patternfolders,1)-1; 
files = dir(path1);
% Find highest movie number in the first data directory
cutoffMovieNo = size(dir(fullfile(path1,'status_files')),1)-2;
for p = 1:length(allPatterns);
    patternNo = allPatterns(p);
    
    if ~isdir(fullfile(newDirName,['p' num2str(patternNo)]))
        mkdir(fullfile(newDirName,['p' num2str(patternNo)]));
    end
    % Make sure new directory is empty
    if size(dir(fullfile(newDirName,['p' num2str(patternNo)])),1) == 2
        % Move files directly from the second folder.
        copyfile(fullfile(path2,['p' num2str(patternNo)]),fullfile(newDirName,['p' num2str(p)]));
        % Rename files by incrementing pattern number
        files = dir(fullfile(newDirName,['p' num2str(patternNo)]));
        for f = 3:size(files,1)
            oldname = files(f).name;
            newmovienum = str2double(oldname((find(oldname == '_')+2):end)) + cutoffMovieNo;
            newname = [oldname(1:find(oldname == '_')) 'm' num2str(newmovienum)];
            movefile(fullfile(newDirName,['p' num2str(patternNo)],oldname),...
                fullfile(newDirName,['p' num2str(patternNo)],newname));
        end
        
        % Check if there was overlap.
        movieNos_dir1 = findMovieNos(path1,patternNo);
        [finalAmp_dir1, ~, ~, ~, ~, ~] = getStimAmps(path1,patternNo, movieNos_dir1(end));
        movieNos_dir2 = findMovieNos(path2,patternNo);
        [firstAmp_dir2, ~, ~, ~, ~, ~] = getStimAmps(path2, patternNo, movieNos_dir2(1));
        
        if finalAmp_dir1 == firstAmp_dir2
            % copy all files except for the final movie from the first
            % directory
            disp(['p' num2str(patternNo) ' has ' num2str(movieNos_dir1(end)) ...
                ' movies; movie cutoff is ' num2str(cutoffMovieNo)]);
            files = dir(fullfile(path1,['p' num2str(patternNo)]));
            for f = 3:size(files,1)-1;
                [status,message] = movefile(fullfile(path1,['p' num2str(patternNo)],files(f).name),fullfile(newDirName,['p' num2str(p)]));
                disp(message);
            end
        else
            % copy all files from first directory
           movefile(fullfile(path1,['p' num2str(patternNo)]),fullfile(newDirName,['p' num2str(patternNo)]));
        end
    end
    disp(['copied pattern no. ' num2str(patternNo)]);
end
%%
if ~exist(fullfile(newDirName,'pattern_files'),'dir')
    mkdir(fullfile(newDirName,'pattern_files'))
end
if ~exist(fullfile(newDirName,'status_files'),'dir')
    mkdir(fullfile(newDirName,'status_files'))
end
%%
for p = 1:length(allPatterns);
    patternNo = allPatterns(p);

        % Copy pattern files directly from the second folder.
        movefile(fullfile(path2,'pattern_files',['pattern' num2str(patternNo) '_*']),...
            fullfile(newDirName,'pattern_files')); 
        % Rename files by incrementing pattern number
        patternfiles = dir(fullfile(path2,['pattern_files/pattern' num2str(patternNo) '_*']));
        for pp = 1:size(patternfiles,1)
            oldname = patternfiles(pp).name;
            newmovienum = str2double(oldname((find(oldname == '_')+2):end-4)) + cutoffMovieNo;
            movefile(fullfile(newDirName,'pattern_files',oldname),...
                fullfile(newDirName,'pattern_files',['pattern' num2str(patternNo) '_m' num2str(newmovienum) '.mat'])); 
        end
        
        % Check if there was overlap.
        movieNos_dir1 = findMovieNos(path1,patternNo);
        [finalAmp_dir1, ~, ~, ~, ~, ~] = getStimAmps(path1,patternNo, movieNos_dir1(end));
        movieNos_dir2 = findMovieNos(path2,patternNo);
        [firstAmp_dir2, ~, ~, ~, ~, ~] = getStimAmps(path2, patternNo, movieNos_dir2(1));
        
        if finalAmp_dir1 == firstAmp_dir2
            % copy all files except for the final movie from the first
            % directory
            patternfiles = dir(fullfile(path1,['pattern_files/pattern' num2str(patternNo) '*']));
            for pp = 1:size(patternfiles,1)-1
                movefile(fullfile(path1,'pattern_files',patternfiles(pp).name),...
                    fullfile(newDirName,'pattern_files'));
            end
            
        else
            movefile(fullfile(path1,'pattern_files',['pattern' num2str(patternNo) '_*']),...
                fullfile(newDirName,'pattern_files'));
        end
    
    disp(['copied pattern no. ' num2str(patternNo)]);
end
%% Status file copying.
% Copy status files directly from the second folder.
movefile(fullfile(path2,'status_files','status_m*'),...
    fullfile(newDirName,'status_files'));
% Rename files by incrementing pattern number
statusfiles = dir(fullfile(path2,'status_files/status_*'));
for st = 1:size(statusfiles,1)
    oldname = statusfiles(st).name;
    newmovienum = str2double(oldname((find(oldname == '_')+2):end-4)) + cutoffMovieNo;
    movefile(fullfile(newDirName,'status_files',oldname),...
        fullfile(newDirName,'status_files',['status_m' num2str(newmovienum) '.mat']));
end
movefile(fullfile(path1,'status_files','status_m*'),...
    fullfile(newDirName,'status_files'));

disp('copied status files  ' );

