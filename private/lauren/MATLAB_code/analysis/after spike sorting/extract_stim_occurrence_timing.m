function [pattern_application_times] = ...
    extract_stim_occurrence_timing(FileName, RawDataDirectory, nPatterns, iMov)

% Extracts the timing of stimulation occurrences from raw electrical
% stimulation data.  
%
% This is a modified version of the standard preprocessing function.
%
% Author: Geoff Weiner, 2013-02-27
%
% **** note: assumes that pattern numbers are 1:nPatterns (as opposed to
% e.g. [2 5 27])


%reading in the movie parameters:
ChunkData=NS_MovieData(FileName, iMov, 'fullPath', [RawDataDirectory filesep 'movie' FileName]); %stimulus information: times at which patterns are played

% interprets first 6 values in chunkdata
% see NS_DecodeMovieDataChunk for descriptions
[~, ~, ~, ~, MovieData]=NS_DecodeMovieDataChunk(ChunkData);

pattern_application_times = cell(nPatterns,2);

nEvents=length(MovieData)/3; %number of pattern applications in movie

% For each event in the movie...
for ii = 1:nEvents
    % Check to see which pattern was applied...
    for jj = 1:nPatterns
        % Then save that pattern number in the first column of
        % pattern_applcation_times and save the application time in the
        % second column.
        if MovieData(3*(ii-1)+2) == jj
            % If the pattern ID hasn't been recorded yet, put it in...
            if isempty(pattern_application_times{jj,1})
                pattern_application_times{jj,1} = jj;
            else
            end
            % Store the pattern application time.
            pattern_application_times{jj,2} = [pattern_application_times{jj,2} MovieData(3*(ii-1)+1)];
        else
        end
    end
end



















