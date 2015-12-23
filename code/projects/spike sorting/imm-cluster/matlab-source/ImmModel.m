%IMMMODEL Function to generate a model file with clusters from IMM clustering.
%   If clustering information from the Infinite Mixture of Gaussians (IMM)
%   algorithm already exists, it can be serialized to a model file using
%   this function. This function assumes there are params.nProcess number
%   of distinct clustering files that need to be parsed
%
%   All parameters are defined in the config.m function. This code should
%   not ever need to be modified.
%
%   tamachado@salk.edu 4/7/08
%
function ImmModel(params)


% This function assumes that the prefix for the .mat files containing
% clustering information are named in the form "part#clusters.mat". This is
% correct if the clustering script (matlab-auto-cluster) is used. 
% However, if it is not, the following constant needs to be changed.
% You could look at params.itName to try and figure this out dynamically. 
% This is still imperfect, unfortunately, because you only have the 
% params struct for one of the n processes.
filePrefix = 'part';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read clustering information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first index corresponds to the MAP.
% This is here because eventually we'll probably want to add
% the ability to save out multiple indices (samples from posterior)
% at once. This is hard due to memory constraints, and because cluster
% information for other indices is not currently saved out for indices
% besides the first one.
index = 1;

% Store all clusters in this cell array
c = cell(params.nElectrodes, 1);
    

% Load each file
for k = 1:params.nProcess
    
    empty  = ones(params.nElectrodes, 1);
    temp   = true;
    first  = -1;
    
    fPath = [params.oPath '/' filePrefix sprintf('%d',k) 'clusters.mat'];
    load(fPath);
    disp(sprintf('Opening part %d...',k))
    
    for e = 1:params.nElectrodes
        empty(e) = isempty(cls{e});
        
        % Last populated electrode in this piece 
        if empty(e) == false
           last = e;
        end
        
        % First populated electrode in this piece
        if first == -1
            if temp == true && empty(e) == false
                first = e;
            end
        end
        temp = empty(e);  
    end
    
    for e = first:last
        % Check if this is a bad electrode
        if ~isempty(cls{e})
            % Get the index of the nth best index
            if isfield(cls{e}, 'bestIndices')
                j = cls{e}.bestIndices(index);
            else
                j = cls{e}.mapIndex;
            end

            nClusters = cls{e}.nClusters(index);

            cls{e} = rmfield(cls{e}, 'weights');
            cls{e} = rmfield(cls{e}, 'nSpikes');

            % Calculate weights and nSpikes per cluster
            if ~isfield(cls{e}, 'weights')
                for n=1:nClusters
                    nSpikes = length(find(cls{e}.classId(:,j) == n));

                    cls{e}.weights(n) = params.nPoints / nSpikes;
                    cls{e}.nSpikes(n) = nSpikes;
                end
            end

            % Get rid of unneeded fields
            if isfield(cls{e}, 'assignments')
               cls{e} = rmfield(cls{e}, 'assignments');
            end
            
            cls{e} = rmfield(cls{e}, 'classId');
            cls{e} = rmfield(cls{e}, 'lpRecord');

            c{e} = cls{e};

        else
            disp(sprintf('Electrode %d is bad. Skipping...', e))
        end
    end

    % Clear the temporary file
    clear cls;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write model file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to write model file to
mPath = [params.oPath '/model/' 'output.model'];

% Create/open model file at path
[mFile, err] = ModelFile(mPath, params.minSpikes, params.modelPath);

% Check for error
if err
    return;
end


% Add all completed electrodes to the model file
for e = 1:params.nElectrodes
    try
        % Only proceed if electrode has been clustered
        if isfield(c{e}, 'nClusters')
            % Add cls to the neurons file
            mFile = add(mFile, e, c{e});
        end
    catch
        continue;
    end
end

% Close the model file
add(mFile);

% Free matlab object now that we are done with it
clear mFile;
