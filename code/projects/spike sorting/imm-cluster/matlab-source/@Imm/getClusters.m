function [cls, imm]  = getClusters(imm, e, index)
%GETCLUSTERS Fit out of sample points and return cluster struct
% [cls, imm] = getClusters(imm, e, index) Fit out of sample datapoints to 
% clusters generated for given electrode. Choose sample from posterior 
% using the index parameter, where index is an integer from 1 to max
% iterations. It is not advised to give an index s.t. index < nBurnIn,
% a warning will be displayed if this is done.
%
%   tamachado@salk.edu 1/28/08

% If we haven't analyzed the electrode, throw an error
if imm.error(e) == 0
    disp(sprintf('Error: Clustering was never done on electrode %d!', e))
    rethrow(lasterror);
end

% If index is not specified, get the MAP
if nargin < 3
    index = imm.clusters{e}.mapIndex;
end

% Check if index is in burn in period
if index <= imm.params.nBurnIn
   disp(sprintf('Warning: Index %d is in the burn in period!', index))
end

% Get spike count from electrode e
d = toStruct(imm.data, e);
nSpikes     = d.spikeCount;
projections = d.projections;
clear d;

nProjections = length(imm.params.nDimensions);

% Display warning if not all dimensions from the projection are being used
if size(projections,1) > length(imm.params.nDimensions)
    disp('Warning: Not all projection dimensions being used!')
end

% Display warning if fewer than desired dimensions 
if size(projections,1) < length(imm.params.nDimensions)
    disp(sprintf('Warning: Projection file contains fewer than %d dimensions!', length(imm.params.nDimensions)))
    nProjections = size(projections,1);
end

% Change to path that contains crp code.
path = imm.params.crpPath;
cwd  = pwd;
cd(path);

% Create the vector to store cluster assignments
assignments = zeros(nSpikes, 1);

try
    % Unless the user is doing something strange (i.e. does not want the
    % map estimate saved out), save out best k solutions--as specified by the
    % parameter file
    if index == imm.clusters{e}.bestIndices(1)
        upperBound = length(imm.clusters{e}.bestIndices);
    else
        upperBound = 1;
    end

    for b = 1:upperBound
        
        % Index of bth best index
        j = imm.clusters{e}.bestIndices(b);


        % Total number of clusters present in solution[index]
        nClusters = imm.clusters{e}.nClusters(b);
        
        lpa = zeros(nClusters,1);
        
        
        % Calculate the cluster weights and sizes
        for k=1:nClusters
           nSpk = length(find(imm.clusters{e}.classId(:,j) == k));
            
           imm.clusters{e}.weights{k} = imm.params.nPoints / nSpk;
           imm.clusters{e}.nSpikes{k} = nSpikes;
        end

        % Fit out of sample points if we're going to calculate neurons
        if ~strcmp(imm.params.modelOnly, 'true')
            % Assign each spike to a cluster
            for spike=1:nSpikes
                
                imm.clusters{e}.entropy(spike) = 0;
                
                for k=1:nClusters

                    sigma = imm.clusters{e}.covRecord{j}(:,:,k);
                    mu    = imm.clusters{e}.meanRecord{j}(:,k);

                    lp(k) =  double(fvnlp(projections(1:nProjections,spike), mu, sigma));
                    
                    lpa(k) = log(imm.clusters{e}.weights{k}) + fvnlp(projections(1:nProjections,spike), mu, sigma);
                end
                
                
                % Calculate entropy of cluster assignments
                if lp ~= 0
                    num  = exp(lp);
                    norm  = num/sum(num);
                    lnorm = log2(norm);
                    % By definition of entropy, 0 log 0 is 0
                    lnorm(find(lnorm == -Inf)) = 0;
                    norm (find(norm == 0)) = 0;
                    imm.clusters{e}.entropy(spike) = -sum((norm).*lnorm);
                end

                % Print current state to the display if we're in verbose mode
                if (imm.params.verbose == true) && (mod(spike, 10000) == 0)
                    disp(sprintf('%d of %d points assigned', spike, nSpikes));
                end
                
                
                % Save best cluster found
                lpa(lpa == 0) = NaN;
                a = find(max(lpa) == lpa);

                % Sometimes two values are equal, only save one
                assignments(spike) = a(1);

                % Reset the log probability vector
                lpa = zeros(imm.clusters{e}.kRecord(j),1);
                lpa(lpa == 0) = NaN;
            end

            % Update the imm object
            imm.clusters{e}.assignments{b} = assignments;
        end

    end

    % Return the cluster struct for the current electrode
    cls = imm.clusters{e};

    % Electrode e has completed out of sample fitting successfully
    imm.error(e) = imm.SUCCESS_FIT;

    % Restore the old working directory
    cd(sprintf(cwd));

catch

    cd(sprintf(cwd));
    disp('Error: Out of Sample Fitting Failed!');
    rethrow(lasterror)

end

