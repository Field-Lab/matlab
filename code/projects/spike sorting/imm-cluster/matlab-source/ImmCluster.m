function ImmCluster(params)
%IMMCLUSTER Function to run IMM clustering algorithm.
%   User must specify various parameters below, and provide a path
%   to a projections file. If this information is valid, the Infinite
%   Mixture Model (IMM) clustering algorithm will generate clusters.
%   The resultant neurons file will appear in a path provided below once
%   the operation has completed successfully.
%
%   All parameters are defined in the config.m function. This code should
%   not ever need to be modified.
%
%   tamachado@salk.edu 1/23/08
%

% NOTE: This code should work on Windows if the
% folder separator character is changed from '/' to '\'


% This option exists so if computation crashes we can manually
% load in the clusters.mat file and write out all existing clusters
% to a neurons file.
writeExistingNeurons = 'false';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse parameters for clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0
    % Load config.mat if a config struct is not passed in
    % Make sure to generate a config struct by calling the config function
    if exist('config.mat', 'file')
        load config
    else
        error('Run CONFIG before executing ImmCluster. Type "help config" for more information.')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the projections file
data = Data(params.prjPath, params);

% Have we already done clustering?
% We could just check for this, but this option is used to manually print
% out neurons if the code crashed unexpectedly for some reason.
if strcmp(writeExistingNeurons,'false')

    % Skip this step if clusters already exist in a mat file
    if params.clustersExist == false

        % Create a cell array to store clusters for each electrode
        cls = cell(1, params.nElectrodes);

        % Instantiate the Imm object to prepare for clustering
        imm  = Imm(params, data, params.nElectrodes);

        % Cluster all electrodes and determine which ones are bad
        for e = params.range
            try
                if params.verbose
                    tic
                end

                % Cluster electrode e
                [imm, mapIndex]  = clusterElectrode(imm, e);
                % Fit out of sample points to electrode e
                [cls{e}, imm]    = getClusters(imm, e, mapIndex);

                % Display time remaining
                if params.verbose
                    t     = (toc/60);
                    eLeft = params.range(end) - e;
                    disp(sprintf('Finished electrode %d in %0.1f minutes\n', e, t))
                    disp(sprintf('%.2f hours to go on %d more electrode(s)\n', (t*eLeft)/60, eLeft))
                end

                % Save clusters out to a file
                path = [params.oPath '/' params.itName 'clusters.mat'];
                save( path, 'cls');

            catch
                continue;
            end
        end

        % How many electrodes have been analyzed?
        [err, badElectrodes] = checkError(imm);

    else

        % Load clusters from file instead
        path = [params.oPath '/' params.itName 'clusters.mat'];
        load(path)

    end

    % Save out intermediate clustering information so out of sample fits can be
    % recomputed differently another time.

    path = [params.oPath '/' params.itName 'imm.mat'];
    save( path, 'imm');

    path = [params.oPath '/' params.itName 'error.mat'];
    save( path, 'err');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write neurons file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do not do this if we only want a model file
if ~strcmp(params.modelOnly, 'true')

    % Path to write neurons file to
    nPath = [params.oPath '/neurons/' params.itName '.neurons'];

    % Create/open neurons file at path
    nFile = NeuronsFile(nPath, data, params.minSpikes, params.nBest);

    % Add all completed electrodes to the neurons file
    for e = params.range
        try
            % Only proceed if electrode has been clustered and out of sample
            % fitting has occured
            if isfield(cls{e}, 'nClusters') && isfield(cls{e}, 'assignments')
                % Add cls to the neurons file
                nFile = add(nFile, e, cls{e});
            end
        catch
            continue;
        end
    end

    % Close the neurons file

    % Free matlab object now that we are done with it
    clear nFile;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure to use
f = 1; g = 2;

if params.savePlots
    for e = params.range

        % Check if clusters exist for electrode e
        if(isempty(cls{e}))
            continue;
        end;
        
        % Clear figures
        figure(f); clf(f);

        % Create a plot object and plot projections from electrode e
        plot = Plot(data);
        plot = toScreen(plot, f, e, cls{e});

        % Path to write projection plot e
        pPath = [params.oPath '/plots/' params.itName sprintf('-%d', e) '.pdf'];

        % Write figure out to file
        print(f, '-dpdf', pPath)
        
        % Plot entropy plot
        figure(g); clf(g);
        pointSize = 25;
        d = toStruct(data, e);
        scatter( d.projections(1,:), d.projections(2,:), pointSize, cls{e}.entropy(:))
        title(sprintf('Normalized Mean Entropy on Electrode %d is %.5g', e, -(mean(cls{e}.entropy(:)) - log(cls{e}.nClusters))))
        colorbar('EastOutside')
        
        % Write entropy plots to disk
        ePath = [params.oPath '/plots/' params.itName sprintf('-%d', e) '-entropy.png'];
        print(g, '-dpng', ePath)
        
    end
end

end

% The model file is written using the ImmModel function because all MATLAB
% processes must complete before the model can be written.


