function generate_clumped_mosaic_using_ci(path, name, varargin)
% GENERATE_CLUMPED_MOSAIC_USING_CI  Generate artificially clumped cone
% mosaics using the clumping index. There are multiple ways this function
% can generate mosaics: 'save_permuted_only' permutes the cone labels in
% dataset and saves out the resulting mosaics (exiting afterwards),
% if 'fit_permuted' is false this function will generate plots comparing
% the observed clumping index to the clumping index of n permutions 
% (where n is defined by the parameter 'iterations'). Finally, if 
% 'fit_permuted' is true, the function will permute the mosaic once 
% and then swap pairs of cones until the clumping index of the permuted 
% mosaic exceeds the clumping index of the data mosaic. If it does not
% converge within 'iterations' it will stop. Since the default value for
% iterations is orders of magnitude larger than the amount of iterations
% that should be needed, this shouldn't ever be a problem in practice.
% Only neighbor cones will be swapped if 'neighbor_clumped' is enabled.
%
% usage:  generate_clumped_mosaics_using_ci(path, name, <params>)
%
% arguments:  path - path to cones.txt for dataset to permute
%             name - name of dataset to permute
%             varargin - struct or list of optional parameters (see below)
%
% outputs:    result - plots to screen and/or text files containing
%                      permuted cone labels
%
% optional params, their default values, and what they specify:
%
% save_permuted_only      false         save permuted mosaic for each seed and break.
% fit_permuted            false         should we fit the permuted mosaic or just compare the data to the permuted clumping index.
% neighbor_clumped        false         use all pairs of cones to fit_permuted or just neighbor cones (closer than 2 median NND).
% save_text               false         save text files for each generated mosaic.      
% make_plots              false         display diagnostic plots to the screen
% start_seed              11111         random number seed.
% iterations              10000000      max number of iterations to converge.
% num_mosaics             64            repeat the analysis this many times, incrementing start_seed each time.   
% num_bins                20            number of bins to use if drps are plotted   
% bin_width               1             bin size to use if drps are plotted 
% folder                  []            name of existing subfolder to put text files in
% base_path    ../one/permuted-mosaics  path to folder
% 6/18/2009 tamachado
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

% specify list of optional parameters
p.addParamValue('save_permuted_only',false);
p.addParamValue('neighbor_clumped', false);
p.addParamValue('save_text', false);
p.addParamValue('start_seed', 11111);
p.addParamValue('radius_column', 8);
p.addParamValue('cell_column', 2);
p.addParamValue('label_column', 4);
p.addParamValue('cone_likelihood_column', 10);
p.addParamValue('iterations', 10000000);
p.addParamValue('num_mosaics',64);
p.addParamValue('fit_permuted', false);
p.addParamValue('make_plots', false);
p.addParamValue('num_bins', 20);
p.addParamValue('bin_width', 1);
p.addParamValue('on_midget', 3);
p.addParamValue('off_midget', 4);
p.addParamValue('folder', []);
p.addParamValue('base_path', '/snle/lab/Experiments/Array/Shared/one/simulations');
p.addParamValue('clumping_multiplier', 1);
%p.addParamValue('base_path', pwd);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% plot output at a given offset
if params.make_plots, f = figure; offset = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the clumping statistic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the single cone data for the current dataset
rgcPath = path;
cones = importdata(path);

% put the text file into a matrix
coneData = dlmread(path,'\t');

% load rgc connectivity data
rgcPath(strfind(path,'cones'):end) = [];
rgcs = importdata([rgcPath 'rgcs.txt']);

% find mean(robustMean(offMidgetFitRadius),robustMean(onMidgetFitRadius))
radiusOnMidget  = rgcs(rgcs(:,params.cell_column) == params.on_midget,params.radius_column);
radiusOffMidget = rgcs(rgcs(:,params.cell_column) == params.off_midget,params.radius_column);
ron  = radiusOnMidget(((radiusOnMidget > prctile(radiusOnMidget,25)) + (radiusOnMidget < prctile(radiusOnMidget,75))) > 1);
roff = radiusOffMidget(((radiusOffMidget > prctile(radiusOffMidget,25)) + (radiusOffMidget < prctile(radiusOffMidget,75))) > 1);
fd = mean([mean(ron) mean(roff)]);

% annihilate s and u cones
uInd = find(cones(:,params.label_column) == 0);
sInd = find(cones(:,params.label_column) == 3);
lmInd = setdiff(1:length(cones),union(uInd,sInd));
cones(uInd,:) = [];
cones(cones(:,params.label_column) == 3,:) = [];

% compute the median nearest neighbor distance
locations = [cones(:,2) cones(:,3)];
distances = squareform(pdist(locations));

% store nearest neighbor distances
nearest = zeros(length(cones(:,1)),1); 
nCones = length(cones(:,1));
for di = 1:nCones
    indices = setdiff(1:length(cones(:,1)),di);
    nearest(di) = min(distances(di,indices));
end

% this is the median of the nnd distribution
medianNnd = median(nearest(nearest > 0));

% get indices of each cone type
c = cell(2,1);
for cone = [1 2]
    c{cone}=find(cones(:,params.label_column) == cone);
end

conesO = cones; cO = c;
SEED = params.start_seed;
lastError = Inf;

disp(name)


%% added by gdf
% if needed compute the mean clumping statistic for a fully random cone
% mosaic
if params.clumping_multiplier ~= 1;
    num_iters = 50;
    clumping_stat_random = zeros(num_iters,1);
    for iter = 1:num_iters
        temp_cones = cones;
        nLabels = length(temp_cones(:,params.label_column));
        temp_cones(:,params.label_column) = temp_cones(randperm(nLabels),params.label_column);
        
        s = 0;
        cone_type_indices{1} = find(temp_cones(:, params.label_column) == 1);
        cone_type_indices{2} = find(temp_cones(:, params.label_column) == 2);
        for type = 1:2
            % get the distances pairwise
            locations = [temp_cones(cone_type_indices{type},2) temp_cones(cone_type_indices{type},3)];
            distances = pdist(locations);
            % weight them using a gaussian
            distances = exp(-((distances/fd).^2)/2);
            % add up the relevant distances
            s = s + sum(sum(distances));
        end        
        % accumulate clumping statistic over random mosaics
        clumping_stat_random(iter) = s;
    end
    
    mean_clumping_stat = mean(clumping_stat_random);
end
%%        


for nTimes = 1:params.num_mosaics
    
    % reset values
    c = cO;
    cones = conesO;
    histVector = zeros(params.iterations-1,1);
    
    for iteration=1:params.iterations
        % permute cones (if we're fitting a permuted mosaic, permute once only)
        if iteration > 1
            if (iteration > 2 && ~params.fit_permuted) || iteration == 2
                % set up the random number generator.
                rand('twister', SEED);
                SEED = SEED + 1;
                % permute the cone labels (column 4).
                nLabels = length(cones(:,params.label_column));
                cones(:,params.label_column) = cones(randperm(nLabels),params.label_column);
                conesP = cones;
                % don't do anything except permute once and save the results
                if params.save_permuted_only
                    coneData(lmInd,params.label_column) = cones(:,params.label_column);
                    % write out the mosaic, we've already incremented seed
                    if ~isempty(params.folder)
                       dlmwrite(sprintf([params.base_path '/%s/%s-%d.txt'],params.folder,name,SEED-1),coneData,'delimiter','\t');
                    else
                       dlmwrite(sprintf([params.base_path '/%s-%d.txt'],name,SEED-1),coneData,'delimiter','\t');
                    end
                    break;
                end
                
            elseif params.fit_permuted
                lastCones = cones;
                % once we've permuted the mosaic just swap two cones
                swap  = unidrnd(length(cones),2,1);
                % if the labels are the same, this will do nothing so continue
                if cones(swap(1),params.label_column) == cones(swap(2),params.label_column)
                    histVector(iteration-1) = histVector(iteration-2);
                    continue;
                end
                % if the cones are not closer than 2 NND, then continue
                if params.neighbor_clumped
                    locations = [cones(swap(1),2) cones(swap(1),3);...
                                 cones(swap(2),2) cones(swap(2),3)];
                    histVector(iteration-1) = histVector(iteration-2);
                    distanceNeighbor = pdist(locations);
                    if distanceNeighbor > 2*medianNnd
                        continue; 
                    end
                end
                % perform the swap
                label = cones(swap(1),params.label_column);
                cones(swap(1),params.label_column) = cones(swap(2),params.label_column);
                cones(swap(2),params.label_column) = label;
            end
        end
        % re-get cone indices
        for cone = [1 2], c{cone}=find(cones(:,params.label_column) == cone); end
        % save the permuted mosaic before we fit it
        if iteration == 2, cP = c; end
        s = 0;
        % compute distances
        for type = 1:2
            % get the distances pairwise
            locations = [cones(c{type},2) cones(c{type},3)];
            distances = pdist(locations);
            % weight them using a gaussian
            distances = exp(-((distances/fd).^2)/2);
            % add up the relevant distances
            s = s + sum(sum(distances));
        end
        %%
        if iteration == 1
            % save the clumping index of the real data
            dataValue = s;
            
            % if the clumping index is set to be different than observed in
            % the data, compute the difference between randomly permuted
            % mosaics and the observed mosaic, multiply this difference by
            % the clumping_multiplier, and use the product as the threshold
            % of generating clumped mosaics.
            if params.clumping_multiplier ~= 1
                clumping_diff = dataValue - mean_clumping_stat;
                added_clumping = clumping_diff * params.clumping_multiplier;
                dataValue = dataValue + added_clumping;
            end

        %%  
        else
            error = (dataValue - s);
            % accept the swap iff it reduces the error
            if params.fit_permuted && iteration > 3
                % if it doesn't reduce the error reject the swap
                if lastError*signE - (error)*signE <= 0
                    cones = lastCones;
                    histVector(iteration-1) = histVector(iteration-2);
                else
                    lastError = error;
                    histVector(iteration-1) = s;
                    % convergence occurs when we've swapped past the observed ci
                    if (signE == 1 && s > dataValue) || (signE == -1 && s < dataValue)
                        disp(sprintf('error = %.2f; i = %d',error,iteration))
                        d.fitValue = s;
                        if params.save_text
                            coneData(lmInd,params.label_column) = cones(:,params.label_column);
                            % write out the mosaic, we've already incremented seed
                            if ~isempty(params.folder)
                                dlmwrite(sprintf([params.base_path '/%s/%s-%d.txt'],params.folder,name,SEED-1),coneData,'delimiter','\t');
                            else
                                dlmwrite(sprintf([params.base_path '/%s-%d.txt'],name,SEED-1),coneData,'delimiter','\t');
                            end
                        end
                        break;
                    end
                end
                % plot the error
                % figure(12); plot(histVector(2:iteration-1)); hold on;
                % plot(repmat(dataValue,length(2:iteration-1),1),'r'); hold off;
            else
                % save the clumping index of the permuted data
                histVector(iteration-1) = s;
                lastError = error;
                % see if the permuted mosaic is more clumped or less clumped than the data
                if error < 0, signE = -1; else signE = 1; end;
            end
        end
    end
end


if ~params.save_permuted_only
    [drp, centers] = density_recovery_profile([cones(cO{1},2) cones(cO{1},3)], params.num_bins, params.bin_width);
    [drpFit, centersFit] = density_recovery_profile([cones(c{1},2) cones(c{1},3)], params.num_bins, params.bin_width);
    [drpP, centersP] = density_recovery_profile([conesP(c{1},2) cones(cP{1},3)], params.num_bins, params.bin_width);

    d.drp = drp; d.centers = centers;
    d.drpFit = drpFit; d.centersFit = centersFit;
    d.drpP = drpP; d.centersP = centersP;

    if params.make_plots
        if params.fit_permuted
            figure(f); subplot(2,3,1); hold on;
            plot(cones(cO{1},2),conesO(cO{1},3),'r.');% L
            plot(cones(cO{2},2),conesO(cO{2},3),'g.');% M
            title(sprintf('n = %d',length(cO{1}) + length(cO{2})));
            axis image;
            subplot(2,3,2);  bar(centers,drp); axis square; colormap bone;
            subplot(2,3,3);  plot(histVector(2:iteration-1)); hold on;
            plot(repmat(dataValue,length(2:iteration-1),1),'r');
            axis square;
            subplot(2,3,4); hold on;
            plot(cones(c{1},2),cones(c{1},3),'r.');% L
            plot(cones(c{2},2),cones(c{2},3),'g.');% M
            title(sprintf('n = %d',length(c{1}) + length(c{2})));
            axis image;
            subplot(2,3,5);  bar(centersFit,drpFit); axis square; colormap bone;
        else
            % increment the offset for this dataset
            offset = offset+0;
            % plot the cone mosaic
            figure(f); subplot(1,4,1+offset); hold on;
            plot(cones(c{1},2),cones(c{1},3),'r.');% L
            plot(cones(c{2},2),cones(c{2},3),'g.');% M
            title(sprintf('n = %d; %d(n), %d(u)',length(c{1}) + length(c{2})));
            axis image;
            % plot histogram
            figure(f); hold on;
            subplot(1,4,2+offset); hold on; hist(histVector);
            limits = ylim;
            plot([dataValue dataValue], [limits(1) limits(2)], 'r', 'LineWidth', 2);
            axis square;
            z = (dataValue - mean(histVector)) / std(histVector);
            [h, p] = ttest2(dataValue,histVector,[],'right');
            title(sprintf('v =%.2f; z = %.2f; p = %.2f',dataValue,z,p));
        end
    end
end
