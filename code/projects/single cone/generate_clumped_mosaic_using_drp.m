function [drp_data,cone_data,error_vector] = generate_clumped_mosaic_using_drp(path, name, varargin)
% GENERATE_CLUMPED_MOSAIC_USING_DRP  Generate artificially clumped cone
% mosaics by fitting to four density recovery profiles: DRP(L,L), DRP(M,M),
% DRP(L,M), and DRP(M,L). The algorithm permutes the cone labels once
% and then on each subsequent iteration it accepts a swap of two cone
% labels if the swap decreases the summed squared error between the
% observed drps and the drps of the permuted mosaic. When the error falls
% below 'threshold,' the algorithm stops.
%
% usage:  [drp_data, cone_data] = generate_clumped_mosaics_using_drp(path, name, <params>)
%
% arguments:  path - path to cones.txt for dataset to permute
%             name - name of dataset to permute
%             varargin - struct or list of optional parameters (see below)
%
% outputs:    result    - plots to screen and/or text files containing
%                         permuted cone labels
%             drp_data  - struct containing information about permuted mosaic
%             cone_data - struct containing information about data mosaic
%
% optional params, their default values, and what they specify:
%
% save_text               false         save text files for each generated mosaic.      
% threshold               2E-6          convergence when error is below threshold.
% start_seed              11111         random number seed.
% iterations              100000        max number of iterations to converge.
% num_mosaics             64            repeat the analysis this many times, incrementing start_seed each time.   
% num_bins                20            number of bins to use in each drp. 
% bin_width               1             bin size to use in each drp.
% folder                  []            name of existing subfolder to put text files in.
% base_path   .../one/permuted-mosaics  path to folder.
% 6/18/2009 tamachado
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up optional arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;

% specify list of optional parameters
p.addParamValue('label_column', 4);
p.addParamValue('threshold', 2E-6);
p.addParamValue('iterations', 100000);
p.addParamValue('num_mosaics',64);
p.addParamValue('save_text', false);
p.addParamValue('start_seed', 11111);
p.addParamValue('cone_likelihood_column', 10);
p.addParamValue('make_plots', false);
p.addParamValue('num_bins', 20);
p.addParamValue('bin_width', 1);
p.addParamValue('folder', []);
p.addParamValue('base_path', '/snle/lab/Experiments/Array/Shared/one/simulations');
%p.addParamValue('base_path', pwd);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare to fit to the drp
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the single cone data for the current dataset
cones = importdata(path);

% put the text file into a matrix
coneData = dlmread(path,'\t');

% get cone type indices
uInd  = find(cones(:,params.label_column) == 0);
sInd  = find(cones(:,params.label_column) == 3);
lmInd = setdiff(1:length(cones),union(uInd,sInd));
cones(uInd,:) = [];
cones(cones(:,params.label_column) == 3,:) = [];

% get cone indices
c = cell(2,1);
for cone = [1 2], c{cone}=find(cones(:,params.label_column) == cone); end

% set up values before iterating
conesO = cones; cO = c;
SEED = params.start_seed;
lastError = Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin permuting
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('%s (fitting to drp)',name))

for  nTimes = 1:params.num_mosaics
    
    % reset values
    c = cO;
    cones = conesO;
    error_vector = zeros(params.iterations-1,1);
    
    for iteration=1:params.iterations
        % permute cones once only (on iteration 2; first iteration is real data)
        if iteration > 1
            if iteration == 2
                % set up the random number generator.
                rand('twister', SEED);
                SEED = SEED + 1;
                % permute the cone labels (column 4).
                nLabels = length(cones(:,params.label_column));
                cones(:,params.label_column) = cones(randperm(nLabels),params.label_column);
            else
                lastCones = cones;
                % once we've permuted the mosaic just swap two cones
                swap  = unidrnd(length(cones),2,1);
                label = cones(swap(1),params.label_column);
                cones(swap(1),params.label_column) = cones(swap(2),params.label_column);
                cones(swap(2),params.label_column) = label;
                % if both labels are the same then don't do anything
                if length(unique(cones(swap,params.label_column))) == 1
                    error_vector(iteration-1) = error_vector(iteration-2);
                    continue;
                end
            end
        end
        
        % get updated cone indices
        for cone = [1 2], c{cone}=find(cones(:,params.label_column) == cone); end
        
        if iteration < 3
            % generate the drps on the first iteration
            [drpL,centers] = density_recovery_profile([cones(c{1},2) cones(c{1},3)], params.num_bins, params.bin_width);
            [drpM] = density_recovery_profile([cones(c{2},2) cones(c{2},3)], params.num_bins, params.bin_width);
            [drpA] = density_recovery_profile([cones(c{1},2) cones(c{1},3)], params.num_bins, params.bin_width,'reference_centers',[cones(c{2},2) cones(c{2},3)]);
            [drpB] = density_recovery_profile([cones(c{2},2) cones(c{2},3)], params.num_bins, params.bin_width,'reference_centers',[cones(c{1},2) cones(c{1},3)]);

            if iteration == 1
                % save the clumping index of the real data
                dataValue = [drpL drpM drpA drpB];
            elseif iteration == 2
                % save out error information for permuted mosaic
                drp = [drpL drpM drpA drpB];
                error_vector(iteration-1) = sum((dataValue - drp).^2);
                lastError = error_vector(iteration-1);
                % compute all distances pairwise
                allDistances = squareform(pdist([cones(:,2) cones(:,3)]));
            end
        else
            % update each drp quickly for each proposed swap
            L = 1; M = 2;
            drpL = update_drp_fast(c,allDistances,params.num_bins,params.bin_width,L,L);
            drpM = update_drp_fast(c,allDistances,params.num_bins,params.bin_width,M,M);
            drpA = update_drp_fast(c,allDistances,params.num_bins,params.bin_width,L,M);
            drpB = update_drp_fast(c,allDistances,params.num_bins,params.bin_width,M,L);
            
            % update the error
            drp = [drpL drpM drpA drpB];
            error = sum((dataValue - drp).^2);
            
            if lastError - error <= 0
                % if the swap doesn't reduce the error reject it
                cones = lastCones;
                error_vector(iteration-1) = error_vector(iteration-2);
            else
                % otherwise accept the swap
                lastError = error;
                error_vector(iteration-1) = error;
                currDrp = drp;
                if error < params.threshold
                    disp(sprintf('error = %.2f; i = %d',error,iteration))
                    d.fitValue = drp;
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
            if params.make_plots
                % plot the error
                figure(12); plot(error_vector(2:iteration-1)); hold on;
                plot(repmat(params.threshold,length(2:iteration-1),1),'r'); hold off;
            end
        end
    end
    
    if params.make_plots
        figure(88); clf;
        
        % data cone mosaic
        subplot(3,3,1); hold on;
        plot(cones(cO{1},2),conesO(cO{1},3),'r.');% L
        plot(cones(cO{2},2),conesO(cO{2},3),'g.');% M
        title(sprintf('n = %d',length(c{1}) + length(c{2})));
        axis image;
        
        % data drp
        subplot(3,3,2); bar(centers,dataValue(1:params.num_bins)); axis square; colormap bone;
        title('L drp of data mosaic (drp_d)')
        
        % error vs. iteration
        subplot(3,3,3); plot(error_vector(2:iteration-1)); hold on;
        plot(repmat(params.threshold,length(2:iteration-1),1),'r'); hold off;
        xlabel('iteration'); ylabel('error');
        title(sprintf('%d params.iterations',iteration-1));
        axis square;
        
        % fitted cone mosaic
        subplot(3,3,4); hold on;
        plot(cones(c{1},2),cones(c{1},3),'r.');% L
        plot(cones(c{2},2),cones(c{2},3),'g.');% M
        title(sprintf('n = %d; seed = %d',length(c{1}) + length(c{2}),SEED-1));
        axis image;
        
        % fitted drp
        subplot(3,3,5); 
        bar(centers,currDrp(1:params.num_bins)); axis square; colormap bone;
        title('L drp of fitted mosaic (drp_f)')
        
        % difference of drps
        subplot(3,3,6);
        plot(dataValue-currDrp, 'k'); axis square;
        title('drp_d - drp_f')
        
        % L to M drp of data mosaic
        subplot(3,3,7);
        bar(centers,dataValue(3*params.num_bins+1:4*params.num_bins));
        title('L to M drp of data mosaic');
        axis square; colormap bone;
        
        % L to M drp of fitted mosaic
        subplot(3,3,8);
        bar(centers,currDrp(3*params.num_bins+1:4*params.num_bins));
        title('L to M drp of fitted mosaic');
        axis square; colormap bone;
        
        % L to M difference
        subplot(3,3,9);
        plot(dataValue(3*params.num_bins+1:4*params.num_bins)-currDrp(3*params.num_bins+1:4*params.num_bins), 'k'); axis square;
        title('drp_d - drp_f');
        
        % print it
        if params.save_text
            set(gcf,'PaperUnits','centimeters');
            xSize = 20; ySize=28;
            xLeft = (22-xSize)/2; yTop = (30-ySize)/2;
            set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
            if ~isempty(params.folder)
                print(88,'-dpdf', [params.base_path sprintf('/%s/%d-%s.pdf',params.folder,SEED-1,name)]);
            else
                print(88,'-dpdf', [params.base_path sprintf('/%d-%s.pdf',SEED-1,name)]);
            end
        end
    end
end

% save out drp info
drp_data.data    = dataValue;
drp_data.fitted  = currDrp;
drp_data.centers = centers;
drp_data.nBins = params.num_bins;

% save out mosaic info
cone_data.cones      = cones;
cone_data.dataLabels = cO;
cone_data.fitLabels  = c;
cone_data.seed = SEED-1;
