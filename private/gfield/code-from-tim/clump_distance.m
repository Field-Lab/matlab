clear all;

% set parameters here
START_SEED = 11111;
LABEL_COLUMN = 4;
N_ITERATIONS = 100;
NUM_BINS = 20;
PLOT_BINS = 20;
BIN_WIDTH = 1;

% true means we include S and U cones in scale constant computation
EJ_MODE = true;

% add to end of plot filenames
% format will be a 2 page eps file with 
% the name 'datasetname-suffix.eps
suffix = 'dot75';

% get datasets to use
ALL = 0;
[path, name] = clump_get_dataset(18);

% choose the datset
for ds = 1:length(name)
    f = figure(1); clf;
    g = figure(2); clf;
    for roi=0:0
        cones = importdata(path{ds});
        drp = cell(N_ITERATIONS,2);
        centers = cell(2,1);

        % restrict to ROI
        if roi==1
            cones = cones(cones(:,8) == 1,:);
            %cones = cones(c4.roi,:);
            offset = 2;
        else
            offset = 0;
        end
        
        % compute nnd with all cone types (like ej)
        if EJ_MODE
            c = cell(3,1);
            for cone = [1 2 3 4]
                c{cone}=find(cones(:,LABEL_COLUMN) == cone-1);
            end
            fd = clump_scale_constant(cones);
        end

        % annihilate s and u cones
        cones(cones(:,LABEL_COLUMN) == 0,:) = [];
        cones(cones(:,LABEL_COLUMN) == 3,:) = [];

        % get cones
        c = cell(3,1);
        for cone = [1 2]
            c{cone}=find(cones(:,LABEL_COLUMN) == cone);
        end
        % estimate fundamental distance between cones
        % excluding the s and u cones
        if ~EJ_MODE
            fd = fundamental_distance(cones,c);
        end

        % plot the cone mosaic
        figure(f)
        subplot(1,4,1+offset); hold on;
        plot(cones(c{1},2),cones(c{1},3),'r.');% L
        plot(cones(c{2},2),cones(c{2},3),'g.');% M
        axis image;

        histVector = zeros(N_ITERATIONS-1,1);
        SEED = START_SEED;
        for iteration=1:N_ITERATIONS
            % permute cones
            if iteration > 1
                % set up the random number generator.
                rand('twister', SEED);
                SEED = SEED + 1;
                % permute the cone labels (column 4).
                nLabels = length(cones(:,LABEL_COLUMN));
                cones(:,LABEL_COLUMN) = cones(randperm(nLabels),LABEL_COLUMN);
                % re-get cone indices
                for cone = [1 2]
                    c{cone}=find(cones(:,LABEL_COLUMN) == cone);
                end
            end
            s = 0;
            % compute distances
            count = 0;
            for type = 1:2
                locations = [cones(c{type},2) cones(c{type},3)];
                for type2 = 1:2
                    count = count + 1;
                    % compute drp
                    if type == type2
                        if iteration > 1
                            drp{iteration,count} = density_recovery_profile(locations, NUM_BINS, BIN_WIDTH);
                        else
                            [drp{iteration,count}, centers{count}] = density_recovery_profile(locations, NUM_BINS, BIN_WIDTH);
                        end
                    else
                        locationsPair = [cones(c{type2},2) cones(c{type2},3)];
                        if iteration > 1
                            drp{iteration,count} = density_recovery_profile(locationsPair, NUM_BINS, BIN_WIDTH, 'reference_centers', locations);
                        else
                            [drp{iteration,count}, centers{count}] = density_recovery_profile(locationsPair, NUM_BINS, BIN_WIDTH, 'reference_centers', locations);
                        end
                    end
                end
                distances = pdist(locations);

                % normalize distances
                distances = exp(-((distances/fd).^2)/2);

                % add up matrix
                s = s + sum(sum(distances));
            end
            if iteration == 1
                dataValue = s;
            else
                histVector(iteration-1) = s;
            end
        end
        
        % plot permutations behind data
        figure(g);
        v = cell(2,1); r = cell(2,1);
        v{1} = 'r'; r{1} = '(all)';
        v{2} = 'g'; r{2} = '(roi)';
        t = 'lm';
        count = 0;
        for type = 1:2
            for type2 = 1:2
                count = count + 1;
                if roi == 1
                    ff = 5;
                else
                    ff = 0;
                end
                subplot(2,5,count+ff); hold on;
                for iteration = 2:N_ITERATIONS-1
                    plot(centers{count}(1:PLOT_BINS),drp{iteration,count}(1:PLOT_BINS),v{type});
                end
                plot(centers{count}(1:PLOT_BINS),drp{1,count}(1:PLOT_BINS),'b');
                title(sprintf('%s to %s %s',t(type),t(type2),r{roi+1}))
            end
        end
        
        % plot subtracted plot
        subplot(2,5,5+ff); hold on;
        for iteration = 2:N_ITERATIONS-1
            plot(centers{count}(1:PLOT_BINS),drp{iteration,1}(1:PLOT_BINS)-drp{iteration,2}(1:PLOT_BINS),'r');
        end
        plot(centers{count}(1:PLOT_BINS),drp{1,1}(1:PLOT_BINS)-drp{1,2}(1:PLOT_BINS),'b');
        
        disp(dataValue)
        figure(f); hold on;
        subplot(1,4,2+offset); hold on; hist(histVector);
        limits = ylim;
        plot([dataValue dataValue], [limits(1) limits(2)], 'r', 'LineWidth', 2);
        title(sprintf('%s - all neighbors (fd = %.2f)',name{ds},fd));
        print(figure(f), '-dpsc2', sprintf('%s-%s.eps',name{ds},suffix));
        print(figure(g), '-append', '-dpsc2', sprintf('%s-%s.eps',name{ds},suffix));
        
        pause(3)
    end
end




