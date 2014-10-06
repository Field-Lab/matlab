% s-cone efficiency script
%% new cone finding

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';

% blueberry
path_and_name{2,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{2,2} = 'blueberry';

% apricot
path_and_name{3,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{3,2} = 'apricot';

% peach
path_and_name{4,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{4,2} = 'peach';


%kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi';


% grapes
path_and_name{1,1} = '/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = 'grapes';

% apple
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
path_and_name{1,2} = 'apple';

%% old cone finding

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain-old';

% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = 'blueberry-old';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{1,2} = 'apricot-old';

%kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi-old';

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = 'peach-old';

% grapes
path_and_name{1,1} = '/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = 'grapes-old';

% apple
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
path_and_name{1,2} = 'apple-old';

%%


% or just grab all S cones
%s_cone_indices = find(datarun.cones.types == 'S');
% s_cone_centers = datarun.cones.centers(s_cone_indices,:);

cell_types = [1,2,3,4];
thresholds = [0.15 0.2];
radius_fraction = 1;

num_dset = size(path_and_name,1);
num_thresh = length(thresholds);
num_cell_types = length(cell_types);

thresh_v_efficiency = zeros(num_dset, num_thresh, num_cell_types);
thresh_v_encounters = zeros(num_dset, num_thresh, num_cell_types);
hits = zeros(num_dset, num_thresh, num_cell_types);

for dset = 1:num_dset
    
    clear datarun
    
    % load data
    datarun = load_data(path_and_name{dset,1});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = set_polarities(datarun);
    datarun = import_single_cone_data(datarun, path_and_name{dset,2});    
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

    mosaic_s_cone_encounters = zeros(length(cell_types), 1);
    mosaic_s_cone_hits = zeros(length(cell_types), 1);
    mosaic_s_cone_efficiency = zeros(length(cell_types), 1);


    % get S cone centers
    [s_weights, s_selection, s_extras] = select_cone_weights(datarun, {5}, 'thresh', 0.1,...
                                                'radius', [0 inf], 'polarity', 1,...
                                                'contiguity', false, 'scale', 20.0, 'remove_cones', 'LMU');   
    s_cone_centers = s_extras.new_datarun.cones.centers;


    for thresh = 1:length(thresholds)

        weight_fraction = thresholds(thresh);

        for ctype = 1:length(cell_types)
            s_cone_encounter = 0;
            s_cone_hits = 0;
            s_cone_fp = 0;
            cell_indices = get_cell_indices(datarun, {ctype});
            num_rgcs = length(cell_indices);
            %cell_radii = zeros(num_rgcs, 1);
            for rgc = 1:num_rgcs
                cell_fit = datarun.cones.rf_fits{cell_indices(rgc)};
                cell_radius = cell_fit.center_radius;
                cell_center = cell_fit.center;

             %   cell_radii(rgc) = cell_radius;

                s_cone_distances = ipdm(cell_center, s_cone_centers);

                near_s_cone_indices = find(s_cone_distances <= (cell_radius * radius_fraction));

                if ~isempty(near_s_cone_indices)
                    for scn = 1:length(near_s_cone_indices)
                        s_cone_encounter = s_cone_encounter + 1;
                        peak_cone_weight = max(datarun.cones.weights(:,cell_indices(rgc)));

                        s_to_cone_distances = ipdm(s_cone_centers(near_s_cone_indices(scn),:), datarun.cones.centers);
                        s_cone_pointer = find(s_to_cone_distances == 0);

                        s_cone_weight = datarun.cones.weights(s_cone_pointer, cell_indices(rgc));
                        if s_cone_weight >= (weight_fraction * peak_cone_weight)
                            s_cone_hits = s_cone_hits + 1;
                        end
                    end
                end

                % false positive check
                % check cones outside 3 sigma are not crossing threshold
                far_s_cone_indices = find(s_cone_distances > (cell_radius * 3));

                for scn = 1:length(near_s_cone_indices)
                    peak_cone_weight = max(datarun.cones.weights(:,cell_indices(rgc)));
                    s_to_cone_distances = ipdm(s_cone_centers(far_s_cone_indices(scn),:), datarun.cones.centers);
                    s_cone_pointer = find(s_to_cone_distances == 0);
                    s_cone_weight = datarun.cones.weights(s_cone_pointer, cell_indices(rgc));
                    if s_cone_weight >= (weight_fraction * peak_cone_weight)
                        s_cone_fp = s_cone_fp + 1;
                    end
                end
            end

            %hist(cell_radii, 30)
            %pause    
            s_cone_fp;
            mosaic_s_cone_encounters(ctype) = s_cone_encounter;
            mosaic_s_cone_hits(ctype) = s_cone_hits;
            mosaic_s_cone_efficiency(ctype) = s_cone_hits ./ s_cone_encounter;
        end

        thresh_v_efficiency(dset,thresh, 1:4) = mosaic_s_cone_efficiency;
        thresh_v_encounters(dset,thresh, 1:4) = mosaic_s_cone_encounters;
        hits(dset,thresh,1:4) = mosaic_s_cone_hits(1:4);
    end

end

%[phat, pci] = binofit(mosaic_s_cone_hits(1), mosaic_s_cone_encounters(1))
%% cumulate hits and encounters across data sets

total_hits = squeeze(sum(hits,1));
total_encounters = squeeze(sum(thresh_v_encounters,1));
threshold_curves = total_hits ./ total_encounters;

figure
semilogx(thresholds, threshold_curves(:,1), 'ko-')
hold on
semilogx(thresholds, threshold_curves(:,2), 'bo-')
semilogx(thresholds, threshold_curves(:,3), 'ro-')
semilogx(thresholds, threshold_curves(:,4), 'go-')
legend('on parasol','off parasol', 'on midget', 'off midget')
hold off

[phat, pci] = binofit(total_hits(9,4), total_encounters(9,4));
phat
phat - pci(1)
pci(2) - phat

%%
figure
semilogx(1./thresholds, thresh_v_efficiency(1,:) ./ max(thresh_v_efficiency(1,:)), 'ko-')
hold on
semilogx(1./thresholds, thresh_v_efficiency(2,:) ./ max(thresh_v_efficiency(2,:)), 'bo-')
semilogx(1./thresholds, thresh_v_efficiency(3,:) ./ max(thresh_v_efficiency(3,:)), 'ro-')
semilogx(1./thresholds, thresh_v_efficiency(4,:) ./ max(thresh_v_efficiency(4,:)), 'go-')
legend('on parasol','off parasol', 'on midget', 'off midget')
title([path_and_name{1,2},' normalized'])
hold off


coef = [10,5];
fit_coef_offm = nlinfit(1./thresholds', (thresh_v_efficiency(4,:) ./ max(thresh_v_efficiency(4,:)))', 'cumulative_gaussian', coef)
fit_coef_onm = nlinfit(1./thresholds', (thresh_v_efficiency(3,:) ./ max(thresh_v_efficiency(3,:)))', 'cumulative_gaussian', coef)
fit_coef_offp = nlinfit(1./thresholds', (thresh_v_efficiency(2,:) ./ max(thresh_v_efficiency(2,:)))', 'cumulative_gaussian', coef)
fit_coef_onp = nlinfit(1./thresholds', (thresh_v_efficiency(1,:) ./ max(thresh_v_efficiency(1,:)))', 'cumulative_gaussian', coef)





