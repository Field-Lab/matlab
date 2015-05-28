% NB 2015-04-29
% This just picks the neighboring cells
% Will only run this once, and save the output in a mat file in the same
% folder as allcells.mat
exps = ['2012-08-09-3'; '2012-09-27-3'; '2013-08-19-6'; 
'2013-10-10-0'];

load('/Volumes/Lab/Users/akheitman/NSEM_Home/Cell_Selection/allcells.mat');
for i = 2
    exp_info = experimentinfoNSEM(exps(i,:));
    datarun = load_data([exps(i,:) '/' exp_info.dr.mas]);
    datarun = load_params(datarun);
    for type = 1
        if type == 1
            cellgroup = allcells{i}.ONP;
        else
            cellgroup = allcells{i}.OFFP;
        end
        [cell_nums,~,~] = get_cell_indices(datarun, cellgroup);
        n_cells = length(cell_nums);
        sta_centers = zeros(n_cells, 2);
        for j = 1:n_cells
            sta_centers(j,:) = datarun.vision.sta_fits{cell_nums(j)}.mean;
        end
        nearest_neighbor = zeros(n_cells,1);
        cell_dist = zeros(n_cells, n_cells);
        for j = 1:n_cells
            for k = 1:n_cells
                cell_dist(j,k) = norm(sta_centers(j,:) - sta_centers(k,:));
            end
            [~,idx] = sort(cell_dist(j,:));
            nearest_neighbor(j) = cell_dist(j,idx(2)); % 2 because idx(1) will be 0, the cell with itself
        end
        % exp1, type1 = 1.3
        cutoff = 1*max(nearest_neighbor);
        
        %testing
        for j = 1:n_cells
            temp = cell_dist(j,:);
            neighbs = find(temp < cutoff);
            plot(sta_centers(:,1), sta_centers(:,2), 'o');
            hold on
            plot(sta_centers(neighbs,1), sta_centers(neighbs,2), 'o');
            plot(sta_centers(j,1), sta_centers(j,2), 'o');
            hold off
            pause
        end
    end
end
