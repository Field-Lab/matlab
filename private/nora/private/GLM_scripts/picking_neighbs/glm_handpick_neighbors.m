exps = ['2012-08-09-3'; '2012-09-27-3'; '2013-08-19-6'; 
'2013-10-10-0'];
load('/Volumes/Lab/Users/akheitman/NSEM_Home/Cell_Selection/allcells.mat');
for i = 1
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
        for j = 1:n_cells
            plot(sta_centers(:,1), sta_centers(:,2), 'o');
            hold on
            plot(sta_centers(j,1), sta_centers(j,2), 'o');
            [x, y] = ginput;
            pause();
            pairs = length(x);
            neighbor = zeros(pairs,1);
            for k = 1:pairs
               temp = abs(sum(sta_centers - repmat([x(k) y(k)], n_cells, 1),2));
               [~,neighbor(k)] = min(temp);
            end
            eval(['neighbors.cell_' num2str(cellgroup(j)) '.pairs = datarun.cell_ids(cell_nums(neighbor));']);
            eval(['neighbors.cell_' num2str(cellgroup(j)) '.n_couplings = pairs;']);
        end
    end
    save(['neighbors_exp' num2str(i) '_celltype' num2str(type) '.mat'],'neighbors')
end

