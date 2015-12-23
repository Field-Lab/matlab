%%
if ~exist('match_list', 'var')
    load match;
end

%%
datarun.names.rrs_ei_path = fullfile(server_path(), '2007-09-18-4', 'data002-nwpca-duplicates/data002-nwpca-duplicates.ei');
datarun = calc_csd(datarun, 'cellspec', match_list(:,2));

%%
switch 2
    case 1
        switch 2
            case 1
                im = fstitch_arr;
                imxoff = -fstitch_xd(1);
                imyoff = -fstitch_yd(1);
                axon_color = 'w';
                imshow(im);
                axis xy;
                hold on;
            case 2
                imxoff = 0;
                imyoff = 0;
                axon_color = 'k';
        end
        
        for i = 1:length(axons);
            off_axons{i}(:,1) = axons{i}(:,1) + imxoff;
            off_axons{i}(:,2) = axons{i}(:,2) + imyoff;
        end
        
        pos = datarun.ei.position;
        pos(:,1) = pos(:,1) + imxoff;
        pos(:,2) = pos(:,2) + imyoff;
            
        % for each match, plot the axon and EI points
        for mm = 1:size(match_list,1)
            figure; axis equal; 
            hold on;
            
            % add electrodes
            plot(pos(:,1),pos(:,2),'o','color',[.5 .5 0])
            axon_id = match_list(mm,1);
            axon = off_axons{axon_id};
            points = traced_cell_points(axon(1:2,:), axon(2:end,:));
            plot(points(:,1), points(:,2), axon_color);

            cell_id = match_list(mm,2);
            ei = get_ei(datarun, cell_id);
            [poslocs, pos_tracex, pos_tracey] = flipflop_ei_cat2(ei, pos);
%            plot(pos_tracex', pos_tracey');
            
            set(gca, 'ColorOrder', [linspace(1, 0, length(pos_tracex)); linspace(0,0,length(pos_tracex)); linspace(0,1,length(pos_tracex))]');
            plot(pos_tracex, pos_tracey, '.', 'MarkerSize', 20);
            locs = find(cell2mat(collect(poslocs, @(x)(~isempty(x)))));
            first = locs(1);
            second = locs(2);
            plot(pos_tracex(:,first), pos_tracey(:,first), 'ro', 'MarkerSize', 10);
            plot(pos_tracex(:,second), pos_tracey(:,second), 'ms', 'MarkerSize', 10);
        end
        
    case 2        
        % for each match, plot the axon and EI
        for mm = 1:size(match_list,1)
            axon_id = match_list(mm,1);
            cell_id = match_list(mm,2);
            
            switch 4
                case 0
                    csd = get_csd(datarun, cell_id);
                    csd(csd > 0) = 0;
                    figure(cell_id);
                    plot_ei_(csd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', 1.5);
                    
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    hold on;
                    plot(points(:,1),points(:,2),'g');
                case 1
                    ei = get_ei(datarun, cell_id);
                    figure;
                    plot_ei_(ei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0]);
                    
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    hold on;
                    plot(points(:,1),points(:,2),'g');
                case 2
                    ei = get_ei(datarun, cell_id);
                    [firstposes, secondposes, firstnegs, negs, elec_colors] = flipflop_ei_cat(ei);
                    figure;
                    plot_ei_(ei, datarun.ei.position, 0, 'elec_colors', elec_colors);
                    
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    hold on;
                    plot(points(:,1),points(:,2),'g');
                case 3
                    plot_ei_scroll(datarun, cell_id, 'scale', 3, 'axon', axons{axon_id}, 'figure', cell_id);
                case 4
                    scales = [1 2 3];
                    nscales = length(scales);
                    figure;
                    for i = 1:nscales
                        scale = scales(i);
                        
                        ei = get_ei(datarun, cell_id);
                        ei = ei(:,1:end-1);
                        ax = sanesubplot(3, nscales, {1 i});
                        hold on;
                        plot_ei_(ei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale);
                        axon = axons{axon_id};
                        points = traced_cell_points(axon(1:2,:), axon(2:end,:));
%                         plot(points(:,1),points(:,2),'g');
                        title(sprintf('scale %.2f', scale));
                        set(ax, 'XTick', [], 'YTick', []);
                        if i == 1, ylabel('EI'); end
                        
                        csd = get_csd(datarun, cell_id);
                        ax = sanesubplot(3, nscales, {2 i});
                        hold on;
                        plot_ei_(csd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale);
                        axon = axons{axon_id};
                        points = traced_cell_points(axon(1:2,:), axon(2:end,:));
%                         plot(points(:,1),points(:,2),'g');
                        set(ax, 'XTick', [], 'YTick', []);
                        if i == 1, ylabel('CSD'); end

                        negcsd = csd;
                        negcsd(negcsd > 0) = 0;
                        ax = sanesubplot(3, nscales, {3 i});
                        hold on;
                        plot_ei_(negcsd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale);
                        axon = axons{axon_id};
                        points = traced_cell_points(axon(1:2,:), axon(2:end,:));
%                         plot(points(:,1),points(:,2),'g');
                        set(ax, 'XTick', [], 'YTick', []);
                        if i == 1, ylabel('Neg rectified CSD'); end
                    end
                case 5
                    scale = 0.5;
                    roots = [1 0.75];
                    nroots = length(roots);
                    figure;
                    for i = 1:nroots
                        root = roots(i);
                        
                        ei = get_ei(datarun, cell_id);
                        ei = ei(:,1:end-1);
                        rootei = ((abs(ei) + 1).^root - 1) .* sign(ei);
                        ax = sanesubplot(3, nroots, {1 i});
                        hold on;
                        plot_ei_(rootei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', (scale+1).^root - 1);
                        axon = axons{axon_id};
                        points = traced_cell_points(axon(1:2,:), axon(2:end,:));
%                         plot(points(:,1),points(:,2),'g');
                        title(sprintf('root %.2f', root));
                        set(ax, 'XTick', [], 'YTick', []);
                        if i == 1, ylabel('EI'); end
                        
                        csd = get_csd(datarun, cell_id);
                        rootcsd = ((abs(csd) + 1).^root - 1) .* sign(csd);
                        ax = sanesubplot(3, nroots, {2 i});
                        hold on;
                        plot_ei_(rootcsd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', (scale+1).^root - 1);
                        axon = axons{axon_id};
                        points = traced_cell_points(axon(1:2,:), axon(2:end,:));
%                         plot(points(:,1),points(:,2),'g');
                        set(ax, 'XTick', [], 'YTick', []);
                        if i == 1, ylabel('CSD'); end

                        negcsd = rootcsd;
                        negcsd(negcsd > 0) = 0;
                        ax = sanesubplot(3, nroots, {3 i});
                        hold on;
                        plot_ei_(negcsd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', (scale+1).^root - 1);
                        axon = axons{axon_id};
                        points = traced_cell_points(axon(1:2,:), axon(2:end,:));
%                         plot(points(:,1),points(:,2),'g');
                        set(ax, 'XTick', [], 'YTick', []);
                        if i == 1, ylabel('Neg rectified CSD'); end
                    end
                case 6
                    scale = 1;
                    figure;
                    
                    ei = get_ei(datarun, cell_id);
                    ei = ei(:,1:end-1);
                    ax = sanesubplot(3, 2, {1 1});
                    hold on;
                    plot_ei_(ei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    ylabel('EI');
                    title('Uncompressed EI');
                    
                    csd = get_csd(datarun, cell_id);
                    ax = sanesubplot(3, 2, {2 1});
                    hold on;
                    plot_ei_(csd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    title('Uncompressed CSD');
                    
                    negcsd = csd;
                    negcsd(negcsd > 0) = 0;
                    ax = sanesubplot(3, 2, {3 1});
                    hold on;
                    plot_ei_(negcsd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    title('Uncompressed RectCSD');
                    
                    lnei = log(abs(ei) + 1) .* sign(ei);
                    ax = sanesubplot(3, 2, {1 2});
                    hold on; cla
                    plot_ei_(lnei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', log(scale+1), 'cutoff', 0.65);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    ylabel('EI');
                    title('ln(EI)');

                    logei = log10(abs(ei) + 1) .* sign(ei);
                    ax = sanesubplot(3, 2, {2 2});
                    hold on;
                    plot_ei_(logei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', log10(scale+1), 'cutoff', 0.3);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    ylabel('EI');
                    title('log10(EI)');
                    
                case 7
                    figure(cell_id);
                    
                    ei = get_ei(datarun, cell_id);
                    ei = ei(:,1:end-1);
                    ax = sanesubplot(3, 2, {1 1});
                    hold on;
                    plot_ei_(ei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', 1.5);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    title('EI scale 1.5');
                    
                    csd = get_csd(datarun, cell_id);
                    ax = sanesubplot(3, 2, {2 1});
                    hold on;
                    plot_ei_(csd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', 1.5);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    title('CSD scale 1.5');
                    
                    negcsd = csd;
                    negcsd(negcsd > 0) = 0;
                    ax = sanesubplot(3, 2, {3 1});
                    hold on;
                    plot_ei_(negcsd, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', 1.5);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    title('RectCSD scale 1.5');
                    
                    clip = 0.1;
                    scale = 1;
                    ax = sanesubplot(3, 2, {1 2});
                    hold on; cla
                    plot_ei_(ei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale, 'cutoff', clip);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    ylabel('EI');
                    title(sprintf('EI clip %.2f scale %.2f', clip, scale));

                    clip = 0.15;
                    scale = 1;
                    ax = sanesubplot(3, 2, {2 2});
                    hold on; cla
                    plot_ei_(ei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale, 'cutoff', clip);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    ylabel('EI');
                    title(sprintf('EI clip %.2f scale %.2f', clip, scale));
                    
                    clip = 0.16;
                    scale = 0.9;
                    root = 0.75;
                    rootei = ((abs(ei) + 1).^root - 1) .* sign(ei);
                    ax = sanesubplot(3, 2, {3 2});
                    hold on; cla
                    plot_ei_(rootei, datarun.ei.position, 0, 'neg_color', [0 0 1], 'pos_color', [1 0 0], 'scale', scale, 'cutoff', clip);
                    axon = axons{axon_id};
                    points = traced_cell_points(axon(1:2,:), axon(2:end,:));
                    %                         plot(points(:,1),points(:,2),'g');
                    set(ax, 'XTick', [], 'YTick', []);
                    ylabel('EI');
                    title(sprintf('EI^{%.2f} clip %.2f scale %.2f', root, clip, scale));                    
            end
        end
end
