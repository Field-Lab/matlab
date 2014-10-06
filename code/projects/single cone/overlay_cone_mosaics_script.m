% load two cone mosaics, overlay them

if ~exist('cones','var')
    offset = [0 0];
    rf_size = [320 320];
    switch 1
        case 1  % peach old vs peach new
            cones{1} = import_single_cone_data([],'mango');
            cones{2} = import_single_cone_data([],'mango-old');
            rf_size = [320 640];
        case 2  % plantain int 1 vs plantain int 6
            cones{1} = import_single_cone_data([],'plantain');
            cones{2} = import_single_cone_data([],'plantain-1');
            offset = [1 1];
    end
end


% plot cones
if 1
    % parameters
    bg_col = [1 1 1]*48/255;
    symbols = '..';
    sizes = [5 3];
    cols{1} = [120 1 3; 21 123 14; 2 0 219; 190 190 190]/255;
    cols{2} = [255 122 131; 92 255 85; 106 94 255; 190 190 190]/255;

    % set up figure
    figure(1);clf;
    imagesc(reshape(repmat(bg_col,prod(rf_size),1),rf_size(1),rf_size(2),3))
    axis image; hold on

    % plot mosaics
    for cc=1:length(cones)

        centers = cones{cc}.cone_centers;
        if cc == 1
            centers = centers + repmat(offset,size(centers,1),1);
        end

        reds = centers(cones{cc}.cone_types=='L',:);
        greens = centers(cones{cc}.cone_types=='M',:);
        blues = centers(cones{cc}.cone_types=='S',:);
        grays = centers(cones{cc}.cone_types=='U',:);

        plot(reds(:,1),reds(:,2),symbols(cc),'MarkerSize',sizes(cc),'color',cols{cc}(1,:))
        plot(greens(:,1),greens(:,2),symbols(cc),'MarkerSize',sizes(cc),'color',cols{cc}(2,:))
        plot(blues(:,1),blues(:,2),symbols(cc),'MarkerSize',sizes(cc),'color',cols{cc}(3,:))
        plot(grays(:,1),grays(:,2),symbols(cc),'MarkerSize',sizes(cc),'color',cols{cc}(4,:))


    end

end

print(1,'/snle/home/gauthier2/public_html/temp.pdf','-dpdf')

% compute statistics
if 1
    tic;[c,d,m] = cMIX(cones{1}.cone_centers,cones{2}.cone_centers,1, 0.1, 100,true,'icp3');toc
    
    cutoff = 1; % max distance for two cones to be considered the same
    
    % identify matches
    [ii,jj] = find(m>0.4);
    dists = zeros(length(ii),1);
    for cc=1:length(ii)
        dists(cc) = norm(cones{1}.cone_centers(ii(cc),:) - cones{2}.cone_centers(jj(cc),:));
        if dists(cc) < cutoff
            plot([cones{1}.cone_centers(ii(cc),1) cones{2}.cone_centers(jj(cc),1)],...
                [cones{1}.cone_centers(ii(cc),2) cones{2}.cone_centers(jj(cc),2)],'w')
            % highlight mismatch
            if cones{1}.cone_types(ii(cc)) ~= cones{2}.cone_types(jj(cc))
                plot(cones{1}.cone_centers(ii(cc),1),cones{1}.cone_centers(ii(cc),2),'ok')
            end
        end
    end
    
    matches = dists <= cutoff;
    
    same_color = sum(cones{1}.cone_types(ii(matches)) == cones{2}.cone_types(jj(matches)));
    fprintf('same location, same type: %d\n',same_color)
    fprintf('same location, different type: %d\n',sum(matches)-same_color)
    fprintf('cones unique to bayesian: %d\n',length(cones{1}.cone_types)-sum(matches))
    fprintf('cones unique to old method: %d\n',length(cones{2}.cone_types)-sum(matches))
    
    fprintf('\n\n')
    
end

print(1,'/snle/home/gauthier2/public_html/temp.pdf','-dpdf')
