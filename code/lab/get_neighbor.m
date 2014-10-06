function pairs = get_neighbor(datarun, cell_specification, varargin)
%Calculate paired mosaic information using Voronoi tesselation.
%from MOSAIC_INFO shlens
%
% defaults.hull = 0; exlude outside hull 
% defaults.distance = 0; exlude distance > mean(sd)*maxdistance
% defaults.maxdistance = 4; effects only when distance = 1
%
%greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('hull', 0);% 
    p.addParamValue('distance', 0);
    p.addParamValue('normdistance', 0);
    p.addParamValue('maxdistance', 4);
    p.addParamValue('plot', 0);
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;
    
   


% get cell numbers
index=get_cell_indices(datarun,cell_specification);

if params.distance
    [dist_stixel dist_sd normdist] = get_sta_fit_distance(datarun,cell_specification,cell_specification);
    if params.normdistance
        distances=dist_sd;
    else
        distances=normdist;
    end
end

% allocate memory
%distances = zeros(length(index));
types = zeros(length(index));
pairs = zeros(length(index));
mn = zeros(length(index),2);

% determine receptive field centers
mn=zeros(length(index),2);
for i=1:length(index)
    mn(i,:)=datarun.(datarun.default_sta_fits).sta_fits{index(i)}.mean;
end

%sd = mean_sd(datarun,cell_specification);

% voronoi tesselation
[V, C] = voronoin(mn);

% convex hull
hull = convhull(mn(:,1), mn(:,2));

% find vertices outside convex hull
exclude = find(~inpolygon(V(:,1), V(:,2), mn(hull,1), mn(hull,2)));


% for each pair
for i=1:length(index)
  for j=(i+1):length(index)
    
    % determine distance
    %d = datarun.(datarun.default_sta_fits).sta_fits{index(i)}.mean - datarun.(datarun.default_sta_fits).sta_fits{index(j)}.mean;
    %distances(i,j) = norm(d);
    %distances(j,i) = distances(i,j);
        
    % number of overlap
    overlap = intersect(C{i}, C{j});
  
    % remove excluded vertices (outside of mosaic)
    if params.hull
        overlap = setdiff(overlap, exclude);
    end

    % check if length is reasonable
    if (length(overlap) > 2) | (length(overlap) < 0)
      warning(['voronoi failed: ' num2str(length(overlap)) ' intersections']);
    end
    
    if params.distance
        %if distances(j,i)<sd*params.maxdistance
        if distances(j,i)<params.maxdistance
            pairs(i,j) = (length(overlap) > 0);
        end
    else
        pairs(i,j) = (length(overlap) > 0);
    end
    pairs(j,i) = pairs(i,j);
    
  end
end


if params.plot% plot voronoi tesselation
    %figure('Color','white');
    %clf
    if datarun.default_sta_fits=='obvius'
        [VX, VY] = voronoi(mn(:,1), datarun.stimulus.field_height-mn(:,2));
        h = plot(VX,VY,'-k',mn(:,1),datarun.stimulus.field_height-mn(:,2),'.k');
        set(h(1:end-1),'xliminclude','off','yliminclude','off')
        hold on;
        plot(V(exclude,1),datarun.stimulus.field_height-V(exclude,2),'ob');
        plot(mn(hull,1), datarun.stimulus.field_height-mn(hull,2))
        for i=1:length(index)
            for ii=i+1:length(index)
                if pairs(i,ii)
                    plot(mn([i ii],1), datarun.stimulus.field_height-mn([i ii],2),'r')
                end
            end
        end
        axis ij
        axis equal
    else
        [VX, VY] = voronoi(mn(:,1), mn(:,2));
        h = plot(VX,VY,'-k',mn(:,1),mn(:,2),'.k');
        set(h(1:end-1),'xliminclude','off','yliminclude','off')
        hold on;
        plot(V(exclude,1),V(exclude,2),'ob');
        plot(mn(hull,1), mn(hull,2))
        for i=1:length(index)
            for ii=i+1:length(index)
                if pairs(i,ii)
                    plot(mn([i ii],1), mn([i ii],2),'r')
                end
            end
        end
        axis ij
        axis equal        
    end
    
end











