function [best_stacks,dists] = identify_relevant_stacks(stacks,point)
% identify_relevant_stacks      identify which image stacks contain a given point
%
% usage:  [best_stacks,dists] = identify_relevant_stacks(stacks,point)
%
% arguments:     stacks - struct array of image stacks generated in image_catalog
%                 point - 2-element vector of the point (array coordinates)
%
% outputs:     best_stacks - index number of stacks containing the point, most deeply containing first
%                    dists - distance to the nearest edge of each stack (of the stacks that contains the point)
%
%
% 2010-04  gauthier
%
%



%figure(110);clf;hold on
%plot(point(1),point(2),'r.')

% for each stack, note how deeply the point is contained in the image
for ss=1:length(stacks)

    % compute rectangle of the image in array coordinates
    wd = getfield(imfinfo(stacks(ss).im{1}),'Width');
    ht = getfield(imfinfo(stacks(ss).im{1}),'Height');
    array_rect = tforminv(stacks(ss).tforminv,[1 1;1 ht; wd ht; wd 1; 1 1]);

    % identify distance to closest edge
    dist_to_edge = Inf;
    for rr = 1:4
        dist_to_edge = min(dist_to_edge,distance_point_to_segment(point,array_rect(rr,:),array_rect(rr+1,:)));
    end

    % if point is not contained, make distance to edge Inf
    if ~inpolygon(point(1),point(2),array_rect(:,1),array_rect(:,2))
        dist_to_edge = Inf;
    end

    % note this rectangle's min dist
    min_dists(ss) = dist_to_edge;

    % plot outline
    %plot(array_rect(:,1),array_rect(:,2),'k.-')
end
axis equal

% sort stacks based on distance
[dists,best_stacks] = sort(min_dists);

% return those which contain the point
best_stacks = best_stacks(dists<Inf);
dists = dists(dists<Inf);
