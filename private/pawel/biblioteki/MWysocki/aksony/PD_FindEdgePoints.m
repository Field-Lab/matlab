function [ edge ] = PD_FindEdgePoints( wsp, orient )
%  wsp
    edge_points = zeros(4,3);

    % right edge
    x = 945;
    if orient
        y = wsp(1)*x + wsp(2);
    else
        y = (x - wsp(2))/wsp(1);
    end
    if y >= -450 && y <= 450
%         w = 'right'
        edge_points(1,1) = x;
        edge_points(1,2) = y;
        edge_points(1,3) = 1;
    end

    % bottom edge
    y = -450;
    if orient
        x = (y - wsp(2))/wsp(1);
    else
        x = wsp(1)*y + wsp(2);
    end
    if x >= -945 && x <= 945
%         w = 'bottom'
        edge_points(2,1) = x;
        edge_points(2,2) = y;
        edge_points(2,3) = 1;
    end

    % left edge
    x = -945;
    if orient
        y = wsp(1)*x + wsp(2);
    else
        y = (x - wsp(2))/wsp(1);
    end
    if y >= -450 && y <= 450
%         w = 'left'
        edge_points(3,1) = x;
        edge_points(3,2) = y;
        edge_points(3,3) = 1;
    end

    % top edge
    y = 450;
    if orient
        x = (y - wsp(2))/wsp(1);
    else
        x = wsp(1)*y + wsp(2);
    end
    if x >= -945 && x <= 945
%         w = 'top'
        edge_points(4,1) = x;
        edge_points(4,2) = y;
        edge_points(4,3) = 1;
    end

    edge_points(edge_points(:,3)==0,:) = [];
    edge = unique(edge_points(:,1:2),'rows');
    if length(edge) > 2
        edge(3:end,:) = [];
    end
end

