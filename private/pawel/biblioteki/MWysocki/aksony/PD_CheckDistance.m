function near = PD_CheckDistance( x0, y0, xe, ye, dist)
    % This function checks if distance between
    % specified points (x0,y0) and electrode position (xe,ye)
    % is less than specified value (dist).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs:
    %  ----------------
    %   x0   - vector of x coordinates of points we are checking
    %   y0   - vector of y coordinates of points we are checking
    %   xe   - electrode x coordinate
    %   ye   - electrode y coordinate
    %   dist - minimum permissible distance between (x0,y0) points and (xe,ye)
    %
    % Outputs:
    %  ----------------
    %   near - returns 1 if distance between (x0,y0) and (xe,ye)
    %           is less than specified value (dist)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dd = arrayfun(@(x,y) PD_Dist([x, y], [xe, ye]), x0, y0);
    near = dd <= dist;
end

