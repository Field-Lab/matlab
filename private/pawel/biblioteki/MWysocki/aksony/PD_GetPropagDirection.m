function [line, orient, v0, ve, vep, edges, func, success, rms] = PD_GetPropagDirection (DataPath,PatternNumber,MovieNumber,BadElectrodes,InitialDirection)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % This function operates on data from one event (one pattern and one movie).
    % It finds maximum negative signal on every electrode,
    % uses these values to create 2-D interpolation (output: func)
    % and finds direction of axon signal propagation.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inputs:
    %  ----------------
    %   PatternNumber   - event pattern No.
    %   MovieNumber     - event movie No.
    %   DataPath        - data directory
    %   BadElectrodes   - broken electrodes. Data from it shouldn't be analyze
    %
    % Outputs:
    %  ----------------
    %   line    - coefficients of line, representing direction of signal propagation
    %   v0      - stimulating electrode point, projected onto the direction line
    %   ve      - unit vector along the direction line
    %   vep     - unit vector perpendicular to the direction line
    %   edges   - matrix of two row vectors, where each vector represents coordinates of one edge point
    %                (point where the direction line intersects one of MEA edges)
    %   func    - interpolating function that you can evaluate at any point on the MEA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(BadElectrodes(BadElectrodes == PatternNumber))
        line = [];
        orient = [];
        v0 = [];
        ve = [];
        vep = [];
        edges = [];
        func = [];
        success = 0;
        rms = [];
        return;
    end

    [DataTraces0,~,~]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0); % read data
    s=reshape(mean(DataTraces0),512,140); % average 51 traces for every electrode
    tmp = mean(s(:,135:end),2); % evaluate offset as mean of 6 last samples
    s = s - repmat(tmp,1,140);  % subtract the offset
    st = s(:,7:40); % skip 6 first samples (due to stimulus artifacts)
    [minimus,mintimes] = min(st,[],2);  % find minimum values on every electrode...
                                        
    maximus = -minimus;                 % ... and consider its modulus as maximum spike signals

    % get electrode map to find electrodes coordinates
    electrodeMap = edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

    elec = zeros(512,4); % electrode No.,x,y,maximum spike
    % Get electrodes positions:
    for i = 1:512
        elec(i,2)=electrodeMap.getXPosition(i);
        elec(i,3)=electrodeMap.getYPosition(i);
    end
    elec(:,4) = maximus;

    elec(:,1) = (1:512)';   % add information about electrode No.
                            % (in next line we will remove bad electrodes,
                            %  so position in the matrix will no longer
                            %  carry this information)
                            
    ES_pos = elec(PatternNumber,2:3); % stimulating electrode position on the MEA
                            
    neighbors = electrodeMap.getAdjacentsTo(PatternNumber,1)'; % get stimulating electrode neighbors (radius = 1)
    elec([BadElectrodes neighbors],:) = []; % remove BAD electrodes and stimulating electrode neighbors
    elec(elec(:,4) < 0,:) = []; % remove electrodes where maximum spike is negative
    
    tic = 10;
    xx = -945:tic:945;
    yy = -450:tic:450;
    [qx,qy] = meshgrid(xx,yy); % create grid with 'tic' um spacing
    xi = reshape(qx,[],1);  % transform grid into column vector form
    yi = reshape(qy,[],1);
    % Use high quality inerpolation to find points in the created grid:
    zi = PD_BiharmonicSplineInterp2(elec(:,2),elec(:,3),elec(:,4),xi,yi);
    %  (this algorithm gives very good spline interpolation,
    %   but it is very slow, so we use it only once)
    % ...then do fast interpolation, based on 'grid points':
    f = TriScatteredInterp(xi,yi,zi,'natural');
    % Now we evaluate values...
    z = f(qx,qy);
    % ...but we will also use the 'f' function many times later to find intersection plots

    maksimum = max(max(z));
    dist = 60;
    
    % [Direction of propagation finding algorithm]
    if isempty(InitialDirection)
        [vert_max,xpos] = max(z,[],2);
        [hor_max,ypos] = max(z,[],1);
        % vertical
        xvert = xx(xpos);
        yvert = yy;
        yvert(vert_max < maksimum/7) = [];
        xvert(vert_max < maksimum/7) = [];
        check = PD_CheckDistance(xvert,yvert,ES_pos(1),ES_pos(2),dist);
        yvert(check) = [];
        xvert(check) = [];
%         scatter(xvert,yvert)
        if length(yvert) < 5
            vert_rmsd = 1000;
        else
            vert_fit = polyfit(yvert,xvert,1);
            vert_rmsd = PD_RMSD(yvert,xvert,vert_fit);
        end
        % horizontal
        xhor = xx;
        yhor = yy(ypos);
        xhor(hor_max < maksimum/7) = [];
        yhor(hor_max < maksimum/7) = [];
        check = PD_CheckDistance(xhor,yhor,ES_pos(1),ES_pos(2),dist);
        xhor(check) = [];
        yhor(check) = [];
        if length(yhor) < 5
            hor_rmsd = 1000;
        else
            hor_fit = polyfit(xhor,yhor,1);
            hor_rmsd = PD_RMSD(xhor,yhor,hor_fit);
        end
        %scatter(xvert,yvert);
        %axis([-945 945 -450 450])
        %rms = min(hor_rmsd,vert_rmsd);
        if min(hor_rmsd,vert_rmsd) == 1000
            line = [];
            orient = [];
            v0 = [];
            ve = [];
            vep = [];
            edges = [];
            func = [];
            success = 0;
            rms = [];
            return;
        end
        if hor_rmsd < vert_rmsd
            wsp = hor_fit;
            orient = 1;
        else
            wsp = vert_fit;
            orient = 0;
        end
    else
        if abs(InitialDirection(1) > InitialDirection(2))
            orient = 1;
            wsp(1) = InitialDirection(2)/InitialDirection(1); % a = vy/vx
            wsp(2) = ES_pos(2) - wsp(1)*ES_pos(1);
        else
            orient = 0;
            wsp(1) = InitialDirection(1)/InitialDirection(2); % a = vx/vy
            wsp(2) = ES_pos(1) - wsp(1)*ES_pos(2);
        end
    end
    
    edge_points = PD_FindEdgePoints(wsp, orient);
    if length(edge_points(:,1)) < 2
        line = [];
        orient = [];
        v0 = [];
        ve = [];
        vep = [];
        edges = [];
        func = [];
        rms = [];
        success = 0;
        return;
    end
    obrot = [0 -1; 1 0];

    line_dist = 100;
    n = 0;
    while line_dist > 1 && n < 10
        v = [edge_points(2,1)-edge_points(1,1);edge_points(2,2)-edge_points(1,2)];
        v_norm = sqrt(v'*v);
        v = v/v_norm;
        vp = obrot*v;

        rmax = 200;
        dr = 10;
        v_start = edge_points(1,:)';
        points = zeros(floor(v_norm/tic)+1,4);
        for d = 0:1:floor(v_norm/tic)
           sect = PD_GetIntersection( v_start + 10*d*v, vp, rmax, dr, f );
           [~, pos] = max(sect(:,4));
           points(d+1,:) = sect(pos,:);
        end
        points(points(:,4) < maksimum/7,:) = [];
        check = PD_CheckDistance(points(:,2),points(:,3),ES_pos(1),ES_pos(2),dist);
        points(check,:) = [];
        if length(points(:,1)) < 3
            line = [];
            orient = [];
            v0 = [];
            ve = [];
            vep = [];
            edges = [];
            func = [];
            rms = [];
            success = 0;
            return;
        end
        if orient
            wsp = polyfit(points(:,2),points(:,3),1);
        else
            wsp = polyfit(points(:,3),points(:,2),1);
        end

        edge_points2 = PD_FindEdgePoints(wsp, orient);
        if length(edge_points2(:,1)) < 2
            line = [];
            orient = [];
            v0 = [];
            ve = [];
            vep = [];
            edges = [];
            func = [];
            rms = [];
            success = 0;
            return;
        end
        
        v2 = [edge_points2(2,1)-edge_points2(1,1);edge_points2(2,2)-edge_points2(1,2)];
        v_norm = sqrt(v2'*v2);
        v2 = v2/v_norm;
        line_dist = PD_LineDist( v, edge_points, v2, edge_points2 );
        edge_points = edge_points2;
        n = n + 1;
        
%         xx = -945:10:945;
%         yy = -450:10:450;
%         if orient
%            plot(xx,polyval(wsp,xx),'r')
%         else
%            plot(polyval(wsp,yy),yy,'r')
%         end
%         title(num2str(line_dist))
%         axis([-1000 1000 -500 500]);
%         hold on
%         scatter(points(:,2),points(:,3),'r')
%         pause
    end
%     n
    rms = sqrt(sum(points(:,1).^2)/length(points));
    
    if line_dist > 1
        success = 0;
    else
        success = 1;
    end
    
    % setting outputs:
    
    % przes - punkt przeciecia prostej z osia OX
    if orient
        przes = [-wsp(2)/wsp(1);0]; % x = -b/a, y = 0
    else
        przes = [wsp(2);0]; % x = b, y = 0
    end
    
    v = [edge_points(2,1)-edge_points(1,1);edge_points(2,2)-edge_points(1,2)];
    v_norm = sqrt(v'*v);
    ve = v/v_norm;
    vep = obrot*ve;
    v0 = v*(v'*(ES_pos' - przes)) + przes;
    
    line = wsp;
    edges = edge_points;
    %func = f;
    func = zi;
end