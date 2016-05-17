function contours = hex_contour(xCoords, yCoords, zVals, nContours, varargin)

% HEX_CONTOUR      finds contour paths based on values in a hexagonal lattice (electrode array)
%
% The algorithm used is analagous to the algorithm used by the built-in matlab function CONTOUR
%
% note: this function has not been tested on 512/519 datasets (only 61)
%
% note: hexagons specified by coordinates must be close to regular hexagons (i.e. longer sides must
%       be no more than 20% longer than shorter sides length sides
%
% usage:  contours = hex_contour(xCoordinates, yCoordinates, values, nContours)
%
% arguments:
%        xCoordinates - vector of x coordinates for hexagonal lattice
%        yCoordinates - vector of y coordinates for hexagonal lattice
%              values - vector of amplitude values that contours will be calculated from (should be
%                       same length as xCoordinates and yCoordinates)
%           nContours - integer specifying desired number of contour levels
%                cmap - matlab string for a colormap to use, default gray
%
%
% outputs:
%            contours - nContours x 1 struct array containing contour paths:
%                    contours(i).elevation - amplitude at which ith contour is drawn
%                     contours(i).paths{j} - n x 2 matrix of x,y coordinates defining the jth path 
%                                            for contour level i
%
% optional params:
%
% contourSpacing    'linear'            spacing at which contour levels are chosen ('linear' or
%                                           'log')
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure.
%                                           if empty, don't plot.  if -1, plot in current.
%
% author: Lauren Jepson
% date: 2010-02-03

% use this command for easy plotting:
% hex_contour(datarun.ei.position(:,1), datarun.ei.position(:,2), max(abs(get_ei(datarun,cell_id)),[],2), 10)



p = inputParser;

% optional parameters
p.addParamValue('contourSpacing', 'linear', @(x)any(strcmpi(x,{'linear', 'log'})));
p.addParamValue('fig_or_axes', []);
p.addParamValue('plotCoords', true, @islogical) %if true, plots black dots at electrode locations
p.addParamValue('cmap','gray',@ischar); 
% resolve user input and default values
p.parse(varargin{:});

% get params struct
contSpacing = p.Results.contourSpacing;
fig_or_axes = p.Results.fig_or_axes;
plotElecs = p.Results.plotCoords;
cmap = p.Results.cmap;

%convert xCoords and yCoords to column vectors if not already
if size(xCoords, 2) > 1
    xCoords = xCoords';
end
if size(yCoords, 2) > 1
    yCoords = yCoords';
end
bothCoords = [xCoords yCoords];

nElec = length(xCoords);
edgeLength = min(pdist(bothCoords));

%finds nearest-neighbor electrode pairs (edges)
edgePairs = pdist(bothCoords) < 1.2*edgeLength; %1.2 to account for rounding errors
nEdges = sum(edgePairs);

edges = zeros(nEdges, 2); %electrode identities in each pair
iPair = 0;
iEdge = 0;
%increments through electrode pairs in the same order that pdist outputs the distances between them
for i = 1:nElec-1 %first value in pair
    for j = i+1:nElec %second value in pair
        iPair = iPair+1;
        if edgePairs(iPair)
            iEdge = iEdge+1;
            edges(iEdge, :) = [i, j];
        end
    end
end

edgeElectrodesBin = zeros(nElec, 1);
for i = 1:nElec
    if sum(sum(edges == i)) > 0 && sum(sum(edges == i)) < 6
        edgeElectrodesBin(i) = 1;
    end
end

%determines values of contours
zValsNoZeros = zVals(zVals ~= 0); %if there are non-electrodes in the list (e.g. 9, 25, 57)
if strcmpi(contSpacing, 'log')
    zRange = [log(min(zValsNoZeros)) log(max(zValsNoZeros))];
    cInt = (zRange(2) - zRange(1))/nContours;
    contourVals = zRange(1)+0.5*cInt:cInt:zRange(2)-0.4*cInt;
    contourVals = exp(contourVals);
else
    zRange = [min(zValsNoZeros) max(zValsNoZeros)];
    cInt = (zRange(2) - zRange(1))/nContours;
    contourVals = zRange(1)+0.5*cInt:cInt:zRange(2)-0.4*cInt;
end



contourPaths = cell(nContours, 1);
loopFlags = cell(nContours, 1);

%finding the locations of the contours!!!
for i = 1:nContours
    contourPaths{i} = [];
    edgeCrossings = zeros(nEdges, 1);
    for j = 1:nEdges
        %if contour crosses edge
        if (zVals(edges(j,1)) <= contourVals(i) && zVals(edges(j,2)) > contourVals(i)) ||...
                (zVals(edges(j,1)) > contourVals(i) && zVals(edges(j,2)) <= contourVals(i))
            edgeCrossings(j) = 1;
        end
    end
    
    edgeCrossingsUnvisited = edgeCrossings;
    %start at edge that is on edge of array, if one exists
    foundEdge = 0;
    for j = 1:nEdges
        if edgeCrossings(j) && edgeElectrodesBin(edges(j,1)) && edgeElectrodesBin(edges(j,2))
            currentEdge = j;
            foundEdge = 1;
            loopFlags{i} = 0;
            break
        end
    end
    
    if ~foundEdge %just pick any of the edge crossings
        currentEdge = find(edgeCrossings, 1);
        loopFlags{i} = 1;
    end
    
    contourOpen = 1;
    while contourOpen
        edgeElecs = edges(currentEdge,:);
        
        %calculates location along edge where contour crosses
        t = (contourVals(i) - zVals(edgeElecs(1)))/(zVals(edgeElecs(2)) - zVals(edgeElecs(1)));
        xPos = xCoords(edgeElecs(1)) + t*(xCoords(edgeElecs(2)) - xCoords(edgeElecs(1)));
        yPos = yCoords(edgeElecs(1)) + t*(yCoords(edgeElecs(2))-yCoords(edgeElecs(1)));
        contourPaths{i} = [contourPaths{i}; xPos yPos];
        edgeCrossingsUnvisited(currentEdge) = -1;
        
        %locates next edge with contour-crossing
        if sum(edgeCrossingsUnvisited == 1) == 0 %no edge-crossings that haven't already been examined
            contourOpen = 0;
        else
            neighborElecsBoth = [];
            for j = 1:nElec
                if ~isnan(xCoords(j))
                    length1 = norm(bothCoords(j, :) - bothCoords(edgeElecs(1), :));
                    length2 = norm(bothCoords(j, :) - bothCoords(edgeElecs(2), :));
                    if length1 < 1.2*edgeLength && length2 < 1.2*edgeLength && j~=edgeElecs(1) && j~=edgeElecs(2)
                        neighborElecsBoth = [neighborElecsBoth j]; %#ok<AGROW> %should be 1 or 2 electrodes
                    end
                end
            end
            
            
            if length(neighborElecsBoth) == 1 %current edge is on edge of array
                nextEdge1 = find(max((edges(:,1) == edgeElecs(1)) .* (edges(:,2) == neighborElecsBoth),...
                    (edges(:,2) == edgeElecs(1)) .* (edges(:,1) == neighborElecsBoth))); %should only be one of the edges
                nextEdge2 = find(max((edges(:,1) == edgeElecs(2)) .* (edges(:,2) == neighborElecsBoth),...
                    (edges(:,2) == edgeElecs(2)) .* (edges(:,1) == neighborElecsBoth)));
                
                if edgeCrossingsUnvisited(nextEdge1)==1 && edgeCrossingsUnvisited(nextEdge2)==0
                    currentEdge = nextEdge1;
                elseif edgeCrossingsUnvisited(nextEdge2)==1 && edgeCrossingsUnvisited(nextEdge1)==0
                    currentEdge = nextEdge2;
                elseif edgeCrossingsUnvisited(nextEdge1)==-1 || edgeCrossingsUnvisited(nextEdge2)==-1 &&...
                        ~(edgeCrossingsUnvisited(nextEdge1)==-1 && edgeCrossingsUnvisited(nextEdge2)==-1)
                    contourPaths{i} = [contourPaths{i}; inf inf]; %flags location in contourPaths to indicate break in path
                    % look for a new contour start-point on edge of array
                    foundEdge = 0;
                    for j = 1:nEdges
                        if edgeCrossingsUnvisited(j)==1 && edgeElectrodesBin(edges(j,1)) && edgeElectrodesBin(edges(j,2))
                            currentEdge = j;
                            foundEdge = 1;
                            loopFlags{i} = [loopFlags{i} 0];
                            break
                        end
                    end
                    if ~foundEdge %just pick any of the edge crossings
                        currentEdge = find(edgeCrossingsUnvisited == 1, 1);
                        loopFlags{i} = [loopFlags{i} 1];
                    end
                else
                    error('There is a problem with the logic--fix code!!!!')
                end
            elseif length(neighborElecsBoth) == 2 %current edge is not on edge of array
                nextEdge1 = find(max((edges(:,1) == edgeElecs(1)) .* (edges(:,2) == neighborElecsBoth(1)),...
                    (edges(:,2) == edgeElecs(1)) .* (edges(:,1) == neighborElecsBoth(1)))); %should only be one of the edges
                nextEdge2 = find(max((edges(:,1) == edgeElecs(2)) .* (edges(:,2) == neighborElecsBoth(1)),...
                    (edges(:,2) == edgeElecs(2)) .* (edges(:,1) == neighborElecsBoth(1))));
                nextEdge3 = find(max((edges(:,1) == edgeElecs(1)) .* (edges(:,2) == neighborElecsBoth(2)),...
                    (edges(:,2) == edgeElecs(1)) .* (edges(:,1) == neighborElecsBoth(2))));
                nextEdge4 = find(max((edges(:,1) == edgeElecs(2)) .* (edges(:,2) == neighborElecsBoth(2)),...
                    (edges(:,2) == edgeElecs(2)) .* (edges(:,1) == neighborElecsBoth(2))));
                if edgeCrossingsUnvisited(nextEdge1)==1 && edgeCrossingsUnvisited(nextEdge2)==0
                    currentEdge = nextEdge1;
                elseif edgeCrossingsUnvisited(nextEdge2)==1 && edgeCrossingsUnvisited(nextEdge1)==0
                    currentEdge = nextEdge2;
                elseif edgeCrossingsUnvisited(nextEdge3)==1 && edgeCrossingsUnvisited(nextEdge4)==0
                    currentEdge = nextEdge3;
                elseif edgeCrossingsUnvisited(nextEdge4)==1 && edgeCrossingsUnvisited(nextEdge3)==0
                    currentEdge = nextEdge4;
                elseif (edgeCrossingsUnvisited(nextEdge1)==-1 || edgeCrossingsUnvisited(nextEdge2)==-1 &&...
                        ~(edgeCrossingsUnvisited(nextEdge1)==-1 && edgeCrossingsUnvisited(nextEdge2)==-1)) &&...
                        (edgeCrossingsUnvisited(nextEdge3)==-1 || edgeCrossingsUnvisited(nextEdge4)==-1 &&...
                        ~(edgeCrossingsUnvisited(nextEdge3)==-1 && edgeCrossingsUnvisited(nextEdge4)==-1))
                    
                    contourPaths{i} = [contourPaths{i}; inf inf]; %flags location in contourPaths to indicate break in path
                    % look for a new contour start-point on edge of array
                    foundEdge = 0;
                    for j = 1:nEdges
                        if edgeCrossingsUnvisited(j)==1 && edgeElectrodesBin(edges(j,1)) && edgeElectrodesBin(edges(j,2))
                            currentEdge = j;
                            foundEdge = 1;
                            loopFlags{i} = [loopFlags{i} 0];
                            break
                        end
                    end
                    if ~foundEdge %just pick any of the edge crossings
                        currentEdge = find(edgeCrossingsUnvisited==1, 1);
                        loopFlags{i} = [loopFlags{i} 1];
                    end
                else
                    error('There is a problem with the logic--fix code!!!!')
                end
            else %something is wrong--should only be 1 or 2 shared neighbors for any pair of electrodes
                error('logic is broken')
            end
        end
    end
end


contours = struct('elevation', 0, 'paths', cell(nContours, 1));

for i = 1:nContours
    contours(i).elevation = contourVals(i);
    
    pathInd = 1;
    contours(i).paths{pathInd} = [];
    lastPathStart = 1;
    for j = 1:length(contourPaths{i})-1
        if ~any(isinf(contourPaths{i}(j, :))) && ~any(isinf(contourPaths{i}(j+1, :)))
            contours(i).paths{pathInd} = [contours(i).paths{pathInd}; contourPaths{i}(j,1), contourPaths{i}(j,2)];
        elseif any(isinf(contourPaths{i}(j+1, :)))
            contours(i).paths{pathInd} = [contours(i).paths{pathInd}; contourPaths{i}(j,1), contourPaths{i}(j,2)];
            if loopFlags{i}(pathInd)

                contours(i).paths{pathInd} = [contours(i).paths{pathInd}; contourPaths{i}(lastPathStart,1), contourPaths{i}(lastPathStart,2)];
            end
            pathInd = pathInd+1;
            contours(i).paths{pathInd} = [];
            lastPathStart = j+2;
        end
        if j == length(contourPaths{i}) - 1 && loopFlags{i}(pathInd)
            contours(i).paths{pathInd} = [contours(i).paths{pathInd}; contourPaths{i}(j+1,1), contourPaths{i}(j+1,2)];
            contours(i).paths{pathInd} = [contours(i).paths{pathInd}; contourPaths{i}(lastPathStart,1), contourPaths{i}(lastPathStart,2)];
        elseif j == length(contourPaths{i})-1
            contours(i).paths{pathInd} = [contours(i).paths{pathInd}; contourPaths{i}(j+1,1), contourPaths{i}(j+1,2)];
        end
    end
end

eval(sprintf('contColors = flipud(%s(nContours));',cmap)); 

plot_axes = set_up_fig_or_axes(fig_or_axes);
if ~isempty(plot_axes)

    axes(plot_axes)
    hold on
    % for i = 1:64
    %     plot(xCoords(i), yCoords(i), 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')
    % end
    % plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')

    %plots electrode positions
    if plotElecs
        for i = 1:length(xCoords)
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')
        end
    end

    %plots contours
    for i = 1:nContours
        for j = 1:length(contours(i).paths)
            plot(contours(i).paths{j}(:,1), contours(i).paths{j}(:,2), 'color', contColors(i,:))
        end
    end
    hold off

    set(gcf, 'color', [1 1 1]);
    % set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
    axis equal
    axis off

end



