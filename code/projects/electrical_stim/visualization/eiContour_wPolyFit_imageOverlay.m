function eiContour_wPolyFit_imageOverlay(eiAmps,imageData, electrodeCoords, varargin)
% EICONTOUR_WPOLYFIT_IMAGEOVERLAY() plots an EI as a contour plot with a 
% polynomial fit of the axon both overlayed on an image
%  inputs: eiAmps:
%          imageData: image (RGB or grayscale) to show underneath 
%                   stimulation data 
%          electrodeCoords: electrode coordinates with respect to the image
%            
%   optional:  saveFiles
%              printElecs
%              threshForContours apply threshold to eiAmps so that the
%              contour lines arebetter scaled, default is 4 set to 1 for
%              no thresh
%              eiFitThresh minimum value required to include an electrode
%              in the polynomial fit for the ei. default is 6
%              showSoma   default is true, plots the location for the soma
%              figureNum default creates a new figure, else specify the
%              figure number (*not* handle) in which to plot
%              plotCoords %if true, plots black dots at electrode locations
%              % must be logical true or false
%  outputs:  XI: vector of 1-D look-up table "x" points
%            YI: vector of 1-D look-up table "y" points
%
% usage: eiContour_wPolyFit_imageOverlay(eiAmps,registered_tubulin, newXYCoords)
%
% Lauren Grosberg 5/2016 


% Set up defaults for optional parameters
printElecs = 0 ;
saveFiles = 0;
threshForContours = 4; 
numContourLevels = 8; 
eiFitThresh = 6;  % Default threshold used for calculating polynomial fits to the EI
showSoma = 'true'; 
figureNum = 0; 
plotCoords = true; 
ei_thresh = 0;
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'printelecs'
            printElecs = varargin{j*2};
        case 'savefiles'
            saveFiles = varargin{j*2};
        case 'threshforcontours'
            threshForContours = varargin{j*2};
        case 'numcontourlevels'
            numContourLevels = varargin{j*2};
        case 'linfitthresh'
            eiFitThresh = varargin{j*2};
        case 'showsoma'
            showSoma = varargin{j*2}; 
        case 'figurenum'
            figureNum = varargin{j*2}; 
        case 'plotcoords'
            plotCoords = varargin{j*2};
        case 'ei_thresh'
            ei_thresh = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end


if size(electrodeCoords,1) > size(electrodeCoords,2)
    electrodeCoords = electrodeCoords'; 
end
xc = electrodeCoords(1,:);
yc = electrodeCoords(2,:);

eiAmpsT = eiAmps'; %Threshold eiAmps for contour plotting

eiAmpsT(find(eiAmpsT>max(eiAmpsT)/threshForContours)) = max(eiAmpsT)/threshForContours;
if figureNum==0
    f=figure();
    figureNum=f.Number;
else figure(figureNum)
end
imshow(imageData); axis image; hold on;
% Generate contour plot in a new figure window
hex_contour(xc, yc, eiAmpsT, numContourLevels, 'fig_or_axes', figureNum, ...
    'contourSpacing', 'linear','plotCoords',plotCoords,'cmap','summer');
%scatter(xc,yc,eiAmpsT*7,'filled','MarkerFaceColor', 'black');
% Calculate linear fit using only eiAmpls above a given threshold
if size(eiAmps,2) == 1
    eiAmps = eiAmps'; 
end
[~,col,~] = find(eiAmps > eiFitThresh);
aa = round(eiAmps(col))';
yy = xc(col)';
xx = yc(col)';
 
% for plotting the soma location based on the largest recorded amplitudes
sortaa = sort(aa,1,'descend');
largestAmps = sortaa(1:2);
[~,IA,~] = intersect(aa,largestAmps);

COMx = 1/sum(largestAmps) * sum(xx(IA).*aa(IA));
COMy = 1/sum(largestAmps) * sum(yy(IA).*aa(IA));


[xf, yf] = weighted_axon_poly_reg(eiAmps, 'ei_thresh', ei_thresh,'newXYCoords',electrodeCoords');
hold on; plot(xf,yf,'-','Color','black','LineWidth',2);
if showSoma
    hold on; scatter(COMy,COMx,8*mean(largestAmps), 'black','filled');
end
axis image; axis off; 

if printElecs
    for e = 1:512;
        text(xc(e),yc(e),num2str(e),'HorizontalAlignment','center')
    end
end

if saveFiles
    [savename,savepath] = uiputfile('*.*','save images as');
    savingName = [savepath savename];
    saveas(gcf, savingName,'epsc');
    saveas(gcf, savingName,'jpg');
end
end