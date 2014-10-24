function plotSubunitWeights(dat,fit,newFig,rows)

%-------------------------------------------
%
% plotSubunitWeights(fit,newFig,rwos)
%
% plot the subunit weight parameters of an LN cascade model
%
% fit must contain the fields:
% 
% fit.A_sc -- weights from stimulus to subunits, S x C
% or
% fit.I_sc -- binary assignments from stimulus to subunits
% fit.S -- number of subunits
% fit.B_s -- weights from subunits to cell, 1 x S (optional)
%
% additional optional inputs:
% newFig -- 1 to plot in a new figure
% rows -- integer indicating a row to plot into
%         (if we want to plot into an existing figure)
% 
% freeman, 3-19-2012
%-------------------------------------------


% need to get center and dist based on locs...
center = mean(fit.locs_c);
dist = max(getDistances(fit.locs_c,mean(fit.locs_c)));
offset = 1;
mkSz = 75;

if newFig
    figure
    if exist('rows','var')
        set(gcf,'Position',[0 237 1440 180*rows(1)]);
    else
        set(gcf,'Position',[0 237 1440 180]);
    end
    set(gcf,'Color','white');
    colormap(redgrayblue);
end

if isfield(fit,'B_s');
    [~,order] = sort(fit.B_s,'descend');
else
    order = 1:fit.S;
end

if isfield(fit,'A_sc')
    weightsAll = fit.A_sc;
else
    weightsAll = fit.I_sc;
end

nonSingInds = sum(fit.I_sc,2) ~= 1;
nNonSing = sum(nonSingInds);
subClrs = zeros(fit.S,3);
subClrs(nonSingInds,:) = cmap_hsv_dark(nNonSing);


for is=1:fit.S
    
    if exist('rows','var')
        subp(rows(1),fit.S,(rows(2)-1)*fit.S+is,0.01);
    else
        subp(1,fit.S,is,0.01);
    end
    
    % use the same color for each subunit as from above,
    % but adjust the saturation of each cone's
    % color according to the weight
    weights = weightsAll(is,:);
    clrsThisSub = repmat(subClrs(is,:),[length(weights) 1]);
    clrsThisSub_hsv = rgb2hsv(clrsThisSub);
    clrsThisSub_hsv(:,2) = clrsThisSub_hsv(:,2).*weights(:)/max(weights(:));
    clrsThisSub_hsv(:,3) = clrsThisSub_hsv(:,3)+0.2*(1-weights(:)/max(weights(:)));
    clrsThisSub = hsv2rgb(clrsThisSub_hsv);
    clrsThisSub(weights(:)==0,1) = 1;
    clrsThisSub(weights(:)==0,2) = 1;
    clrsThisSub(weights(:)==0,3) = 1;
    % plot the individual cones
    % (use this instead of scatter because scatter
    % fucks stuff up in the rest of the plot)
    hold on
    for ic=1:length(weights)
        h = plot(dat.locs_c(ic,1),dat.locs_c(ic,2),'o');
        set(h,'MarkerFaceColor',clrsThisSub(ic,:));
        set(h,'MarkerSize',12);
        set(h,'MarkerEdgeColor','k');
    end
    hold off
    % wrangle the axes
    xlim(dat.xRange+[-2 2]);
    ylim(dat.yRange+[-2 2]);
    axis equal off;
    set(gca,'CLim',[-max(weights(:)) max(weights(:))]);
    set(gca,'YDir','reverse');
end