function plotPairwiseInt(data,ind1,ind2,stimType,out,fit)

sig1 = data.X_ct(ind1,:)';
sig2 = data.X_ct(ind2,:)';
spks = data.R_t';

% get 2d rate function
switch stimType
    case 'binary'
        xinds = [-1 0 1 inf];
        yinds = [-1 0 1 inf];
        rate = zeros(3);
    case 'ternary'
        xinds = [-1 1 inf];
        yinds = [-1 1 inf];
        rate = zeros(2);
    otherwise
        nBins = 10;
        mxVal = max(vector([sig1(:);sig2(:)]));
        mnVal = min(vector([sig1(:);sig2(:)]));
        xinds = linspace(mnVal,mxVal,nBins+1);
        yinds = linspace(mnVal,mxVal,nBins+1);
        wbinX = xinds(2)-xinds(1);
        wbinY = yinds(2)-yinds(1);
        rate = zeros(nBins,nBins);
end
for ix=1:length(xinds)-1
    for iy=1:length(yinds)-1
        keepinds = logical(sig1>=xinds(ix)) & logical(sig1<xinds(ix+1)) & logical(sig2>=yinds(iy)) & logical(sig2<yinds(iy+1));
        if sum(keepinds)
            rate(iy,ix) = mean(spks(keepinds));
        end
    end
end


if strcmp(stimType,'binary')
    [X Y] = meshgrid(xinds(1:end-1),yinds(1:end-1));
else
    [X Y] = meshgrid(xinds(1:end-1)+wbinX/2,yinds(1:end-1)+wbinY/2);
    wbinX = xinds(2)-xinds(1);
    wbinY = yinds(2)-yinds(1);
end

%[~,~,null] = pairwiseNullErrFun(out.null_params,X(:),Y(:),rand(size(X(:))),fit);
%[~,~,alt] = pairwiseAltErrFun(out.alt_params,X(:),Y(:),rand(size(X(:))),fit);

%null = reshape(null,size(X));
%alt = reshape(alt,size(X));

null = getPairwisePredicNull(out.null_params,X,Y,fit);
alt = getPairwisePredicAlt(out.alt_params,X,Y,fit);

out.null_params
out.alt_params



%%
figure;
set(gcf,'Position',[277    66   243   740]);
colormap(jet);

if strcmp(stimType,'binary') || strcmp(stimType,'ternary');
    xIndVals = xinds(1:end-1);
    yIndVals = yinds(1:end-1);
else
    xIndVals = xinds(1:end-1) + wbinX/2;
    yIndVals = yinds(1:end-1) + wbinY/2;
end
mxValPlot = prctile(rate(:),100);
mnValPlot = prctile(rate(:),0);

subplot(3,1,1); 
contourf(xIndVals,yIndVals,rate,20);
%imagesc(xIndVals,yIndVals,rate);
set(gca,'Clim',[mnValPlot mxValPlot]);
set(gca,'YDir','normal');
axis square;
box off;
set(gca,'TickDir','out');
title(sprintf('test stat is %g',out.alt_err-out.null_err));

subplot(3,1,2); 
contourf(xIndVals,yIndVals,null,10);
%imagesc(xIndVals,yIndVals,null);
set(gca,'Clim',[mnValPlot mxValPlot]);
set(gca,'YDir','normal');
axis square;
box off;
set(gca,'TickDir','out');
title(sprintf('null lik is %g',out.null_err));

subplot(3,1,3); 
contourf(xIndVals,yIndVals,alt,10);
%imagesc(xIndVals,yIndVals,alt);
set(gca,'Clim',[mnValPlot mxValPlot]);
set(gca,'YDir','normal');
axis square;
box off;
set(gca,'TickDir','out');
title(sprintf('alt lik is %g',out.alt_err));

