function [error dataLabels] = linFitSimulError(modelParams, labels, data, plotMod, varargin)

%objective funciton used by linFitSimulMinError


p = inputParser;



%threshold mean squared error (squared distance from point to fit line)
%before additional "regions" are added to account for nonlinearities
p.addParamValue('plotAxes', [], @iscell) %cell array of handles corresponding to indeces in plotMod 

p.parse(varargin{:})

plotAxes = p.Results.plotAxes;



%initialize structures for storing data and model values
for ii = 1:length(data)
    %model intercepts
    model(ii).tp = 0; %intercept of model through primary electrode axis
    model(ii).tsPos = inf; %intercept of model through positive secondary electrode axis (inf signifies no measured sec-alone thresh)
    model(ii).tsNeg = inf; %intercept of model through positive secondary electrode axis (inf signifies no measured sec-alone thresh)

    %lambda values: 1 for each region (primary, secondary positive,
    %secondary negative
    model(ii).pLam = 0;
    model(ii).sLamPos = inf;
    model(ii).sLamNeg = inf;
end


%unpack parameters
for ii = 1:length(modelParams)
    switch labels(ii).fieldName   
        case 'tp'
            for jj = 1:length(model)
                model(jj).tp = modelParams(ii); %same for all secondary electrodes
            end
        case 'tsPos'
            model(labels(ii).sInd).tsPos   = modelParams(ii);
        case 'tsNeg'
            model(labels(ii).sInd).tsNeg   = modelParams(ii);
        case 'pLam'
            model(labels(ii).sInd).pLam    = modelParams(ii);
        case 'sLamPos'
            model(labels(ii).sInd).sLamPos = modelParams(ii);
        case 'sLamNeg'
            model(labels(ii).sInd).sLamNeg = modelParams(ii);
        otherwise
            error('unrecognized field name')
    end
end


%calculate error for current parameter set

pairErrs = zeros(1,length(data));
plotColors = lines(3);
for ii = 1:length(data)
    q1 = inf*ones(2,3); q2 = inf*ones(2,3);
    q3 = inf*ones(2,3); %used for plotting only
    
    r = cell(1,3);

    p = [data(ii).sVal'; data(ii).pVal'];
    dists = inf*ones(3,length(data(ii).pVal)); %row for each potential model line

    %positive primary region: 2 points on line, [sVal; pVal]
    q1(:,1) = [0; model(ii).tp];
    q2(:,1) = [1; model(ii).tp - model(ii).pLam];
    q3(:,1) = [-1; model(ii).tp + model(ii).pLam]; %for plotting only
    
    r{1} = zeros(size(p));
    
    for jj = 1:size(p,2)
        %min distance between each point and line
        dists(1,jj) = abs(det([q2(:,1)-q1(:,1), p(:,jj)-q1(:,1)]))/norm(q2(:,1)-q1(:,1));
        
        if any(plotMod == ii)
            %alternative: find nearest point on line, then calculate distance
            r{1}(:,jj) = (dot(p(:,jj)-q2(:,1), q1(:,1)-q2(:,1))*q1(:,1)+dot(p(:,jj)-q1(:,1), q2(:,1)-q1(:,1))*q2(:,1))/dot(q2(:,1)-q1(:,1), q2(:,1)-q1(:,1));
            if dists(1,jj) - norm(r{1}(:,jj)-p(:,jj)) > 1e-10
                error('difference between distance measures')
            end
        end
    end
    
    %positive secondary region: 2 points on line
    if ~isinf(model(ii).tsPos)
        q1(:,2) = [model(ii).tsPos;                     0];
        q2(:,2) = [model(ii).tsPos - model(ii).sLamPos; 1];
        q3(:,2) = [model(ii).tsPos + model(ii).sLamPos; -1]; %for plotting only
        
        r{2} = zeros(size(p));

        for jj = 1:size(p,2)
            %min distance between each point and line
            dists(2,jj) = abs(det([q2(:,2)-q1(:,2), p(:,jj)-q1(:,2)]))/norm(q2(:,2)-q1(:,2));
            
            if any(plotMod == ii)
                %alternative: find nearest point on line, then calculate distance
                r{2}(:,jj) = (dot(p(:,jj)-q2(:,2), q1(:,2)-q2(:,2))*q1(:,2)+dot(p(:,jj)-q1(:,2), q2(:,2)-q1(:,2))*q2(:,2))/dot(q2(:,2)-q1(:,2), q2(:,2)-q1(:,2));
                if dists(2,jj) - norm(r{2}(:,jj)-p(:,jj)) > 1e-10
                    error('difference between distance measures')
                end
            end
        end
    end
    
    %negative secondary region: 2 points on line
    if ~isinf(model(ii).tsNeg)
        q1(:,3) = [-model(ii).tsNeg                    ; 0];
        q2(:,3) = [-model(ii).tsNeg - model(ii).sLamNeg; 1];
        q3(:,3) = [-model(ii).tsNeg + model(ii).sLamNeg; -1]; %for plotting only

        r{3} = zeros(size(p));

        for jj = 1:size(p,2)
            %min distance between each point and line
            dists(3,jj) = abs(det([q2(:,3)-q1(:,3), p(:,jj)-q1(:,3)]))/norm(q2(:,3)-q1(:,3));
            
            if any(plotMod == ii)
                %alternative: find nearest point on line, then calculate distance
                r{3}(:,jj) = (dot(p(:,jj)-q2(:,3), q1(:,3)-q2(:,3))*q1(:,3)+dot(p(:,jj)-q1(:,3), q2(:,3)-q1(:,3))*q2(:,3))/dot(q2(:,3)-q1(:,3), q2(:,3)-q1(:,3));
                if dists(3,jj) - norm(r{3}(:,jj)-p(:,jj)) > 1e-10
                    error('difference between distance measures')
                end
            end
        end
    end
    
%     
%     %error is sum of squared minimum distances, with the contraint that
%     %each region contains a contiguous set of current ratios (default to
%     inclusion in primary-positive region if not contiguous)
    [~, ratioOrder] = sort(p(1,:)./p(2,:));
    
    [~, iMin] = min(dists,[],1);
    
    % constrain regions to continuous set of current ratios
    firstReg1 = find(iMin(ratioOrder)~=3, 1, 'first');
    lastReg1 =  find(iMin(ratioOrder)~=2, 1, 'last');
    
    iReg = iMin;
    reg1DefaultInds = ratioOrder(firstReg1:lastReg1);
    iReg(reg1DefaultInds) = 1;
    
    regDist = zeros(1,size(dists,2));
    for jj = 1:size(dists,2)
        regDist(jj) = dists(iReg(jj),jj);
    end
    
    dataLabels(ii).p    = iReg == 1;
    dataLabels(ii).sPos = iReg == 2;
    dataLabels(ii).sNeg = iReg == 3;
    
    pairErrs(ii) = sum(regDist.^2);
        
    if any(plotMod == ii)
        if isempty(plotAxes)
            figure; hold on
        else
            axes(plotAxes{plotMod == ii}); cla; hold on
        end
        for jj = 1:3
            plot(p(1,iReg==jj), p(2,iReg==jj), 'o', 'markerEdgeColor', plotColors(jj,:))
            for kk = 1:size(p,2)
                if iReg(kk)==jj
                    plot([p(1,kk) r{jj}(1,kk)], [p(2,kk) r{jj}(2,kk)], '-', 'color', plotColors(jj,:))
                    plot(p(1,kk), p(2,kk), 'o', 'markerFaceColor', plotColors(jj,:), 'markersize', 100*(regDist(kk)^2)+0.1)
                end
            end            
            if ~isinf(q1(1,jj))
                plot([q3(1,jj) q2(1,jj)], [q3(2,jj) q2(2,jj)], '-', 'color', plotColors(jj,:))
            end
        end
        plot([0 0], [min(p(2,:)) max(p(2,:))], 'k-')
        plot([min(p(1,:)) max(p(1,:))], [0 0], 'k-')
        
%         for jj = 1:ii
%             text(0, (jj-1)/10, num2str(pairErrs(jj)))
%         end
        
        axis equal
        set(gca, 'xlim', [-0.1 0.1]+[min(p(1,:)) max(p(1,:))], 'ylim', [-0.1 0.1]+[min(p(2,:)) max(p(2,:))])
        %title(['sElec ' num2str(ii) ', error = ' num2str(sum(regDist.^2))])
        title(['sElec ' num2str(ii) ', error = ' num2str(pairErrs(ii))])
        %pause(0.5)
    end
end
    
error = sum(pairErrs);

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    