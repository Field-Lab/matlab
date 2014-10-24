function U = getClusts(method,v,params)

% get rid of negative values
vTmp = v;
vTmp(vTmp<0) = 0;

if params.rand
    inds = logical(triu(ones(size(vTmp)),1));
    vTmp(inds) = shuffle(vTmp(inds));
end

switch method
    case 'apclusters'
        vTmp(vTmp<params.thresh) = 0;
        vTmp = triu(vTmp) + triu(vTmp)';
        isolated = sum(vTmp,1) == 0;
        nIsolated = sum(isolated);
        vSub = vTmp(~isolated,~isolated);
        idx = apcluster(vSub,median(triuVec(vSub,1)));
        exemp = unique(idx);
        U = zeros(size(v,1),length(exemp));
        for i=1:length(exemp)
            inds = idx == exemp(i);
            U(inds,i) = 1;
        end
        if params.res == 1
            remInds = find(sum(U,2)==0);
            iCount = size(U,2)+1;
            for i=1:length(remInds)
                U(remInds(i),iCount) = 1;
                iCount = iCount + 1;
            end
        end
        
    case 'clusters'
        
        vTmp(vTmp<params.thresh) = 0;
        vTmp = triu(vTmp) + triu(vTmp)';
        isolated = sum(vTmp,1) == 0;
        nIsolated = sum(isolated);
        vSub = vTmp(~isolated,~isolated);
        [S C] = graphconncomp(sparse(vSub>0),'directed','false');
        if S > params.n
            warning('(getClusts) more conn comp than requested clusters');
            params.n = S;
        end
        estLabelsSub = spectralClust(vSub,params.n,'affinity');
        estLabels = zeros(1,size(v,1));
        estLabels(~isolated) = estLabelsSub;
        U = zeros(size(v,1),length(unique(estLabelsSub)));
        for i=1:length(unique(estLabelsSub))
            inds = estLabels == i;
            U(inds,i) = 1;
        end
        if params.res == 1
            remInds = find(sum(U,2)==0);
            iCount = size(U,2)+1;
            for i=1:length(remInds)
                U(remInds(i),iCount) = 1;
                iCount = iCount + 1;
            end
        end
        if params.res == 2
            U = [U eye(size(U,1))];
        end
        %%
    case 'cliques'
        vTmp = v;
        if params.rand
            inds = logical(triu(ones(size(vTmp)),1));
            vTmp = rand(size(vTmp));
            vTmp = vTmp + vTmp';
        end
        c = maximalCliques(vTmp>params.thresh);
        U = zeros(size(v,1),length(c));
        for i=1:length(c)
            inds = c{i};
            U(inds,i) = 1;
        end
        if params.res == 1
            remInds = find(sum(U,2)==0);
            iCount = size(U,2)+1;
            for i=1:length(remInds)
                U(remInds(i),iCount) = 1;
                iCount = iCount + 1;
            end
        end
        if params.res == 2
            U = [U eye(size(U,1))];
        end
    case 'subgraphs-e-nearest'
        vTmp = (v+v')>params.thresh;
        [S C] = graphconncomp(sparse(vTmp),'directed','false');
        U = zeros(size(v,1),S);
        for i=1:S
            inds = C == i;
            U(inds,i) = 1;
        end  
    case 'subgraphs-k-nearest'
        vTmp = kNearestNeighborGraph(v+v',params.thresh);
        [S C] = graphconncomp(sparse(vTmp>0),'directed','false');
        U = zeros(size(v,1),S);
        for i=1:S
            inds = C == i;
            U(inds,i) = 1;
        end  
end

