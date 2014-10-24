function fit = getSTC(data,rnd)

%-------------------------------------------
%
% getSTC(data)
%
% do a simple STC analysis on the spiking
% response to a stimulus. assumes that stimuli
% have been temporally filtered so that we 
% only need to compute the filters at 
% one time point.
%
%
% inputs:
% data -- stucture containing the data
% data.X_ct -- stimulus
% data.R_t -- response
%
% optional inputs:
% rnd -- number of randomizations to use (default = 0)
%
% outputs:
% fit -- structure with the following fields
% fit.k -- sta, normalized to be unit 1
% fit.e -- eigenvalues of spike triggered ensemble
% fit.V -- eigenvectors of spike triggered ensemble
%
% freeman, 4-1-2012
%-------------------------------------------

if ~exist('rnd','var') || isempty(rnd)
    rnd = 0;
end

% compute the mean of the spike triggered ensemble
tmp = data.X_ct(:,data.R_t>0);
tmp = bsxfun(@times,tmp,data.R_t(data.R_t>0));
sta = mean(tmp,2);
%sta = mean(data.X_ct(:,data.R_t>0),2);
staUnit = sta/norm(sta);

% project out the sta from the stimulus ensemble
tmpStimRmv = zeros(size(data.X_ct));
for iStim=1:size(data.X_ct,2)
    stimFrame = data.X_ct(:,iStim);
    tmpStimRmv(:,iStim) = stimFrame - (stimFrame'*staUnit)*staUnit;
end

% compute eigenvectors of spike triggered ensemble
tmpSpkRmv = tmpStimRmv(:,(data.R_t)>0);
tmpSpkRmv = bsxfun(@times,tmpSpkRmv,data.R_t(data.R_t>0));
[U S V] = svd(cov(tmpSpkRmv'));

% collect the results
fit.k = staUnit;
fit.e = diag(S(1:end-1,1:end-1));
fit.V = U(:,1:end-1);

% stuff for a randomization test, need to implement...
if rnd == 1
    %fprintf('(getSTC) doing randomization test\r');
    eigValsRand = zeros(rnd,size(data.X_ct,1)-1);
    for i=1:100
        shiftAmt = ceil(rand*length(data.R_t));
        spkShift = circshift(data.R_t,[0 shiftAmt]);
        tmpSpkRmv = tmpStimRmv(:,spkShift>0);
        tmpSpkRmv = bsxfun(@times,tmpSpkRmv,spkShift(spkShift>0));
        [U S V] = svd(cov(tmpSpkRmv'));
        eigValsTmp = diag(S);
        eigValsRand(i,:) = eigValsTmp(1:end-1);
    end
    fit.e_rnd = eigValsRand;
end

% nested hypothesis test
if rnd == 2
    %fprintf('(getSTC) doing nested randomization test\r');
    %%
    flag = 0;
    cur_eig = 1;
    while flag == 0;
        eigValsRand = zeros(rnd,size(data.X_ct,1)-1);
        for i=1:100
            shiftAmt = ceil(rand*length(data.R_t));
            spkShift = circshift(data.R_t,[0 shiftAmt]);
            tmpSpkRmv = tmpStimRmv(:,spkShift>0);
            tmpSpkRmv = bsxfun(@times,tmpSpkRmv,spkShift(spkShift>0));
            [U S V] = svd(cov(tmpSpkRmv'));
            eigValsTmp = diag(S);
            eigValsRand(i,:) = eigValsTmp(1:end-1);
        end
        %%
        if fit.e(cur_eig) > prctile(eigValsRand(:,cur_eig),99)
            % remove this eigenvector
            v = fit.V(:,cur_eig);
            for iStim=1:size(tmpStimRmv,2)
                stimFrame = tmpStimRmv(:,iStim);
                tmpStimRmv(:,iStim) = stimFrame - (stimFrame'*v)*v;
            end
            cur_eig = cur_eig + 1
        else
            flag = 1;
        end
    end
    fit.n_sig = cur_eig;
end