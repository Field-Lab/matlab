function [v resp] = getPairwiseIntAll(train,test,fit,crossVal)

nCones = size(train.X_ct,1);
v = zeros(nCones,nCones);

for iCone=1:nCones
    for jCone=iCone+1:nCones
        %fprintf('(clusterCones) cone pair (%g,%g)/(%g,%g)\r',iCone,jCone,nCones,nCones);
        out = getPairwiseInt(train,test,iCone,jCone,fit,0);
        if exist('crossVal','var') && crossVal == 1
            testStat = out.alt_err_test - out.null_err_test;
        else
            testStat = out.alt_err - out.null_err;
        end
        v(iCone,jCone) = testStat;
        %fprintf('(clusterCones) test stat is %g\r',testStat);
    end
end
