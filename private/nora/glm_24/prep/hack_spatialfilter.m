
function spfilter = hack_spatialfilter(stringdate, cellname,xdim,ydim)

eval(sprintf('load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/%s/BW_mapPRJ/STAandROI_default/STAandROI_%s', ...
    stringdate,cellname));

klen = length(xdim);
%load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/2012-08-09-3/BW_mapPRJ/STAandROI_default/STAandROI_ONPar_1276.mat
%load /netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/2012-08-09-3/BW_mapPRJ/STAandROI_default/STAandROI_OFFPar_1471.mat
STA = ( STAandROI.STA(xdim,ydim,:) );
%imagesc(squeeze(STAandROI.STA(:,:,5)));
STA = reshape(STA, [klen^2,30])  - mean(STA(:)) ;


isfiniteSTA = isfinite(STA);

if isempty(find(isfiniteSTA == 0 ) )


    [U,S,V]  = svd (STA);
    S = diag(S);


    %imagesc(reshape(U(:,1), [80,40]))

    xx2 = ( S(1)*U(:,1)*V(5,1) + S(2)* U(:,2)*V(5,2) ) / norm( S(1)*U(:,1)*V(5,1) + S(2)* U(:,2)*V(5,2)) ;
    xx = ( S(1)*U(:,1)*V(5,1) ) / norm( S(1)*U(:,1)*V(5,1) ) ;
    xxclean = xx;

    c1 = std(xx(:));
    throw = find(abs(xx) < 1*c1);


    xxclean(throw) = 0;
    xxclean = xxclean / norm(xxclean);


    figure; 
    subplot(1,2,1); imagesc(squeeze(STAandROI.STA(:,:,5))); colorbar
    subplot(1,2,2); imagesc(reshape(xx,[klen,klen] ) ); colorbar;
    %subplot(1,3,3); imagesc(reshape(xxclean,[80,40]) ); colorbar;

    spfilter = xx;
else
    spfilter = zeros(klen^2 , 1);
end



end