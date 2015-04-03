% AKHeitman 2014-0428
% Assumes STA is given in (xcoord, ycoord, frames)
function [leadspfilter leadtimefilter] = spatialfilterfromSTA(STA,xcoord,ycoord)


% Reshape into space,time  2d notation
duration = size(STA,3);
klen = length(xcoord);
STA = (STA(xcoord,ycoord,:) );
STA = reshape(STA, [klen^2,duration])  - mean(STA(:)) ;

% Making sure no wiered NAN stuff
isfiniteSTA = isfinite(STA);

% Singular Value Decomposition
if isempty(find(isfiniteSTA == 0 ) )
    [U,S,V]  = svd (STA);
    S = diag(S);
    % Choosing the V(5,1) lets us make the direction consistent
    % OFFs will look OFF , ONs will look ON
    xx = ( S(1)*U(:,1)*V(5,1) ) / norm( S(1)*U(:,1)*V(5,1) ) ;
    leadspfilter = xx;
    
    yy = S(1)*V(:,1) / norm( S(1)*V(:,1) ); 
    leadtimefilter = yy;
else
    error('STA is not well definted')
end



end