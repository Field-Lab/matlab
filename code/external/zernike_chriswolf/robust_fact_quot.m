% **********************************************************************
% Q = robust_fact_quot (X,Y)
% 
% A slow but robust way to calculate a quotient of two values each of
% which is given as a list of factors. First we cancel all common
% factors and then we do the calculation one by one, alternating between
% multiplications and divisions.
%
% X .... The vector of factors in the numerator
% Y .... The vector of factors in the denominator
% The vectors are not necessarily of the same size
%
% Christian Wolf, christian.wolf@liris.cnrs.fr
% **********************************************************************

% **********************************************************************
% Change log:
% 26.12.2009 cw: -Begin
% **********************************************************************

function R = robust_fact_quot (X,Y)

    % ---- Find all unique elements in the first vector.
    % ---- => candidates for removal
    ca=unique(X);
    
    % ---- Iterate through them and remove the elements common
    % ---- to both vectors
    for k=1:length(ca)
        
        % ---- Find the occurrences of this factor in both vectors
        i1=find(X==ca(k));
        i2=find(Y==ca(k));
        m=min(length(i1),length(i2));
        
        % ---- Remove the common factors
        X(i1(1:m)) = [];
        Y(i2(1:m)) = [];        
    end
    
    % ---- Calculate the quotient, robustly, one factor / divisor after 
    % ---- the other
    R=1;
    l1=length(X);
    l2=length(Y);
    for k=1:min(l1,l2)
        R=R*X(k)/Y(k);
    end
    for k=l2+1:l1
        R=R*X(k);
    end
    for k=l1+1:l2
        R=R/Y(k);
    end
end

