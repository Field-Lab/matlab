% **********************************************************************
% PQ=zernike_orderlist(ORDER)
% PQ=zernike_orderlist(ORDER, WITHNEG)
%
% Create the moment indices (oder p and repetition) satisfying the 
% necessary criteria (p-|q|=even, |q|<p) up to order p=ORDER.
%
% If WITHNEG is 1, then the basis functions with negative repetition are
% included.
%
% PQ is a matrix with 2 columns (p and q).
%
% Christian Wolf, http://liris.chrs.fr/christian.wolf
%
% Permission is granted for anyone to copy, use, modify, or distribute this
% program and accompanying programs and documents for any purpose, provided
% this copyright notice is retained and prominently displayed, along with
% a note saying that the original programs are available from my web page.
%
% The programs and documents are distributed without any warranty, express 
% or implied.  As the programs were written for research purposes only, they 
% have not been tested to the degree that would be advisable in any 
% important application.  
% All use of these programs is entirely at the user's own risk.
%
% **********************************************************************

% **********************************************************************
% Change log:
% 21.12.2009 cw: -Begin
% **********************************************************************

function PQ=zernike_orderlist(ORDER, WITHNEG)

    if nargin<2
        WITHNEG=0;
    end

    PQ=[];    
    if WITHNEG    
        for p=0:ORDER
            for q=-p:p                  
                if mod(abs(p-q),2)==0
                    PQ = [ PQ ; p q ];
                end                    
            end
        end
    else
        for p=0:ORDER
            for q=0:p                  
                if mod(abs(p-q),2)==0
                    PQ = [ PQ ; p q ];
                end                    
            end
        end
    end