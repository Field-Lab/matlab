% **********************************************************************
% Z=zernike_mom(I,ZBFSTR)
% 
% Calculate the Zernike moments from an image. We suppose that
% some amount of first (lower) order moments, not necessarily all moments,
% are available.
%
% Z ........ The COMPLEX Zernike moments 
%
% I ........ The input image 
% ZBFSTR ... The Zernike basis functions, can be computed by zernike_bf
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

function Z=zernike_mom(I,ZBFSTR)    

    if size(I,1)~=size(I,2)
        error ('The image must be of square size!');
    end   

    bf=ZBFSTR.bf;
    pq=ZBFSTR.orders;
    id=ZBFSTR.index;
            
    % ---- Iterate through the moments, beginning from the lowest
    % ---- to the highest    
    len=size(bf,3);
    Z=zeros(len,1);    
    for flat=1:len      
        m=pq(flat,1);
        n=pq(flat,2);
        Z(flat) = (m+1)/pi*...
            sum(sum(I.*conj(bf(:,:,flat))));                    
    end     
    
