% **********************************************************************
% I=zernike_rec(Z,SZ,ZBFSTR)
% I=zernike_rec(Z,SZ,ZBFSTR,OPTSTARTIND)
%
% Reconstruct an image from some of its Zernike moments. We suppose that
% some amount of first (lower) order moments, not necessarily all moments,
% are available.
%
% Z ............. The COMPLEX Zernike moments
% SZ ............ The size of the square image
% ZBFSTR ........ The Zernike basis functions, must be precomputed 
%                 with zernike_bf()
% OPTSTARTIND ... If given, the starting index from which reconstr. is
%                 performed. Allows to ignore a higher or lower number
%                 of lower order polynomials. THIS IS NOT THE ORDER, BUT
%                 ITS INDEX!
%                 Default value: 4 for a full basis,
%                                3 for a basis w/o negative repetitions
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
% 20.12.2009 cw: -Begin
% **********************************************************************

function I=zernike_rec(Z,SZ,ZBFSTR,OPTSTARTIND)  
    len=length(Z);    
    od=ZBFSTR.orders;    
    id=ZBFSTR.index;
    bf=ZBFSTR.bf;   
    maxorder=ZBFSTR.maxorder;
    
    if nargin<4
        if ZBFSTR.withneg
            OPTSTARTIND=4;
        else
            OPTSTARTIND=3;
        end
        
    end
                
    if size(bf,3)~=len
        size(bf,3)
        len %#ok<NOPRT>
        error (['**** ERROR *** in zernike_rec: Zernike basis ' ...
            'functions do not match input vector!' ]);
    end    
    
    I=zeros(SZ);
    % ---- The reconstructed image is a weighted sum of the 
    % ---- Zernike basis functions
    % ---- We discard the first components
    
    % ---- If the set of basis functions is "full", that means it also
    % ---- contains the redundant functions with negative repetitions,
    % ---- then the reconstruction is a simple sum.    
    if ZBFSTR.withneg
        for i=OPTSTARTIND:len                        
            I=I+Z(i).*bf(:,:,i);
        end
        
    % ---- If the set of basis functions is _NOT_ "full", we also need
    % ---- to artificially include the missing redundant values
    else
        for i=OPTSTARTIND:len              
            I=I+Z(i).*bf(:,:,i);            
            
            % ---- Recall that Ap-q = A*pq
            p=od(i,1);
            q=od(i,2);
            if q~=0
                ieq=id(1+maxorder+p,1+maxorder+abs(q));
                if ieq<1
                    error ('Invalid equivalent moment!');
                end

                I=I+(conj(Z(ieq)).*conj(bf(:,:,ieq)));
            end
%             fprintf ('rec it. %d: min=%f,max=%f.\n',i,min(min(I)),max(max(I)));
        end                
    end
    % ---- Take the real part
    I=real(I);
 
    % ---- threshold
    I=im2bw(I,0.5);
    