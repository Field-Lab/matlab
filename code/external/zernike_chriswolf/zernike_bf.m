% **********************************************************************
% ZBFSTR=zernike_bf(SZ,ORDER)
% ZBFSTR=zernike_bf(SZ,ORDER,WITHNEG)
% 
% Calculate the Zernike basis functions up to order ORDER for a
% square image of size SZ (rows = cols!)
% If WITHNEG is 1, then the basis functions with negative repetition are
% included.
%
% The array contains basis function entries for couples of (p,q), where q
% can be negative. 
%
% ZBFSTR is a structure with the following elements:
%  .bf ....... the basis functions (SZ, SZ, LENGTH).
%  .orders ... a matrix (2, LENGTH), where the first column holds the 
%              order p and the second column the repetition q.
%  .index .... a matrix (1+2*ORDER, 1+2*ORDER) which holds for each element
%              the index in .orders and in .bf of this couple of p,q.
%              -1 means that there are no basis functions for this couple
%              ATTENTION!!! since indixes in matlab always start with 1 
%              they should be negative (as well as positive) here,
%              the index starts with -ORDER, which requires adding ORDER+1
%              to both dimensions. So given an order of 12,
%              instead of ZBFSTR.index(-5,3) we write 
%              ZBFSTR.index(8,16);
%  .maxorder . the maximum order, basically the ORDER parameter given
%              to this function.
%  .withneg .. 1 if there are basis functions with negative repetition,
%              so basically the WITHNEG parameter to this function
%
% All pixels out of the disk inscribed into the image square are ignored.
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

function ZBFSTR=zernike_bf(SZ,ORDER,WITHNEG)

    if nargin<3
        WITHNEG=0;
    end

    limitfastcomp=50;
    
    % ---- Precompute factorials 
    % ---- Indexes start with 0!!!!!!!!
    F=factorial(0:ORDER);       

    % ---- Create the moment indices from the flat indices in the 
    % ---- flat moments array
    pq=zernike_orderlist(ORDER, WITHNEG);
    
    len=size(pq,1);
    szh=SZ/2;
    
    % ---- Create the index matrix
    pqind=-1*ones(1+2*ORDER,1+2*ORDER);            
    src=1+ORDER+pq;    
    pqind(sub2ind(size(pqind),src(:,1),src(:,2)))=(1:len)';
       
    % ---- Precalculate the polynomials which do not depend on 
    % ---- the position of the pixel
    % ---- Indexes start with 
    % ---- DIM 1 (=m): 1+ORDER
    % ---- DIM 2 (=n): 1+ORDER
    % ---- DIM 3 (=s): 1
    %
    % ---- The fast but not robust calculation for the lower order
    % ---- polynomials
    Rmns=zeros(1+2*ORDER,1+2*ORDER+1,1+2*ORDER);        
    for flat=1:min(len,limitfastcomp);               
        m=pq(flat,1);
        n=pq(flat,2);      
        mpnh=floor((m+abs(n))/2);
        mmnh=floor((m-abs(n))/2);
        for s=0:mmnh
            Rmns(1+ORDER+m,1+ORDER+n,1+s)=(((-1)^s)*F(1+m-s))/...
                (F(1+s)*F(1+mpnh-s)*F(1+mmnh-s));
        end
    end
    
    % ---- The higher order polynomials are computed more slowly, 
    % ---- but more robustly; else the factorials get too high and
    % ---- and create numerical instabilities    
    for flat=limitfastcomp+1:len
        m=pq(flat,1);
        n=pq(flat,2);
        mpnh=floor((m+abs(n))/2);
        mmnh=floor((m-abs(n))/2);
        for s=0:mmnh            
            Rmns(1+ORDER+m,1+ORDER+n,1+s)=((-1)^s)*robust_fact_quot(...
                1:m-s,...
                [ 1:s 1:mpnh-s 1:mmnh-s ]);              
        end
    end    
    
    % ---- Calculate the basis functions
    ZBF=zeros(SZ,SZ,len);
    for y=1:SZ
        
%         fprintf ('[%3.3d]',y);  

        for x=1:SZ                                           
            % ---- Take cartesian coordinates with the origin in the
            % ---- middle of the image and transform into polar
            rho=sqrt((szh-x)^2+(szh-y)^2);
            theta=atan2(szh-y,szh-x);
            if rho>szh
                continue
            end
            rho=rho/szh;
            if theta<0
               theta=theta+2*pi;
            end
            
            % ---- Iterate through the different orders and repetitions
            for flat=1:len                                                           
                m=pq(flat,1);
                n=pq(flat,2);

                R=0;               
                for s=0:(m-abs(n))/2                               
                    R=R+Rmns(1+ORDER+m,1+ORDER+n,1+s)*(rho^(m-2*s));                    
                end

                ZBF(y,x,flat) = R*exp(n*theta*1j);
            end                
        end
    end
    
    % ---- Package the whole thing into a structure
    ZBFSTR=struct;
    ZBFSTR.maxorder=ORDER;
    ZBFSTR.withneg=WITHNEG;
    ZBFSTR.orders=pq;
    ZBFSTR.index=pqind;
    ZBFSTR.bf=ZBF;

