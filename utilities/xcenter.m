function centered = xcenter(traces, centers, nleft, nright, exval)
% XCENTER    Recenter traces around given points
% usage: centered = xcenter(traces, centers, nleft, [nright, exval])
%
% input:  traces   MxN matrix of traces as column vectors
%         centers  N length vector of centering points, one for each trace
%         nleft    Number of left points in the recentered traces
%         nright   Number of right points, defaults to equal nleft
%         exval    Value to use as extrapolant when centering goes beyond
%                  end of trace data
%
% output: (nleft+nright+1)xN matrix of recentered traces as column vectors
%
% 2010-06 phli
%


error(nargchk(3, 5, nargin));

if nargin < 5
    exval = 0;
end

if nargin < 4
    nright = nleft;
end

[tracelen,numtraces] = size(traces);
centeredlen = nleft + nright + 1;

% ToDo
% This could be made more vectorized, e.g. pad whole TRACES matrix with
% necessary EXVALS, use SUB2IND and RESHAPE to pull all values at once.
centered = ones(centeredlen, numtraces);
centered = centered .* exval;
for i = 1:size(traces, 2)
    center = centers(i);
    ileft  = center - nleft;
    iright = center + nright;
    
    % Fix overflows on either end
    if ileft < 1
        padleft = abs(ileft) + 1 + 1;
        ileft = 1;
    else
        padleft = 1;
    end
    
    if iright > tracelen
        padright = centeredlen - abs(tracelen - iright);
        iright = tracelen;
    else
        padright = centeredlen;
    end
    
    centered(padleft:padright,i) = traces(ileft:iright, i);
end