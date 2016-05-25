% makeaxis(dt,nbins,xstart): Make an evenly spaced axis axis vector of length nbins with
% spacing dt and starting at xstart

function x = makeaxis(dt,nbins,xstart)

if (nargin < 3)
    xstart = 0;
end
x = xstart:dt:(xstart+(nbins-1)*dt);
