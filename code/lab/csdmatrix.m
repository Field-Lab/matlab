function csdM = csdmatrix(positions, varargin)
% CSDMATRIX     Generate EI->CSD transform matrix: csdM * EI = CSD
% usage: csdM = csdmatrix(positions, opts)
%
% opts: badelectrodes   []      Electrodes to ignore for calculation
%
% Bad electrodes are not used to calculate values for neighbors, so their
% columns are zeroed.  
%
% Bad electrodes and convex hull electrodes are simply set to zero in the
% CSD output, so their rows are zeros.  This should be corrected post-hoc
% by setting the CSD output value to NaN if having spurious zero values
% will be confusing; the zeros are left here for sparse matrix
% multiplication efficiency.
%
% Convex hull electrodes are allowed to contribute to neighboring electrode
% values, so their columns are not zeroed.
%
% Be careful that you don't have floating point inaccuracies in the point
% coordinates or you may get funny behavior around the convex hull.
%
% 2013-05, phli
%

opts = inputParser();
opts.addParamValue('badelectrodes', []);
opts.parse(varargin{:});
opts = opts.Results;

if islogical(opts.badelectrodes), opts.badelectrodes = find(opts.badelectrodes); end


dt = DelaunayTri(positions);

% Get convex hull before cutting bad electrodes
ch = dt.convexHull;

% Cut bad electrodes from triangulation
dt.X(opts.badelectrodes,:) = [];

% Create basic csdM crossterms
dtL = dtlengths(dt);
vdL = vdlengths(dt);
csdM = vdL ./ dtL;

% Insert zero rows/cols for bad electrodes
if ~isempty(opts.badelectrodes)
    csdM = insertrows(csdM,  zeros(length(opts.badelectrodes), size(csdM,2)), opts.badelectrodes - (1:length(opts.badelectrodes))');
    csdM = insertrows(csdM', zeros(length(opts.badelectrodes), size(csdM,1)), opts.badelectrodes - (1:length(opts.badelectrodes))')';
end

% Add self terms
sums = sum(csdM,2);
csdM = diag(sums) - csdM;

% Electrodes on the convex hull shouldn't be included as they don't have
% complete Voronoi polygons, but they should be allowed to factor into
% neighboring values, so only clear their rows not columns.
csdM(ch,:) = 0;
