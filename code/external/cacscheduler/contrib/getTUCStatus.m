function [numNodes, numJobs] = getTUCStatus(verbose)

if nargin < 1
    verbose = 1;
end

s = urlread('http://www.cac.cornell.edu/matlab/status/idlecores.aspx');
fprintf(s);
fprintf('\n');