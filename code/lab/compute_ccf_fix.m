function [ccf, time]=compute_ccf_fix(sp1, sp2, options)
%
% wrapper for compute_ccf that makes off-centering artifact teeny-tiny
%
% usage:  [ccf, time]=compute_ccf_fix(sp1, sp2, options)
%
%    arguments and outputs of function are the same as for compute_ccf
%
% originally written by MG, edited by LJ



%% code taken directly from compute_ccf to determine parameter values this wrapper function requires
% (if compute_ccf has changed since 2010-09, this may also need to be updated!)

% default number of samples in CCF
N = 63;

% default arguments
default_options = struct('offset', 100e-3, 'shuffle', 'none', 'trial', 1);
if (nargin == 2)
  options = default_options;
end

% check if individual items are not specified
if ~isfield(options,'offset')
  options.offset = default_options.offset;
end
if ~isfield(options,'dt')
  options.dt = options.offset / N;
end

%% 

dt=options.dt;

bin=dt/.00005;
options.dt=.00005;

tccf=compute_ccf(sp1, sp2, options);

ccf=[];
m=ceil(length(tccf)/2);
for i=m+bin/2+1:bin:length(tccf)-bin
    ccf=[ccf mean(tccf(ceil(i:i+bin)))];
end

ccf=[mean(tccf(ceil(m-bin/2:m+bin/2))) ccf];

for i=m-bin/2+1:-bin:bin+1
    ccf=[mean(tccf(ceil(i-bin:i))) ccf];
end

time=[dt:dt:options.offset-dt];
time=sort([-time 0 time]);

end


