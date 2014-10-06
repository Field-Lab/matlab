function val = erfLikelihood(params, data)
%
% used with fminsearch to estimate erf parameters based on maximum
% likelihood in cases where different amplitudes have different numbers of
% trials
%
% arguments:
%   params: parameters a and b to be fit
%
%   data: 2-d array, with data(1,:) stimulus amplitudes, data(2,:) success
%   rates and data(3,:) number of trials
%
% returns: value is the negative of the portion of the log-likelihood that
% is dependent on the parameters
%
% written by LHJ, Sept. 2011

alpha = params(1);
beta = params(2);

a = data(1,:); %amplitudes
x = data(2,:); %proportion of trials that are successes
n = data(3,:); %number of trials

% calculate expected response probabilities at each current amplitude based on
% current parameters
% p = zeros(1, length(a));
% oneMinP = zeros(1, length(a));

p = 0.5 + 0.5*erf(alpha*a + beta);
oneMinP = 0.5 - 0.5*erf(alpha*a + beta); %so that precision near p == 1 is the same as near p == 0

%set values that have been rounded to 0 to the lowest number that can result from the
%erf without having been rounded to 0
p(p==0) = 5.5511e-17;
oneMinP(oneMinP==0) = 5.5511e-17;



% plug expected response probabilities into part of log likelihood function
% that is dependent on them
ll = 0;
for i = 1:length(a)
    ll = ll + n(i)*(x(i)*log(p(i)) + (1-x(i))*log(oneMinP(i)));
end

%invert sign because fminsearch minimizes and we want the log-likelihood to
%be maximized

val = -ll;