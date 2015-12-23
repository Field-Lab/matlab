function total_sum = cacColSum()
%cacColSum() - example parallel job worker function.  See examples for
%usage.
%
% See also cacsched, LittleJohn
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

fprintf('I am %d', labindex);

if labindex == 1
    A = labBroadcast(1,magic(numlabs));
else
    A = labBroadcast(1);
end

%pause(600);
colsum = sum(A(:,labindex));

total_sum = gplus(colsum);

