tic; A = zeros(10000); toc
% Elapsed time is 0.425886 seconds.

tic; A = zeros(10000); toc
% Elapsed time is 0.396264 seconds.
tic; B = A; toc
% Elapsed time is 0.000018 seconds.

clear B;
tic; A(1) = 1; toc
% Elapsed time is 0.000006 seconds.

B = A;
tic; A(1) = 0; toc
% Elapsed time is 0.511793 seconds.