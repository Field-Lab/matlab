% Problem 4b
clear all; close all; clc

r = [1.05 0.06 1.06 0.07 0.07 1.07]';
alpha1 = [0.35 0.35 0 0.3 0 0]';
alpha2 = [0.6 0.2 0 0.2 0 0]';
M = [0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;];
A1 = M + alpha1*r';
A2 = M + alpha2*r';

lambda1 = eig(A1);
G1 = lambda1(lambda1==max(lambda1));

lambda2 = eig(A2);
G2 = lambda2(lambda2==max(lambda2));

% Problem 4c
beta1 = [1 0 1 0 0 1]';
beta2 = [0 1 0 0 1 0]';
beta3 = [0 0 0 1 0 0]';
[v1_1,~] = eigs(A1,1);
[v1_2,~] = eigs(A2,1);

L1_1 = beta1'*v1_1/sum(v1_1);
L2_1 = beta2'*v1_1/sum(v1_1);
L3_1 = beta3'*v1_1/sum(v1_1);

L1_2 = beta1'*v1_2/sum(v1_2);
L2_2 = beta2'*v1_2/sum(v1_2);
L3_2 = beta3'*v1_2/sum(v1_2);

% Problem 4d
% get the left eigenvectors using A transpose
[w1,~] = eigs(A1',1);