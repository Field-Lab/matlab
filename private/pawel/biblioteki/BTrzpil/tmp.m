function [ p, d ] = tmp()
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
a=[1,2,3];
b=[1,2];
p=a;
d=horzcat(a,b);
a=[-3,2,1];
b=[4,5,-6];
c = 0;
lista = {a, b, c};
d = cell2mat(lista(1))
end


