function [A,B]=NS_GeneratePairs(N);
% N must be even!!
NumberOfPairs=floor(N/2);
if NumberOfPairs*2~=N
    error('N is not even');
end
A=zeros(N-1,NumberOfPairs,2); %N runds, 

PairsPre=zeros(1,N);
for i=1:N-1 %for each rund...
    x=[1 N+2-i:N 2:N-i+1];
    A(i,:,1)=x(1:NumberOfPairs);
    A(i,:,2)=x(N:-1:NumberOfPairs+1);
end

B=zeros((N-1)*NumberOfPairs,2);
for i=1:N-1
    B(NumberOfPairs*(i-1)+1:NumberOfPairs*(i-1)+NumberOfPairs,:)=A(i,:,:);
end