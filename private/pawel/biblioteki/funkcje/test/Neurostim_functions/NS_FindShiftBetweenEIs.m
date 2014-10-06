function [Shift,MaxCorr]=NS_FindShiftBetweenEIs(EI1,EI2);
%The function finds the shift between EIs in time. The procedure is based
%on cross correlation, which is calculated for each channel, and then
%summaraized for all the channels. Positive shift means that EI2 is delayed
%in reference to EI1.
SEI1=size(EI1);
SEI2=size(EI2);

if SEI1(1)~=SEI2(1)
    error('Different number of channels for EIs');
end
c=zeros(1,2*max(SEI1(2),SEI2(2))-1);

for i=1:SEI1(1);
    a=xcorr(EI1(i,:),EI2(i,:));
    c=c+a;
end
MaxCorr=max(c);
b=find(c==MaxCorr);
b=b(1); % in case the length(b) is more than 1 - unlikely

Shift=b-SEI2(2);