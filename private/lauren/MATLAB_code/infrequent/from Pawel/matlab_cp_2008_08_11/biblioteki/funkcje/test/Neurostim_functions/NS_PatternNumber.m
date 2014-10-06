function [PatternNumber,number]=NS_PatternNumber(Pattern,Indexes);
%number - number of patterns in Indexes that fit the Pattern variable.
%Should not be greater than 1, basucally.

number=0;
SI=size(Indexes);
PatternNumber=0;
for i=1:SI(1)
    %Size(
    Indexes(i,:);
    Pattern;
    if Indexes(i,:)==Pattern
        PatternNumber=i;
        number=number+1;
    end
end        