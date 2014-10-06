function [ output_args ] = ShiftStartIndexes_PR( TimeCorelSpikes, UniSpikesIndic )
%UNTITLED Summary of this function goes here
%Funkcja nie skonczona!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   Detailed explanation goes here
g1=find(TimeCorelSpikes==1);
if length(TimeCorelSpikes)~=length(UniSpikesIndic)
    error('The size of TimeCorelSpikes array does not match length of UniSpikesIndic array');
end else if size(g1)
            g2=min(UniSpikesIndic(g1));
            if g2>4
                UniSpikesIndic=UniSpikesIndic-4;
            else if g2>3 & g2<5
                %UniSpikesIndic=UniSpikesIndic-g2+1;
                UniSpikesIndic=UniSpikesIndic-3;
                else if g2>2 & g2<4
                        UniSpikesIndic=UniSpikesIndic-2;
                    else if g2>1 & g2<3
                            UniSpikesIndic=UniSpikesIndic-1;
                        end
                    end
                end
            end
        end
end

