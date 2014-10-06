function [Artifact]=NS512_ArtifactEstimationSomatic(Traces,N);
%
%This function estimates the artifact shape. For the function to work 
%properly, the "Traces" should include at
%least "N" events with no spikes (artifacts only), however, if the number
%of "artifacts only" is just slightly lower, the estimated artifact shape
%will be only slightly distorted.
%The function works as following:
%1) Average all the Traces and subtract the average from all Traces
%2) Find minimal value for each trace
%3) Sort all the traces from the point of view of the minimal value (most
%negative to most positive) - hopefully, the "artifact only" traces are now
%indexed by 1 to something. The sorting is done independently for each
%channel!
%4) Take N first traces, average and this is the estimated artifact
%N - how many traces are taken for estimation of the artifact shape

Traces0=Traces;
ST=size(Traces);
b=mean(Traces);

for i=1:ST(1)
    Traces(i,:,:)=Traces(i,:,:)-b;
end
c1=min(Traces,[],3);
c2=max(Traces,[],3);

% 1. Find 10 artifacts bsed on most negative samples
[d,e]=sort(c1+c2,1); %sort traces for each channel independently, from smallest to largest maximum value to idnetify "artifact only" traces

Artifact=zeros(1,ST(2),ST(3));




Artifact=zeros(1,ST(2),ST(3));

for i=1:ST(2)
    
    NumberOfArtifacts=N;
    
    for n=N:3
        if  max(abs(mean(Traces0(e(ST(1)-N+1:ST(1),i),i,:))))<3*std(mean(Traces0(e(ST(1)-N+1:ST(1),i),i,:))) %najwieksze odchylenie nie wiêksze ni¿ 3*odchylenie standardowe
            break
        else 
            NumberOfArtifacts=NumberOfArtifacts-1;  %je¿eli wiêksze to zmniejszamy ilosc artefaktow branych do sredniej
        end
     
    end
    
    Artifact(1,i,:)=mean(Traces0(e(ST(1)-NumberOfArtifacts+1:ST(1),i),i,:));
end

