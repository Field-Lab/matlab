function ArtF=FilterArtifactLocalStim(Kers,Art,x,cmax,Q,Qt,krondiaginv)
%Filter the Artifact (mean of spike substracted traces accross repetitions
%of the same stimulus amplitudes) using the model and noise level. Only for
%Stimulating electrodes
%Here, x(end) contains log variance sigma0.
%Gonzalo Mena, 03/2016


varsmax=length(x);
logs=0;
vec=[];

Kers{2}=Kers{2}(cmax,cmax);
Art=Art(cmax,:);
Q{2}=1;
Qt{2}=1;

prod1=KronProd(Q,Art);
prod1=krondiaginv.*prod1;

alpha=KronProd(Qt,prod1);

prod2=KronProd(Q,Art).*krondiaginv.^(1/2);


Filtered=exp(x(varsmax-1))*KronProd(Kers,alpha);
ArtF=reshape(Filtered,size(Art,1),size(Art,2));
