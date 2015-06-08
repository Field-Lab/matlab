function Gibbs  =   sampleSpikes(Gibbs)
%YJ{j}=Data with colums=timesx electrodes, i=trials, j=condition
%Xj=Design matrix for condition J
%beta=current estimate of beta
%alpha=current estimate of logistic regression

dataVecJ    = Gibbs.params.dataVecJ;
TfindRel0   = Gibbs.params.TfindRel0;
Tfind0      = Gibbs.params.Tfind0;
lambdas     = Gibbs.params.lambdas;
Xj          = Gibbs.params.Xj;
X           = Gibbs.params.X;    
Kn          = Gibbs.params.Kn;   
I           = Gibbs.params.I;
E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
nNeurons    = Gibbs.params.nNeurons;
Imax        = max(I);

sigma            = Gibbs.variables.sigma;
Artifact         = Gibbs.variables.Artifact; 
Probs            = Gibbs.variables.Probs;     
ActionPotentials = Gibbs.variables.ActionPotentials;




for j=1:J
    for e=1:E
        sigmaVecs(j,1+(e-1)*T:e*T) = sigma(e,j);
    end
end


for j=1:J
    
    ArtifactRep     =  repmat(Artifact(j,:),I(j),1);
    ResArtifactJ{j} =  dataVecJ{j}-ArtifactRep;
    
end

Probs   = min(Probs,1);
Probs   = max(Probs,0);

   

for n=1:nNeurons
       
      inactiveNeurons = setdiff([1:nNeurons],n);

      for j=1:J
        
        ResidualsJ{j} = ResArtifactJ{j};
      
        prob(j)      = Probs(n,j);
        probTimes    = repmat(prob(j)/(length(TfindRel0)-1),1,length(TfindRel0),1);
        probTimes(1) = 1-prob(j);

        for ni=1:length(inactiveNeurons)

            ResidualsJ{j}=ResidualsJ{j}-ActionPotentials{inactiveNeurons(ni),j};
        
        end
        
        for i=1:I(j)
            
            for t=1:length(TfindRel0)
                normexp(t)=-norm((ResidualsJ{j}(i,:)-Kn{n}(:,t)')./(sigmaVecs(j,:)),2).^2./(2);
                logLiks{j,n}(i,t)=normexp(t)-T*E/2*log((2*pi))-sum(log(sigmaVecs(j,:)));
           
            end
            
            logPost                    = log(probTimes)+logLiks{j,n}(i,:);
            logPosts{j,n}(i,:)         = logPost;
            logPostShifted             = logPost-max(logPost); 
            postUnnormalized           = exp(logPostShifted);
            post                       = postUnnormalized/sum(postUnnormalized);
            r                          = mnrnd(1,post);
            indLatency                 = find(r==1);
            spikes{n}(j,i)             = indLatency>1;
            latencies{n}(j,i)          = Tfind0(indLatency);
            ActionPotentials{n,j}(i,:) = (Kn{n}(:,indLatency))';
            
        end
        
        if(I(j) < Imax)
            spikes{n}=double(spikes{n});
            latencies{n}=double(latencies{n});
            spikes{n}(j,I(j)+1:Imax)    = NaN;
            latencies{n}(j,I(j)+1:Imax) = NaN;
        end
        
      end
end


[ResSpikesE ResSpikesEJ]=substractActionPotentials(dataVecJ,ActionPotentials,E);
  
  


  
  Gibbs.variables.spikes           = spikes;
  Gibbs.variables.latencies        = latencies;
  Gibbs.variables.ActionPotentials = ActionPotentials;
  Gibbs.variables.ResSpikesEJ      = ResSpikesEJ;
  Gibbs.variables.ResSpikesE       = ResSpikesE;
  
  