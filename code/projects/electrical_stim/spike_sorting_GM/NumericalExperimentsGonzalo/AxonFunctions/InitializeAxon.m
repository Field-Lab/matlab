function initial=InitializeAxon(input)



data           = input.tracesInfo.data;
templates      = input.neuronInfo.templates;
TfindRange     = input.params.initial.TfindRange;  
TfindRangeRel  = input.params.initial.TfindRangeRel;
TfindRel       = [TfindRangeRel(1):TfindRangeRel(2)];
Tdivide        = input.params.initial.Tdivide;
rho            = input.params.initial.rho;

E = input.tracesInfo.E;
T = input.tracesInfo.T;
I = input.tracesInfo.I;
J = input.tracesInfo.J;
nNeurons = input.neuronInfo.nNeurons;


      
lengthFindRange = length(TfindRel);
lengthSpikes    = sum(I)*nNeurons*lengthFindRange;
lengthDataVec   = sum(I)*E*T;


dataVec=[];

for j=1:J
    for i=1:I(j)
        for e=1:E
            dataVec=[dataVec;data{j,e}(i,:)'];
        end
    end
end







%HardWare Artifact Model


for e=1:E
    breakRanges{e}  = [0 input.tracesInfo.breakRecElecs{e} J];
    
    for m=1:length(breakRanges{e})-1
        lengthRange    = breakRanges{e}(m+1)-breakRanges{e}(m);
        aux            = find(lengthRange<=input.params.initial.degPolRule)-1;
        degPols{e}(m)  = aux(1);
    end
end


for e = 1:E
    sizeAe(e)       = T*sum(degPols{e}+1);
    sizeAeDegCum{e} = [0 T*cumsum(degPols{e}+1)];
end
sizeAeCum = [0 cumsum(sizeAe)];
sizeA     = sum(sizeAe);


Xpol=sparse(0,0);
for j = 1:J

    XjPol{j}=sparse(T*E,sizeA);
    
    for e = 1:E

    XjAlle=sparse(T,sizeAe(e));
    rangeIndex=find(breakRanges{e}<j);
    rangeIndex=rangeIndex(end);
    covPol=(j-breakRanges{e}(rangeIndex));
    covsPol=[];

        for p = 1:degPols{e}(rangeIndex)+1
            covsPol=[covsPol covPol^(p-1)];
        end

    Xje{j,e}=sparse(0,0);
    
        for t=1:T
            indt         = sparse(1,T);
            indt(t)      = 1;
            indtcovs     = sparse(kron(indt,covsPol));
            Xje{j,e}  = sparse([Xje{j,e};indtcovs]);
        end
        
        XjAlle(:,1+sizeAeDegCum{e}(rangeIndex):sizeAeDegCum{e}(rangeIndex+1)) = Xje{j,e};
        XjPol{j}(T*(e-1)+1:T*e,1+sizeAeCum(e):sizeAeCum(e+1)) = XjAlle;
       
    end
    
        Xjipol = sparse(repmat(XjPol{j},I(j),1));
        Xpol  = sparse([Xpol;Xjipol]);

end
    
for e=1:E
    for j=1:J
        XePol{e}((j-1)*T+1:j*T,:) = XjPol{j}((e-1)*T+1:e*T,1+sizeAeCum(e):sizeAeCum(e+1));
    end
end




%Axonal Bundle Model



for e=1:E
    breakRangesAxon{e}  = [input.tracesInfo.breakAxon{e} J];
    
    for m=1:length(breakRangesAxon{e})-1
        lengthRangeAxon        = breakRangesAxon{e}(m+1)-breakRangesAxon{e}(m);
        aux                = find(lengthRangeAxon<=input.params.initial.degPolRule)-1;
        degPolsAxon{e}(m)  = aux(1);
    end
end

degAxon=2;
degPolsAxon{e}=degAxon;
rangeIndex = 1;
for e = 1:E
    sizeAeAxon(e)       = T*sum(degPolsAxon{e}+1);
    sizeAeAxonDegCum{e} = [0 T*cumsum(degPolsAxon{e}+1)];
end
sizeAeAxonCum = [0 cumsum(sizeAeAxon)];
sizeAAxon     = sum(sizeAeAxon);


XpolAxon=sparse(0,0);
for j = 1:J

    XjPolAxon{j}=sparse(T*E,sizeAAxon);
    
    for e = 1:E

    XjAlleAxon=sparse(T,sizeAeAxon(e));

    covPol=max(0,(j-breakRangesAxon{e}(1)-1));
    covsPol=[];
%change

        for p = 1:degAxon+1 
            if(j<=breakRangesAxon{e}(1))
                covsPol =[covsPol 0];
            else
            covsPol=[covsPol covPol^(p-1)];
            end
        end

    XjeAxon{j,e}=sparse(0,0);
    
        for t=1:T
            indt         = sparse(1,T);
            indt(t)      = 1;
            indtcovs     = sparse(kron(indt,covsPol));
            XjeAxon{j,e}  = sparse([XjeAxon{j,e};indtcovs]);
        end
        
        XjAlleAxon(:,1+sizeAeAxonDegCum{e}(rangeIndex):sizeAeAxonDegCum{e}(rangeIndex+1)) = XjeAxon{j,e};
        XjPolAxon{j}(T*(e-1)+1:T*e,1+sizeAeAxonCum(e):sizeAeAxonCum(e+1)) = XjAlleAxon;
       
    end
    
        XjipolAxon = sparse(repmat(XjPolAxon{j},I(j),1));
        XpolAxon  = sparse([XpolAxon;XjipolAxon]);

end
    
for e=1:E
    for j=1:J
        XePolAxon{e}((j-1)*T+1:j*T,:) = XjPolAxon{j}((e-1)*T+1:e*T,1+sizeAeCum(e):sizeAeCum(e+1));
    end
end






K  = sparse(0,0);
K0 = makeToeplitz(templates,TfindRel,T);

cumsumSpikes = [0 cumsum(I*nNeurons*lengthFindRange)];

for j = 1:J
    Kj    = sparse(E*T*I(j),lengthSpikes);
    Kjaux = sparse(kron(speye(I(j)),K0));
    Kj(:,cumsumSpikes(j)+1:cumsumSpikes(j+1)) = Kjaux;
    K = sparse([K;Kj]);
end



indt = ones(1,lengthFindRange);
%at most one spike per neuron per trial
Asp = sparse(0,0);

for j=1:J
    
    for i=1:I(j)
    
        for n=1:nNeurons
            
            indn    = zeros(1,nNeurons);
            indn(n) = 1;
            indtn   = kron(indn,indt);
            indi    = zeros(1,I(j));
            indi(i) = 1 ;
            indtni  = sparse(kron(indi,indtn));
            ind     = zeros(1,lengthSpikes);
            ind(1,1+cumsumSpikes(j):cumsumSpikes(j+1)) = indtni;
            Asp     = sparse([Asp;ind]);
        end
      
    end
end


cumsumI = [0 cumsum(I*nNeurons*lengthFindRange)];

Ai = sparse(0,0);
%increasing spike probabilities
for n = 1:nNeurons

    for j = 1:J-1
         indi1      = ones(1,I(j))/I(j);
         indi2      = -ones(1,I(j+1))/I(j+1);
         indn       = zeros(1,nNeurons);
         indn(n)    = 1;
         indtn      = kron(indn,indt);
         indtni1    = sparse(kron(indi1,indtn));
         indtni2    = sparse(kron(indi2,indtn));
         ind        = zeros(1,lengthSpikes);
         ind(1,1+cumsumI(j):cumsumI(j+1)) = indtni1;
         ind(1,1+cumsumI(j+1):cumsumI(j+2))=indtni2;
         Ai = sparse([Ai;ind]);
    end
    
end

%
Ag = sparse([speye(lengthSpikes);-speye(lengthSpikes);Asp;Ai]);
b  = sparse([ones(lengthSpikes,1);zeros(lengthSpikes,1);ones(size(Asp,1),1);zeros(size(Ai,1),1)]);




c = ones(lengthSpikes,1);
cvx_begin
variable BetaPol(sizeA)
variable BetaPolAxon(sizeAAxon)
variable s(lengthSpikes) 
if(rho>0)
    minimize (rho*c'*s+quad_form((dataVec-Xpol*BetaPol-XpolAxon*BetaPolAxon-K*s),speye(length(dataVec))));
else
    minimize (quad_form((dataVec-Xpol*BetaPol-XpolAxon*BetaPolAxon-K*s),speye(length(dataVec))));
end
Ag*s <= b;
cvx_end

for e=1:E
    BetaPolE{e}=BetaPol(1+sizeAeCum(e):sizeAeCum(e+1));
end

for e=1:E
    BetaPolAxonE{e}=BetaPolAxon(1+sizeAeAxonCum(e):sizeAeAxonCum(e+1));
end



for n = 1:nNeurons
    for j = 1:J
        for i = 1:I(j)
            ind                         = cumsumI(j)+lengthFindRange*nNeurons*(i-1)+lengthFindRange*(n-1);
            GeneralizedSpikes{n,j}(:,i) = s(ind+1:ind+lengthFindRange);
        end
    end
end

for j = 1:J
    for n = 1:nNeurons
            Probs(n,j) = nanmean(nansum(GeneralizedSpikes{n,j})); 
    end
end



for j = 1:J
    A(j,:) = XjPol{j}*BetaPol;
    
end

for j = 1:J
    Axon(j,:) = XjPolAxon{j}*BetaPolAxon;
end


for j = 1:J
    AxonE{E} = Axon(:,1+(e-1)*T:e*T);
    
end

for e = 1:E
    AE{e} = A(:,1+(e-1)*T:e*T);
end


for j=1:J
    for n=1:nNeurons
        indn    = zeros(1,nNeurons);
        indn(n) = 1;
        ActionPotentials{n,j} = (K0*kron(indn',GeneralizedSpikes{n,j}))';
    end
end    

for j = 1:J
    sumActionPotentials = 0;
    
    for n = 1:nNeurons
        sumActionPotentials = sumActionPotentials+ActionPotentials{n,j};
        
        for e=1:E
           
            Residuals{j,e}  = data{j,e} - repmat(AE{e}(j,:),I(j),1) - sumActionPotentials(:,T*(e-1)+1:T*e);
            sigma(e,j)      = sqrt(nansum(nansum((Residuals{j,e}.^2)))/(I(j)*T));
        
        end
    end
end

 

[X Xj]        = makeArtifactCovariates(T,J,I);
[matricesReg] = makeRegularizationMatrices(breakRanges,Tdivide);


[matricesRegAxon]=makeRegularizationMatricesAxon(input);
for e=1:E
    lReg(e) =length(matricesReg(e).Prods);
    lRegAxon(e) =length(matricesRegAxon(e).Prods);
    for l=1:lReg(e)
        matricesRegAll(e).Prods{l}=matricesReg(e).Prods{l};
    end
    for l=1:lRegAxon(e)
        matricesRegAll(e).Prods{l+lReg(e)}=matricesRegAxon(e).Prods{l};
    end
    
end

for e=1:E
    clear quad  
    Beta(:,e)  = reshape(AE{e}',T*J,1);
    for l = 1:length(matricesReg(e).Prods)
       quad(l) = trace(Beta(:,e)'*matricesReg(e).Prods{l}*Beta(:,e));
    end
    BetaAxon(:,e)  = reshape(AxonE{e}',T*J,1);
    
    for l = 1:length(matricesRegAxon(e).Prods)
       quadAxon(l) = trace(BetaAxon(:,e)'*matricesRegAxon(e).Prods{l}*BetaAxon(:,e));
    end
    quads =[quad quadAxon];
    
    lambda0    = ones(lReg(e)+lRegAxon(e),1);
    lambda     = NewtonMaxLogDet(input,lambda0,matricesRegAll(e).Prods,quads);
    lambda     = exp(lambda);
    lambdas1{e} = lambda(1:lReg(e));
    lambdas2{e} = lambda(lReg(e)+1:end);
    lambdas{e}  =  lambda;


Lambdas{e} = 0;


for r = 1:length(matricesReg(e).Prods)
    
    Lambdas{e} = Lambdas{e}+lambdas{e}(r)*matricesReg(e).Prods{r};

end

Lambdas{e} = sparse(Lambdas{e});
LambdasInv{e} = sparse(inv(Lambdas{e}));
end

initial.BetaPolE          = BetaPolE;
initial.BetaPol           = BetaPol;
initial.BetaPolAxonE      = BetaPolAxonE;
initial.BetaPolAxon       = BetaPolAxon;
initial.ArtifactVariables = Beta;
initial.Artifact          = A;
initial.AxonArtifact      = Axon;
initial.AxonArtifactE     = AxonE;
initial.ArtifactE         = AE;
initial.Residual          = Residuals;
initial.sigma             = sigma;
initial.Probs             = Probs;
initial.GeneralizedSpikes = GeneralizedSpikes;
initial.ActionPotentials  = ActionPotentials;

initial.params.Xj                 = Xj;
initial.params.X                  = X;
initial.params.Xpol               = Xpol;
initial.params.XePol              = XePol;
initial.params.XjPol              = XjPol;
initial.params.matricesReg        = matricesReg;
initial.params.lambdas            = lambdas;
initial.params.lambdas1           = lambdas1;
initial.params.lambdas2           = lambdas2;
initial.params.matricesRegAxon    = matricesRegAxon;
initial.params.Lambdas            = Lambdas;
initial.params.LambdasInv         = LambdasInv;
initial.params.a0                 = -0.5*ones(1,J);
initial.params.b0                 =  zeros(1,J);
initial.params.lambdaLogReg       = 0.0001;
initial.params.alphaLogReg        = 0.0001;