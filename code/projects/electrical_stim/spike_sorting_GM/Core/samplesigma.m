function Gibbs = samplesigma(Gibbs)

b0 = Gibbs.params.b0; 
a0 = Gibbs.params.a0;

J  = Gibbs.params.J;
E  = Gibbs.params.E;
I  = Gibbs.params.I;
T  = Gibbs.params.T;
Tdivide   = Gibbs.params.Tdivide;
Residuals = Gibbs.variables.Residuals;

for e = 1:E
    for j = 1:J

        a(e,j) = a0(j)+I(j)*T/2;
        b(e,j) = b0(j)+1/2*sum(sum(Residuals{j,e}.^2));
    end
  sigma(e,:)   = sqrt(1./gamrnd(a(e,:),1./b(e,:)));  
end
 
Gibbs.variables.sigma = sigma;

for t=1:length(Tdivide)-1
    for e = 1:E
        for j=1:J
        
            Interval = Tdivide(t)+1:Tdivide(t+1);
            at(e,j) = a0(j)+I(j)*length(Interval)/2;
            bt(e,j) = b0(j)+1/2*sum(sum(Residuals{j,e}(:,Interval).^2));    
        end
        sigmaT{t}(e,:) = sqrt(1./gamrnd(a(e,:),1./b(e,:))); 
    end
end
Gibbs.variables.sigmaT = sigmaT;
