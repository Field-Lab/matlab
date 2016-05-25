function [Ker KerD]=evalKernels(Dif,x,vars,type,varargin)
%Dif = Diferences between covariates (matrix)
%x = non-stationarity vector (diagonal of Dif usually) (also, for stimulating electrode, indicate
%breakpoints)
%vars = values of the variables, depends on the type
%type of Kernel (see below)
% Gonzalo Mena, 03/2016

if(type==1)
    %Most usual, for 3 Kernels with same structure
    lambda=exp(vars(1));
    beta=exp(vars(2));
    alpha=exp(vars(3));
    
    K=1;
    
    
    
    Dg=K*exp(-x*beta).*x.^alpha;
    DgK=Dg;
    Dga=DgK.*alpha.*log(x);
    Dga(isinf(Dga))=0;
    Dga(isnan(Dga))=0;
    Dgb=-DgK*beta.*x;
    Ker0=(1+sqrt(3)*(Dif)*lambda).*exp(-sqrt(3)*(Dif)*lambda);
    Ker=(Dg*Dg').*Ker0;
    
    Ker0l=-3*lambda*exp(-sqrt(3)*(Dif)*lambda).*lambda.*(Dif).^2;
    KerD{1}=(Dg*Dg').*Ker0l;
    KerD{2}=(Dg*Dgb'+Dgb*Dg').*Ker0;
    KerD{3}=(Dg*Dga'+Dga*Dg').*Ker0;
    
elseif(type==2)
    %for Stimulating electrode, include information of brakpoints in x
    lambda=exp(vars);
    
    lambdas=zeros(size(Dif,1));
    ind=zeros(size(Dif,1));
    for k=1:length(lambda)
        ma=zeros(size(Dif,1));
        %ma(x(k)+1:x(k+1),x(k)+1:x(k+1))=lambda(k);
        ma(1:x(2)-x(1),1:x(2)-x(1))=lambda(k);
        lambdas=lambdas+ma;
        %ind(x(k)+1:x(k+1),x(k)+1:x(k+1))=1;
        ind(1:x(2)-x(1),1:x(2)-x(1))=1;
        KerD{k}=-3*ma.*exp(-sqrt(3)*(Dif).*ma).*ma.*(Dif).^2;
        
    end
    Ker=(ind+sqrt(3)*(Dif).*lambdas).*exp(-sqrt(3)*(Dif).*lambdas);
    
    
elseif(type==3)
    %same as 1, but allows for additional variance
    lambda=exp(vars(1));
    beta=exp(vars(2));
    alpha=exp(vars(3));
    sigma=exp(vars(4));
    K=1;
    
    Dg=K*exp(-x*beta).*x.^alpha;
    DgK=Dg;
    Dga=DgK.*alpha.*log(x);
    Dga(isinf(Dga))=0;
    Dga(isnan(Dga))=0;
    Dgb=-DgK*beta.*x;
    Ker0=(1+sqrt(3)*(Dif)*lambda).*exp(-sqrt(3)*(Dif)*lambda);
    Ker=(Dg*Dg').*Ker0+sigma*eye(size(Ker0,1));
    
    Ker0l=-3*lambda*exp(-sqrt(3)*(Dif)*lambda).*lambda.*(Dif).^2;
    KerD{1}=(Dg*Dg').*Ker0l;
    KerD{2}=(Dg*Dgb'+Dgb*Dg').*Ker0;
    KerD{3}=(Dg*Dga'+Dga*Dg').*Ker0;
    KerD{4}=sigma*eye(size(Ker0,1));
    
    
elseif(type==4)
    %for post-bundle Kernels
    lambda=exp(vars(1));
    beta=exp(vars(2));
    sigma=exp(vars(3));
    K=1;
    
    Dg=x.^(1/2);
    
    Ker0=(1+sqrt(3)*(Dif)*lambda).*exp(-sqrt(3)*(Dif)*lambda);
    Ker=(Dg*Dg').*Ker0;
    KerD{1}=Ker0;
    KerD{2}=Ker0;
    KerD{3}=Ker0;
    
 
elseif(type==5)
    KerPrev=varargin{1};
    %for post-bundle Kernels
    lambda=exp(vars(2));
    K=exp(vars(1));
    sigma=exp(vars(3));
    
    Dg=x.^(1/2);
    
    Ker0=(1+sqrt(3)*(Dif)*lambda).*exp(-sqrt(3)*(Dif)*lambda);
    
    Ker0l=-3*lambda*exp(-sqrt(3)*(Dif)*lambda).*lambda.*(Dif).^2;
   
    Ker=KerPrev+K*(Dg*Dg').*Ker0;
    KerD{1}=K*(Dg*Dg').*Ker0;
    KerD{2}=K*(Dg*Dg').*Ker0l;
    KerD{3}=Ker0;
    
    
elseif(type==6)
    KerPrev=varargin{1};
    %for post-bundle Kernels
    lambda=exp(vars(2));
    K=exp(vars(1));
    sigma=exp(vars(3));
    
    Dg=x;
    
    %%Ker0=(1+sqrt(3)*(Dif)*lambda).*exp(-sqrt(3)*(Dif)*lambda);
    
    %Ker0l=-3*lambda*exp(-sqrt(3)*(Dif)*lambda).*lambda.*(Dif).^2;
   
    Ker=KerPrev+K*diag(Dg);
    KerD{1}=K*diag(Dg);
    %KerD{2}=K*(Dg*Dg').*Ker0l;
    %KerD{3}=Ker0;
    
end


