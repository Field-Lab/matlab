function [val,grad]=logDetKronStimElec(Art,var,Difs,varsactive,types,Diags,nvar)
%Gonzalo Mena, 3/2016
nvarcum=cumsum([0 nvar]);

nTotal=prod(size(Art));

logs=0;

varsmax=length(var);
if(max(exp(var))==Inf)
    val = Inf;
    grad = Inf*ones(length(var),1);
    return
end

logs=0;
vec=[];
krondiag=1;


for k=1:2
    

[Ker KerD]=evalKernels(Difs{k},Diags{k},var(nvarcum(k)+1:nvarcum(k+1)),types(k));
    
[a b]=eig(Ker);
Q{k}=a';
Qt{k}=a;
dL{k}=diag(b);
krondiag=kron(krondiag,dL{k});
KerDs{k}=KerD;
Kers{k}=Ker;
KersQ{k}=diag(Q{k}*Kers{k}*Qt{k});

end

krondiag=exp(var(varsmax-1))*krondiag+exp(var(varsmax));
logs=logs+sum(log(krondiag));

prod1=KronProd(Q,Art);
krondiaginv=krondiag.^(-1);
prod1=krondiaginv.*prod1;

alpha=KronProd(Qt,prod1);


alpha2=dot(Art(:),alpha);
val=alpha2+logs;

for j=1:varsmax-1
    g{j}=Kers;
    diagGQ{j}=KersQ;
end


for k=1:2
    
    for i=1:nvar(k)
        if(~isempty(intersect(varsactive,nvarcum(k)+i)))
        g{nvarcum(k)+i}{k}=exp(var(varsmax-1))*KerDs{k}{i};
        diagGQ{nvarcum(k)+i}{k}=exp(var(varsmax-1))*diag(Q{k}*KerDs{k}{i}*Qt{k});
        end
    end
    
end

g{varsmax-1}{1}=g{varsmax-1}{1}*exp(var(varsmax-1));
diagGQ{varsmax-1}{1}=exp(var(varsmax-1))*diagGQ{varsmax-1}{1};

g{varsmax}{1}=exp(var(varsmax))*speye(size(Kers{1},1)*size(Kers{2},1));

diagGQ{varsmax}{1}=diag(exp(var(varsmax))*speye(size(Kers{1},1)));
diagGQ{varsmax}{2}=diag(exp(var(varsmax))*speye(size(Kers{2},1)));



for j=1:length(varsactive)
    varac=varsactive(j);
    d=1;
    
    for k=1:length(diagGQ{varac})
        d=kron(d,diagGQ{varac}{k});
    end

    grad(j)=dot(krondiaginv,d)-alpha'*KronProd(g{varac},alpha);
    
end


grad=grad';



end
