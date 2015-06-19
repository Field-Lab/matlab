function [matricesRegAxon]=makeRegularizationMatricesAxon(input)

J=input.tracesInfo.J;
E=input.tracesInfo.E;
T=input.tracesInfo.T;
Tdivide = input.params.initial.Tdivide;
breakAxon =input.tracesInfo.breakAxon;
for e=1:E
    lengthRangeAxon(e) = J-breakAxon{e};
end

Tdivide =[0 T];


for m=1:length(Tdivide)-1
        
    for e = 1:E
   
        matricesReg(e).RegT{m}=sparse(0,0);
        
        Tr=[Tdivide(m)+1:Tdivide(m+1)];

            for j = 1:lengthRangeAxon
  
            indj=zeros(1,lengthRangeAxon);
            indj(j)=1;
    
                for t=1:length(Tr)-1
                
                    indt = zeros(1,T);
                    indt(Tr(t)) = -1;
                    indt(Tr(t)+1) = 1;
                    indtj = kron(indj,indt);
                    
                    matricesReg(e).RegT{m} = sparse([matricesReg(e).RegT{m};sparse(indtj)]);
                end
            end
            matricesReg(e).RegT{m} = sparse(blkdiag(sparse(T*breakAxon{e},T*breakAxon{e}),matricesReg(e).RegT{m}));
    end
end

for e=1:E
    
    
    if(lengthRangeAxon(e)<=1)
        
        continue
    else
        
        matricesReg(e).RegJ{1}=[];
        
        for j = 1:lengthRangeAxon(e)-1
            indj=sparse(1,lengthRangeAxon(e));
            indj(j+1)=1;
            indj(j)=-1;
            
            
            for t=1:T
                indt=sparse(1,T);
                indt(t)=1;
                indjt=kron(indj,indt);
                matricesReg(e).RegJ{1}=sparse([matricesReg(e).RegJ{1};sparse(indjt)]);
                
            end
        end
        
    end
    matricesReg(e).RegJ{1}=sparse(blkdiag(sparse(T*breakAxon{e},T*breakAxon{e}),sparse([matricesReg(e).RegJ{1}])));
end


sizeMatrices = T*lengthRangeAxon(e);

for e=1:E
    matricesReg(e).RegOverall{1} = speye(sizeMatrices);
    matricesReg(e).RegOverall{1} =sparse(blkdiag(sparse(T*breakAxon{e},T*breakAxon{e}),matricesReg(e).RegOverall{1}));
end


for e=1:E
    
    matricesReg(e).Prods{1} = matricesReg(e).RegOverall{1};
    
    for m=1:length(matricesReg(e).RegJ)
        
        matricesReg(e).Prods{m+1} = matricesReg(e).RegJ{m}'*matricesReg(e).RegJ{m};
    end
    
    m0 = length(matricesReg(e).Prods);
    
    for t=1:length(matricesReg(e).RegT)
        
        matricesReg(e).Prods{t+m0} = matricesReg(e).RegT{t}'*matricesReg(e).RegT{t};
        
    end
end

matricesRegAxon = matricesReg;