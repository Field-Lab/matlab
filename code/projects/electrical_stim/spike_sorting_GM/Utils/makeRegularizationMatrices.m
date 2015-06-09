function [matricesReg]=makeRegularizationMatrices(breakRanges,Tdivide)

E=length(breakRanges);
J=breakRanges{1}(end);



T=Tdivide(end);

for m=1:length(Tdivide)-1
        
    for e = 1:E
   
        matricesReg(e).RegT{m}=sparse(0,0);
        
        Tr=[Tdivide(m)+1:Tdivide(m+1)];

            for j = 1:J
  
            indj=zeros(1,J);
            indj(j)=1;
    
                for t=1:length(Tr)-1
                
                    indt = zeros(1,T);
                    indt(Tr(t)) = -1;
                    indt(Tr(t)+1) = 1;
                    indtj = kron(indj,indt);
                    
                    matricesReg(e).RegT{m} = sparse([matricesReg(e).RegT{m};sparse(indtj)]);
                end
            end
    end
end

for e=1:E

    conta=1;

    for b=1:length(breakRanges{e})-1
    
        interval=breakRanges{e}(b)+1:breakRanges{e}(b+1);
    
        if(length(interval)<=1)
        
          continue
        else
        
        matricesReg(e).RegJ{conta}=[];
    
        for j = 1:length(interval)-1
            indj=sparse(1,J);
            indj(interval(j+1))=1;
            indj(interval(j))=-1;       
                
            
            for t=1:T
                indt=sparse(1,T);
                indt(t)=1;
                indjt=kron(indj,indt);
                matricesReg(e).RegJ{conta}=sparse([matricesReg(e).RegJ{conta};sparse(indjt)]);
            end
        end
        conta=conta+1;
        end
    
    end
end

sizeMatrices = size(matricesReg(1).RegT{1},2);

for e=1:E
    matricesReg(e).RegOverall{1} = speye(sizeMatrices);
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

