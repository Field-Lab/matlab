function lambda = NewtonMaxLogDet(lambda0,Prods,quad)

beta  = 0.5;
a     = 0.5;
dif   = Inf;
thres = 0.001;

val = logDetExp(log(lambda0),Prods,quad);

if(val == Inf)
    
    lambda=10^(-8)*ones(size(Prods,2),1);
    
    while(true)
        
        val = logDetExp(log(lambda),Prods,quad);
        if(val == Inf)
             lambda = lambda*2;
        else
             lambda = log(lambda);
             break
        end
    end
else
    
    lambda = log(lambda0);
    
end


val0 = val;



[val grad]=logDetExp(lambda,Prods,quad);

cont = 1;
t    = 1;

while(true)
  
    lambdanew = lambda-t*grad;
    val2      = val-a*t*norm(grad)^2; 

    if(logDetExp(lambdanew,Prods,quad) > val2)
        t    = beta*t;
        cont = cont+1;
    else
        lambda     = lambdanew;
        [val grad] = logDetExp(lambda,Prods,quad);
        dif        = val0-val;
        val0       = val;
        t          = 1;
        
        disp(['Iteration ' num2str(cont) ' finished'])
        disp(['Diference in function evaluation with previous iteration = ' num2str(dif)])
        disp(['Gradient = '])
        grad
        cont       = cont + 1;
        if(abs(dif) < thres)
          return
        end

    end
end
