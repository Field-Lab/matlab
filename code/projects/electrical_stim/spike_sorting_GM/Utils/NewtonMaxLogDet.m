function lambda = NewtonMaxLogDet(input,lambda0,Prods,quad)
% NewtonMaxLogDet finds the hyperparameters lambda needed for the regularized model for the artifact using gradient descent 
% it solves the problem Artifact regularized model.  Actually it solves the problem  min_lambda_i quad'*lambda- log(det(sum(lambda_i*Prods{i}))
% Input:  -input: input structure, 
%         -lambda0: initial lambda estimate (same length as Prods)
%         -quad: quadratic product, quad(i)=Artifact'*Prods{i}*Artifact where  Artifact is the vectorized form of the Artifact computed in the initialization (via quadratic programing)    
%         -Prods: a cell array containing the products D'*D where the D'*D's are the regularization matrix in all the features (time and condition) 
%Output:  -lambda: the estimated hyperparameters
% Gonzalo Mena 06/15
beta   = input.params.initial.Newton.beta;
a      = input.params.initial.Newton.a;
thres  = input.params.initial.Newton.thres;
dif    = Inf;

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

    if(logDetExp(lambdanew,Prods,quad) > val2 || isnan(logDetExp(lambdanew,Prods,quad)))
        t    = beta*t;
        cont = cont+1;
    else
        lambda     = lambdanew;
        [val grad] = logDetExp(lambda,Prods,quad);
        dif        = val0-val;
        val0       = val;
        t          = 1;
        if(isnan(dif))
           5
           break;
        end
        %disp(['Iteration ' num2str(cont) ' finished'])
        %disp(['Diference in function evaluation with previous iteration = ' num2str(dif)])
        %disp(['Gradient = '])
        %grad
        cont       = cont + 1;
        if(abs(dif) < thres)
          return
        end

    end
end
