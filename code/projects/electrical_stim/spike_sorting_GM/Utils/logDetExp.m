function [val grad]=logDetExp(lambda,Prods,quad)

if(max(exp(lambda))==Inf)
    val = Inf;
    grad = Inf*ones(length(lambda),1);
    return
end

sumprod = 0;

for j= 1:length(lambda)
  sumprod=sumprod+exp(lambda(j))*Prods{j};
end


[r p] = chol(sparse(sumprod));

if(p > 0)
    val = Inf;
else
val = dot(exp(lambda),quad)-2*sum(log(diag(r)));
end
if(nargout==1)
    return
end
if(val == Inf)
    grad = Inf*ones(length(lambda),1);
    return
end

in = inv(sparse(sumprod));
for j = 1:length(lambda)
    grad(j) = quad(j)*exp(lambda(j))-trace(exp(lambda(j))*Prods{j}*in);

    
end
grad = grad';
if(p > 0)
    val  = Inf;
    grad = Inf*ones(length(lambda),1);
    grad = grad';
    return
end
