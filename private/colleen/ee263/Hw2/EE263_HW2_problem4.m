deltaT = 0.001;
T = max(beta) + deltaT:deltaT:5;

for i = T;
    U = diag((beta-i)./gamma) + repmat(alpha, 1, n).* F;
    x = U\(-Cext);
    if(max(x)<=xmax & min(x)>=xmin & sum(a.*x)<=Amax)
        break
    end
end