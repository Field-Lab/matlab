function PhiTilde = makePhiTilde(Phi, m, n, lambda)
    PhiTilde = nan(m^2+m*n,m*n);
    for i = 1:m
        PhiTilde((i-1)*m+1:i*m,:) = [zeros(m,(i-1)*n) Phi zeros(m,(m-i)*n)];
    end
    for i = 1:m
        PhiTilde(m^2+(i-1)*n+1:m^2+i*n,:) = [zeros(n,(i-1)*n) sqrt(lambda)*eye(n) zeros(n,(m-i)*n)];
    end
end