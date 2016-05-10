K_range = -2:0.1:0.1;
C = 0.01;

n = length(K_range);
S = rand(10,1);
W = rand(10,1);
val = zeros(n);

for K1 = 1:n
    for K2 = 1:n
        K = [K_range(K1) K_range(K2)];
        temp = W'*exp(conv(S,K,'same'));
        val(K1,K2) = temp - C*exp(temp);
    end
end

surf(val)