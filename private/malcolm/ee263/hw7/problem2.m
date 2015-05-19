% problem 2b

% define node-incidence matrix B
B = zeros(m,n);
counter = 1;
for i = 1:n
    for j = i:n
        if A(i,j)
            B(counter,i) = sign(j-i);
            B(counter,j) = sign(i-j);
            counter = counter+1;
        end
    end
end

% one solution for the optimal x and y is the right singular vectors
% of B with the smallest singular values
[u, s, v] = svd(B);
xopt = v(:,end-2);
yopt = v(:,end-1);

% calculate optimal value of J
X = [xopt yopt];
Jopt = trace(X'*(B')*B*X);

% draw circle graph
figure();
[i,j] = find(triu(A));
plot(x_circ , y_circ , 'o'); 
hold on;
for k = 1:m
    plot(x_circ([i(k) j(k)]) , ...
    y_circ([i(k) j(k)]) , '-');
end
hold off;
axis equal;
xlabel('x'); ylabel('y');

% draw optimal graph
figure();
[i,j] = find(triu(A));
plot(xopt , yopt , 'o'); 
hold on;
for k = 1:m
    plot(xopt([i(k) j(k)]) , ...
    yopt([i(k) j(k)]) , '-');
end
hold off;
axis equal;
xlabel('x'); ylabel('y');