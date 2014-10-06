cd /snle/home.1173742073.11117/lhruby/MATLAB_code/

load prob6_3_data.mat

%% computing log probabilities of most likely t-step paths l_it*

l_star = zeros(26, length(o));
s_star = zeros(length(o), 1);

%base case
for i = 1:26
    l_star(i, 1) = log(pi(i)) + log(b(i,o(1)+1));
end

%rest of cases
for t = 2:length(o)
    for j = 1:26
        l_star(j, t) = max(l_star(:,t-1) + log(a(:,j))) + log(b(j, o(t)+1));
    end
end


%end case
s_star(length(o)) = find(l_star(:, length(o)) == max(l_star(:, length(o))));

%back-tracking
for t = length(o)-1:-1:1
    temp = l_star(:, t) + log(a(:,s_star(t+1)));
    s_star(t) = find(temp == max(temp));
end


break

%removing repeats
for i = length(s_star):-1:2
    i
    if s_star(i) == s_star(i-1);
        s_star(i) = [];
    end
end

disp(['' 96+s_star])