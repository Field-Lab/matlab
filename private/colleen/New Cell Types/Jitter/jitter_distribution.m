x = jitter_x_new;
y = jitter_y_new;

results = zeros(16,16);

for i = 1:length(x)
        results(x(i)+9, y(i)+9) = results(x(i)+9, y(i)+9) + 1;
end



