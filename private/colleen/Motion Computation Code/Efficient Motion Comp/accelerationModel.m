% acceleration
noise = random('Normal',0,1,[1000,1])
x_raw = linspace(1,10,1000);

y = (x_raw'+ noise).^2;
figure;
plot(y,x_raw', '.')
set(gcf, 'color', 'w')

extra_noise = 10*rand([500,1])
hold on
y = y(1 : 2 : end);  % => 1 4 7 10
plot(y,extra_noise, '.')
axis off
