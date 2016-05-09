x = 1;
y = glm_SU_optimizationfunction(x, rand(100,10), ones(10, 1), cumsum(round(10*rand(10,1))), 1, ones(1, 100));
plot(x,y)