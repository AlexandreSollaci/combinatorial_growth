clear 
close all

cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial-growth/simulations/simulate model')


x = zeros(10000,1);
y = zeros(10000,1);
tau = 400;
m = 50;
lambda = 2;
kappa = 7;

for i = 1:10000
	[x(i), y(i)] = modeldraws(m, tau, lambda, kappa);
end

[x_dens, x_grid] = ksdensity(x);
[y_dens, y_grid] = ksdensity(y);

figure(1)
plot(x_grid, x_dens)
title('logistic')

figure(2)
plot(y_grid, y_dens)
title('weibull')
