clear
close all
clc

%cd('/Users/alexandresollaci/Documents/UChicago/RA/Combinatorial growth/combinatorial-growth/simulations/parameter search/')

tic
% 
% options = optimset('Display', 'iter', 'PlotFcns', @optimplotx);
% options.MaxFunEvals = 10000;
% options.MaxIter = 250;
% 
% % Intial guesses
% etaH = 0.2; % NT step size
% etaM = 0.04; % NC step size
% tau = 400; % shape parameter for ideas distribution
% lambda = 2.1; % alpha (cost) has weibull distribution with scale parameter lambda ...
% kappa = 7; % and shape parameter kappa
% xi = 75;   % 1/xi is the fraction of feasible combinations
% 
% x0 = [etaH; etaM; tau; lambda; kappa; xi];

%[sol, fval0] = fminsearch( @(x) MMobj(x), x0, options );

etaH = [.2];
etaM = [.05];
tau = [800, 1000];
lambda = [2];
kappa = [7];
xi = [150, 200, 300];
nrofinv = [800, 1000];

fval0 = 20;

for i = 1:length(etaH)
    for j = 1:length(tau)
        for k = 1:length(lambda)
            for l = 1:length(kappa)
                for m = 1:length(xi)
                    for n = 1:length(etaM)
                        for o = 1:length(nrofinv)
                            x0 = [etaH(i); etaM(n); tau(j); lambda(k); kappa(l); xi(m); nrofinv(o)];
                            fval = MMobj(x0);
                            if fval < fval0
                                fval0 = fval;
                                sol = x0;
                            end
                        end
                    end
                end
            end
        end
    end
end

T = table;
T.etaH = sol(1);
T.etaM = sol(2);
T.tau = sol(3);
T.lambda = sol(4);
T.kappa = sol(5);
T.xi = sol(6);
T.num_inv = sol(7);
T.fval = fval0;

display(T)




toc