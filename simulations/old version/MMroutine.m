clear all
close all
clc
tic

% lb = [0; 0; 0 ; 1]; % lower bound
% options = optimoptions(@fmincon,'Algorithm','interior-point');
% options.MaxFunEvals = 10000;
% options.MaxIter = 500;

options = optimset('Display', 'iter', 'PlotFcns', @optimplotx);
options.MaxFunEvals = 10000;
options.MaxIter = 250;

% % tau = [80; 100; 120; 150; 200];
% % etaH = [.1; .15; .2; .25];
% % xi = [1; 1.5; 2];
% % In this range, the minimum value for F was 0.028, with
% % taustar = 125, etaHstar = 0.1, etaMstar = 0.05 and xistar = 2
% % (approximately).


% New set of parameters.
tau = [115; 120; 125; 130; 135];
etaH = [.05; .07; .1; .12];
xi = [1.8; 2; 2.2];

etaHstar = [];
etaMstar = [];
taustar = [];
xistar = [];
funval = [];

for i = 1:length(tau)
    for j = 1:length(etaH)
        for k = 1:length(xi)

            x0 = [etaH(j), etaH(j)-0.05, tau(i), xi(k)]; %etaH, etaM, tau, xi

            [x, fval] = fminsearch( @(x) MMobj(x), x0, options );
            % x = fmincon( @(y) MMobj(y) , y0 , [], [], [], [], lb, [] , [], options);
            x
            etaHstar(i,j,k) = x(1);
            etaMstar(i,j,k) = x(2);
            taustar(i,j,k) = x(3);
            xistar(i,j,k) = x(4);
            funval(i,j,k) = fval;
        end
    end
end


toc