clear all
close all
clc

tau = 100:100:500;
etaH = .05:.1:.45;
lambda = linspace(.2, 2, 5);
kappa = linspace(.02, .2, 5);

fstar = 5;
iter = 1;

for hh = 1:length(etaH)
    for tt = 1:length(tau);
        for ll = 1:length(lambda);
            for kk = 1:length(kappa);
                
                k = rand(1);
                etaM = k*etaH(hh);
                x0 = [etaH(hh), etaM, tau(tt), 1, lambda(ll), kappa(kk)];
                fval = MMobj(x0);

                if fval < fstar;
                    fstar = fval;
                    xstar = x0;
                end

                display(['iteration number = ' num2str(iter)])
                display(['current function value = ' num2str(fval)])
                display(['current arguments = ' num2str(x0)])
                display(['current minimum function value = ' num2str(fstar)])
                display(['argmin of function = ' num2str(xstar)])
                display('################## END OF ITERATION ##################')
                iter = iter + 1;
            end
        end
    end
end

