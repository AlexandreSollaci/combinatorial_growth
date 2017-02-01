function [x,z] = modeldraws(m, tau, lambda, kappa, seed)
    
    % [X, Z] = modeldraws(M, TAU, LAMBDA, KAPPA, SEED).
    %
    % This function provides random draws from (1) a truncated logistic
    % distribution with location parameter M and shape parameter TAU; and
    % (2) a weibull distribution with scale parameter LAMBDA and shape 
    % parameter KAPPA.
    %
    % All parameters must be scalar.
    %
    % SEED is an optinal seed for the random number generator.


    if nargin == 5
    rng(seed)
    end
    
    if (isscalar(m) ~= 1 || isscalar(tau) ~= 1 || isscalar(lambda) ~= 1)
        display('Parameters must be scalars!')
    else
        y = rand(2,1); 
        x = tau*log( (1 + y(1)*exp(m/tau)) / (1-y(1)) );                    % Inverse of G
        z = lambda*(- log(1 - y(2)) )^(1/kappa) ;                           % Inverse of weibull CDF
    end
end