function [logP] = psiLogPrior( mus, sigmas, pi_kk, settings)
    
    alpha = settings.concentrationParameter;
    N = settings.nInputDimensions;
    K = settings.numberOfExperts;
    
    alphaVector = (alpha / K) * ones( 1, K);
    nGridPoints = K^( 1 / N );
    
    mu0 = linspace( 0, 1, nGridPoints + 2);
    
    mu0_vector = mu0(2:end-1);
    mu0 = mu0_vector;
    
    for ii = 2:N
    
        mu0 = combvec( mu0, mu0_vector);
    end
    
    logPi_kk = log( pi_kk );
    logP = sum( alphaVector .* logPi_kk ) - sum( pi_kk );
    
    counter = 1;
    
    for kk = 1:K
        for ii = 1:N
            
            mu_kk = mus(counter);
            sigma_kk = sigmas(counter);
            mu0_ii_kk = mu0( ii, kk);

            logPrior_ii_kk = gatingMuSigmaLogPrior( mu_kk, mu0_ii_kk, sigma_kk, settings);
            logP = logP + logPrior_ii_kk;
            
            counter = counter + 1;
        end
    end
end