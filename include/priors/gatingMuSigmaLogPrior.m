function [logP] = gatingMuSigmaLogPrior( mu, mu0_d, sigma, settings)
    
    N = settings.nInputDimensions;
    K = settings.numberOfExperts;

    nGridPoints = K^( 1 / N );
    
    mu0 = linspace( 0, 1, nGridPoints + 2);
    dMu = mu0(2) - mu0(1);
    
    muSigma = settings.prior.psi.muSigma * dMu;
    sigmaSigma = settings.prior.psi.sigmaSigma * dMu;

    p = ( mu >= 0 ) & ( mu <= 1 ) & ( sigma >= 0 );
    p = double( p );
    
    p = p .* normpdf( mu, mu0_d, muSigma);
    p = p .* normpdf( sigma, 0, sigmaSigma);
    
    logP = log( p );
end