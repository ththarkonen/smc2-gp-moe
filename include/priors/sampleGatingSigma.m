function [sigmas] = sampleGatingSigma( J, kk, K, settings)
    
    nInputDimensions = settings.nInputDimensions;
    nGridPoints = K^( 1 / nInputDimensions );
    
    mu0 = linspace( 0, 1, nGridPoints + 2);
    dMu = mu0(2) - mu0(1);
    
    gatingSigmaSigma = settings.prior.psi.sigmaSigma * dMu;
    sigmas = gatingSigmaSigma * randn( J, 1);
    sigmas = abs( sigmas );
end