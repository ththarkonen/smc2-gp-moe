function [mus] = sampleGatingMu( J, kk, K, settings)

    nInputDimensions = settings.nInputDimensions;
    nGridPoints = K^( 1 / nInputDimensions );
    
    mu0 = linspace( 0, 1, nGridPoints + 2);
    dMu = mu0(2) - mu0(1);
    
    gatingNetworkMuSigma = settings.prior.psi.muSigma * dMu;
    
    mu0_vector = mu0(2:end-1);
    mu0 = mu0_vector;
    
    for ii = 2:nInputDimensions
    
        mu0 = combvec( mu0, mu0_vector);
    end
        
    mu0 = mu0(:,kk);
    mus = zeros( J, nInputDimensions);
    
    for ii = 1:nInputDimensions
        
        mu0_ii = mu0(ii);
        mus(:,ii) = normrnd( mu0_ii, gatingNetworkMuSigma, J, 1);
    end
    
    mus = abs( mus );
end