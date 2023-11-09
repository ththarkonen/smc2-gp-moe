function [psi] = initializePsi( settings )
    
    J = settings.outerParticleAmount;
    N = settings.nInputDimensions;
    K = settings.numberOfExperts;

    gatingMus = zeros( J, N * K);
    gatingSigmas = zeros( J, N * K);
    
    counter = 1;
    
    for kk = 1:K
        
        mu_kk = sampleGatingMu( J, kk, K, settings);
        
        inds_kk = counter:( counter + N - 1 );
        gatingMus( :, inds_kk) = mu_kk;
        
        for n = 1:N
        
            sigma_kk = sampleGatingSigma( J, kk, K, settings);
            gatingSigmas( :, counter) = sigma_kk;

            counter = counter + 1;
        end
    end
    
    gatingPis = sampleGatingPi( J, [], K, settings);
    psi = [ gatingMus, gatingSigmas, gatingPis];
end