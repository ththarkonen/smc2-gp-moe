function [theta, ll_K, logPrior_K, acceptedAmount] = thetaMcmcUpdate( theta, proposalMats, ll_K, logPrior_K, kappa, C, R_K, settings)

    K = size( C, 2);
    N = settings.nInputDimensions;

    [nParticles, nParameters] = size(theta);
    propThetaParticles = zeros( nParticles, nParameters);
        
    thetaRandInds = randperm( nParticles );
    
    for kk = 1:K
        
        inds_kk = computeThetaInds( kk, K, N);
        
        nParametersPerGP = length(inds_kk);

        randMat = randn( nParametersPerGP, nParticles);
        proposalMat_kk = proposalMats(:,:,kk);
        
        thetaParticles_kk = theta(:,inds_kk);
        thetaRandInds = 1:nParticles;
        p0 = thetaParticles_kk( thetaRandInds, :);
        
        propThetaParticles_kk = ( proposalMat_kk * randMat)' + p0;
        propThetaParticles(:,inds_kk) = propThetaParticles_kk;
    end
    
    accepted = zeros( nParticles, K);
    
    for ind = 1:nParticles
        
        propTheta = propThetaParticles( ind, :);
        
        ll_K_m = ll_K(ind,:);
        logPrior_K_m = logPrior_K( ind, :);
        
        [prop_ll_K, prop_ll_K_m] = gpExpertsLogLikelihood( C, R_K, propTheta, N);
        
        for kk = 1:K

            inds_kk = computeThetaInds( kk, K, N);
            
            propTheta_kk = propTheta(inds_kk);
            logPropPrior_kk = thetaLogPrior( propTheta_kk, kk, K, settings);
            
            ll_kk_m = ll_K_m(kk);
            prop_ll_kk_m = prop_ll_K_m(kk);
        
            logLikelihoodRatio_kk = kappa * prop_ll_kk_m;
            logLikelihoodRatio_kk = logLikelihoodRatio_kk + logPropPrior_kk;
            
            logLikelihoodRatio_kk = logLikelihoodRatio_kk - kappa * ll_kk_m;
            logLikelihoodRatio_kk = logLikelihoodRatio_kk - logPrior_K_m(kk);

            logU = log(rand());

            if( logU < logLikelihoodRatio_kk )
                
                theta( ind, inds_kk) = propTheta_kk;
                ll_K( ind, kk) = prop_ll_kk_m;
                logPrior_K( ind, kk) = logPropPrior_kk;

                accepted( ind, kk) = 1;
            end
        end
    end
    
    acceptedAmount = sum(accepted);
end