function [innerSmcObj] = singleStepInnerSMC( innerSmcObj, settings)

    K = innerSmcObj.K;
    M = settings.innerParticleAmount;
    
    thetaParticles = innerSmcObj.thetaParticles;
    
    ll_K = innerSmcObj.innerLogLikelihoods_K;
    logPrior_K = innerSmcObj.innerLogPrior_K;
    logTarget = innerSmcObj.innerLogTarget;
    
    innerWeights = innerSmcObj.innerWeights;
    
    Z_tM = innerSmcObj.Z_tM;
    logZ_tM = innerSmcObj.logZ_tM;
    
    mcmcMP_K = innerSmcObj.mcmcMP_K;
    ar = innerSmcObj.ar;
    
    C = innerSmcObj.C;
    R_K = innerSmcObj.R_K;
    
    minMcmcSteps = settings.innerMinMcmcSteps;
    maxMcmcSteps = settings.innerMaxMcmcSteps;
    
    ESS = computeESS( thetaParticles, innerWeights);
    
    kappaVector = innerSmcObj.kappaVector;
    
    kappa = kappaVector( end );
    dKappa = kappaVector( end ) - kappaVector( end - 1 );
        
    max_ll_K = repmat( max(ll_K), M, 1);
    innerWeights_K = exp( dKappa * ( ll_K - max_ll_K ) );

    weightSums = repmat( sum(innerWeights_K), M, 1);
    innerWeights_K = innerWeights_K ./ weightSums;
        
    tempESS = -1;
    [thetaParticles, ~, innerWeights_K, ll_K, logPrior_K, proposalMats, mcmcMP_K] = innerProposalAndResample( thetaParticles, thetaParticles, innerWeights_K, ll_K, logPrior_K, tempESS, ar, mcmcMP_K, settings);
    amountAccepted = zeros( 1, K);

    mcmcDistance = 0;

    flagStopMCMC = false;
    mcmcCounter = 0;
    proposalUpdateCounter = 0;
    
    prevThetaParticles = thetaParticles;
    
    prevLogTarget = logTarget;
    
    while( ~flagStopMCMC || (mcmcCounter < minMcmcSteps) )
            
        mcmcCounter = mcmcCounter + 1;
        proposalUpdateCounter = proposalUpdateCounter + 1;

        [thetaParticles, ll_K, logPrior_K, tempAccepted] = thetaMcmcUpdate( thetaParticles, proposalMats, ll_K, logPrior_K, kappa, C, R_K, settings);
        amountAccepted = amountAccepted + tempAccepted;
            
        acceptanceRate = tempAccepted / M;
        ar = acceptanceRate;
        
        [~, ~, ~, ~, ~, proposalMats, mcmcMP_K] = innerProposalAndResample( thetaParticles, prevThetaParticles, innerWeights_K, ll_K, logPrior_K, ESS, ar, mcmcMP_K, settings);

        ll = sum( ll_K, 2);
        logTarget = ll + sum( logPrior_K, 2);
        
        [flagStopMCMC, mcmcDistance] = stopInnerMCMC( logTarget, prevLogTarget, mcmcDistance, settings);
        
        if( mcmcCounter > maxMcmcSteps )
            flagStopMCMC = true;
        end
    end
        
    ll = sum( ll_K, 2);

    meanKappaL = mean( exp( dKappa * ll ) );
    Z_tM = Z_tM * meanKappaL;
    logZ_tM = logZ_tM + log( meanKappaL );
    
    acceptanceRate = amountAccepted / (M * mcmcCounter * K);
     
    innerSmcObj.thetaParticles = thetaParticles;
    
    innerSmcObj.innerLogLikelihoods = ll;
    innerSmcObj.innerLogLikelihoods_K = ll_K;
    
    innerSmcObj.innerLogPrior_K = logPrior_K;
    innerSmcObj.innerLogTarget = logTarget;
    
    innerSmcObj.innerWeights = innerWeights;
    innerSmcObj.innerWeights_K = innerWeights_K;
    
    innerSmcObj.Z_tM = Z_tM;
    innerSmcObj.logZ_tM = logZ_tM;

    innerSmcObj.Z_tM_hist = [ innerSmcObj.Z_tM_hist, Z_tM];
    innerSmcObj.outerWeight = meanKappaL;
    innerSmcObj.outerLogLikelihood = mean( ll );
    innerSmcObj.likelihoodEvaluations = innerSmcObj.likelihoodEvaluations + mcmcCounter * M;
        
    innerSmcObj.ar = acceptanceRate;
    innerSmcObj.mcmcMP = mcmcMP_K;
        
    innerSmcObj.C = C;
    innerSmcObj.K = K;
    innerSmcObj.R_K = R_K;
end