function [innerSMC] = innerSMC( innerSMC, settings)
    
    N = settings.nInputDimensions;
    M = settings.innerParticleAmount;
    K = settings.numberOfExperts;
    
    NK = N * K;

    kappaVector = innerSMC.kappaVector;
    nKappa = length( kappaVector );
    kappa = kappaVector(1);
    
    C = innerSMC.C;
    R_K = innerSMC.R_K;
    
    innerSMC = {};
    
    gpMeans = zeros( M, K);
    gpNoiseSigmas = zeros( M, K);
    gpLengthScales = zeros( M, NK);
    gpSigmalSigmas = zeros( M, K);
    
    counter = 1;
    
    for kk = 1:K

        gpMu_kk = sampleExpertMean( M, kk, K, settings);
        gpSigma_kk = sampleExpertNoiseSigma( M, kk, K, settings);
        gpSignalSigma_kk = sampleExpertSignalSigma( M, kk, K, settings);

        gpMeans(:,kk) = gpMu_kk;
        gpNoiseSigmas(:,kk) = gpSigma_kk;
        gpSigmalSigmas(:,kk) = gpSignalSigma_kk;
        
        for n = 1:N
            
            gpLengthScale_kk_n = sampleExpertLengthScale( M, kk, K, settings); 
            gpLengthScales( :, counter) = gpLengthScale_kk_n;
            
            counter = counter + 1;
        end 
    end
    
    theta = [ gpMeans, gpNoiseSigmas, gpLengthScales, gpSigmalSigmas];
    
    ll = zeros( M, 1);
    ll_K = zeros( M, K);
    logPrior_K = zeros( M, K);
    
    logTarget = zeros( M, 1);
    
    for m = 1:M
       
        theta_m = theta( m, :);
        [ ll_m, ll_m_K] = gpExpertsLogLikelihood( C, R_K, theta_m, N);
        
        for kk = 1:K

            inds_kk = computeThetaInds( kk, K, N);
            theta_kk = theta_m( inds_kk );
                
            logPrior_K( m, kk) = thetaLogPrior( theta_kk, kk, K, settings);
        end
        
        ll(m) = ll_m;
        ll_K(m,:) = ll_m_K;
        
        thetaTotalLogPrior = sum( logPrior_K(m,:), 2);
        logTarget(m) = kappa * ll_m + thetaTotalLogPrior;
    end

    innerWeights = ones( M, 1) / M;
    innerWeights_K = ones( M, K) / M;
    
    Z_tM = 1;
    logZ_tM = 0;
    likelihoodEvaluations = M;
    
    mcmcMP_K = settings.initialInnerProposalMultiplier * ones( 1, K);
    targetAR = settings.innerTargetAcceptanceRate;
    
    mcmcAcceptanceRateHist = zeros( nKappa, K);
    mcmcAcceptanceRateHist(1,:) = targetAR * ones( 1, K);
    
    minMcmcSteps = settings.innerMinMcmcSteps;
    maxMcmcSteps = settings.innerMaxMcmcSteps;
        
    innerSMC.thetaParticles = theta;
    innerSMC.propThetaParticles = [];

    innerSMC.innerLogLikelihoods = ll;
    innerSMC.innerLogLikelihoods_K = ll_K;

    innerSMC.innerLogPrior_K = logPrior_K;
    innerSMC.innerLogTarget = logTarget;

    innerSMC.innerWeights = innerWeights;
    innerSMC.innerWeights_K = innerWeights_K;

    innerSMC.Z_tM = Z_tM;
    innerSMC.logZ_tM = logZ_tM;
    innerSMC.Z_tM_hist = [Z_tM];
    innerSMC.outerWeight = 1.0;
    innerSMC.outerLogLikelihood = mean( ll );

    innerSMC.ar = targetAR;
    innerSMC.mcmcMP_K = mcmcMP_K;
    innerSMC.likelihoodEvaluations = likelihoodEvaluations;

    innerSMC.C = C;
    innerSMC.R_K = R_K;
    innerSMC.K = K;

    if( nKappa < 2 )
        return;
    end
    
    ESS = computeESS( theta, innerWeights);
    
    for ii = 2:nKappa
        
        ar = mcmcAcceptanceRateHist( ii - 1, :);
        
        kappa = kappaVector( ii );
        dKappa = kappaVector( ii ) - kappaVector( ii - 1 );
        
        max_ll_K = repmat( max(ll_K), M, 1);
        innerWeights_K = exp( dKappa * ( ll_K - max_ll_K ) );
        
        weightSums = repmat( sum(innerWeights_K), M, 1);
        innerWeights_K = innerWeights_K ./ weightSums;
        
        tempESS = -1;
        [theta, ~, innerWeights_K, ll_K, logPrior_K, proposalMats, mcmcMP_K] = innerProposalAndResample( theta, theta, innerWeights_K, ll_K, logPrior_K, tempESS, ar, mcmcMP_K, settings);
        
        ll = sum( ll_K, 2);
        logTarget = kappa * ll + sum( logPrior_K, 2);
        
        amountAccepted = zeros( 1, K);
            
        mcmcDistance = 0;

        flagStopMCMC = false;
        mcmcCounter = 0;
        proposalUpdateCounter = 0;
        
        prevTheta = theta;
        
        prevLogTarget = logTarget;

        while( ~flagStopMCMC || (mcmcCounter < minMcmcSteps) )
            
            mcmcCounter = mcmcCounter + 1;
            proposalUpdateCounter = proposalUpdateCounter + 1;

            [theta, ll_K, logPrior_K, tempAccepted] = thetaMcmcUpdate( theta, proposalMats, ll_K, logPrior_K, kappa, C, R_K, settings);
            amountAccepted = amountAccepted + tempAccepted;
            
            acceptanceRate = tempAccepted / (M);
            ar = acceptanceRate;
            
            [~, ~, ~, ~, ~, proposalMats, mcmcMP_K] = innerProposalAndResample( theta, prevTheta, innerWeights_K, ll_K, logPrior_K, ESS, ar, mcmcMP_K, settings);

            ll = sum( ll_K, 2);
            expertsLogPrior = sum( logPrior_K, 2);
            
            logTarget = kappa * ll + expertsLogPrior;
        
            [flagStopMCMC, mcmcDistance] = stopInnerMCMC( logTarget, prevLogTarget, mcmcDistance, settings);
            
            ESS = computeESS( theta, innerWeights);
            
            if( mcmcCounter > maxMcmcSteps )
                flagStopMCMC = true;
            end
        end
        
        likelihoodEvaluations = likelihoodEvaluations + mcmcCounter * M;
        ll = sum( ll_K, 2);

        meanKappaL = mean( exp( dKappa * ll ) );
        Z_tM = Z_tM * meanKappaL;
        logZ_tM = logZ_tM + log( meanKappaL );
        
        acceptanceRate = amountAccepted / (M * mcmcCounter);
        mcmcAcceptanceRateHist(ii,:) = acceptanceRate;
    end
     
    innerSMC.thetaParticles = theta;
    
    innerSMC.innerLogLikelihoods = ll;
    innerSMC.innerLogLikelihoods_K = ll_K;
    
    innerSMC.innerLogPrior_K = logPrior_K;
    innerSMC.innerLogTarget = logTarget;
    
    innerSMC.innerWeights = innerWeights;
    innerSMC.innerWeights_K = innerWeights_K;
    
    innerSMC.Z_tM = Z_tM;
    innerSMC.logZ_tM = logZ_tM;
    innerSMC.Z_tM_hist = [ innerSMC.Z_tM_hist, Z_tM];
    innerSMC.outerWeight = meanKappaL;
    innerSMC.outerLogLikelihood = mean( ll );
        
    innerSMC.ar = acceptanceRate;
    innerSMC.mcmcMP = mcmcMP_K;
    innerSMC.likelihoodEvaluations = likelihoodEvaluations;
        
    innerSMC.C = C;
    innerSMC.K = K;
    innerSMC.R_K = R_K;
    innerSMC.kappaVector = kappaVector;
end