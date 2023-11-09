function [ smc2Object, acceptedAmount] = psiMcmcUpdate( smc2Object )
    
    settings = smc2Object.settings;
    
    K = smc2Object.K;
    N = settings.nInputDimensions;

    data = smc2Object.data;

    kappaVector = smc2Object.kappaVector;
    proposalMat = smc2Object.proposalMat;

    psiParticles = smc2Object.psiParticles;
    innerSMCs = smc2Object.innerSmcObjects;
    ll = smc2Object.outerLL;
    
    logTarget = smc2Object.logTarget;
    gatingProbabilities = smc2Object.gatingProbabilities;

    [nParticles, ~] = size( psiParticles );
    
    NK = N * K;
    
    muInds = 1:( NK );
    sigmaInds = ( NK + 1 ):2*NK;
    piInds = ( 2*NK + 1 ):(2*NK + K);
    
    muSigmaInds = [muInds, sigmaInds];
    
    muSigmaProposalMat = proposalMat( muSigmaInds, muSigmaInds);
    piProposalMat = proposalMat( piInds, piInds);
    
    nMuSigma = length( muSigmaInds );
    nPi = length( piInds );
    
    randGatingInds = randperm( nParticles );
    randGatingInds = 1:nParticles;
    
    gatingMuSigma_K = psiParticles( :, muSigmaInds);

    gatingMu_K = psiParticles( :, muInds);
    gatingSigma_K = psiParticles( :, sigmaInds);
    gatingPi_K = psiParticles( :, piInds);
    
    randMatMuSigma = randn( nMuSigma, nParticles);
    randMatPi = randn( nPi, nParticles);
    
    gatingMuSigma_K0 = gatingMuSigma_K( randGatingInds, :);
    gatingPi_K0 = gatingPi_K( randGatingInds, :);

    propGatingMuSigma_K = ( muSigmaProposalMat * randMatMuSigma )' + gatingMuSigma_K0;
    
%     logPropGatingMu_K = logPropGatingMuSigma_K( :, muInds);
%     logPropGatingSigma_K = logPropGatingMuSigma_K( :, sigmaInds);
    logPropGatingPi_K = ( piProposalMat * randMatPi)' + log( gatingPi_K0 );
    
    propGatingMu_K = propGatingMuSigma_K( :, muInds);
    propGatingSigma_K = propGatingMuSigma_K( :, sigmaInds);
    propGatingPi_K = exp( logPropGatingPi_K );
    
    propPsiParticles = [ propGatingMu_K, propGatingSigma_K, propGatingPi_K];
    
    accepted = zeros( nParticles, 1);
    
    smc2Object.propPsiParticles = propPsiParticles;

    kappa = kappaVector(end);
    update = waitbarParfor( nParticles, " ", 'kappa', kappa);

    parfor ind = 1:nParticles
        
        logZ_tM_j = innerSMCs{ind}.logZ_tM;
        logGatingPrior_j = innerSMCs{ind}.gatingLogPrior;
        
        propMu_j = propGatingMu_K( ind, :);
        propSigma_j = propGatingSigma_K( ind, :);
        propPi_j = propGatingPi_K( ind, :);

        [propPart_jj, propGateProb_j, ~, errorFlag] = partition( data, propMu_j, propSigma_j, propPi_j);
        
        
        if( errorFlag )
            prop_logZ_tM_j = -Inf;
        else

            tempInnerSmcObj = struct();
            
            tempInnerSmcObj.K = K;
            tempInnerSmcObj.C = propPart_jj.C;
            tempInnerSmcObj.R_K = propPart_jj.R_K;
            tempInnerSmcObj.kappaVector = kappaVector;
            
            propLogPrior = psiLogPrior( propMu_j, propSigma_j, propPi_j, settings);
            logPropGatingPrior_j = propLogPrior + log( propGateProb_j );

            prop_innerSmcObject_j = innerSMC( tempInnerSmcObj, settings);
            prop_innerSmcObject_j.gatingLogPrior = logPropGatingPrior_j;
            
            prop_logZ_tM_j = prop_innerSmcObject_j.logZ_tM;
        end
        
        logLikelihoodRatio = prop_logZ_tM_j;
        logLikelihoodRatio = logLikelihoodRatio - logZ_tM_j;

        logLikelihoodRatio = logLikelihoodRatio + logPropGatingPrior_j;
        logLikelihoodRatio = logLikelihoodRatio - logGatingPrior_j;
        
        logU = log(rand());
        
        if( logU <= logLikelihoodRatio && isfinite(logLikelihoodRatio) )
            
            propPi_j = propPi_j / sum( propPi_j );
            
            psiParticles(ind,:) = [propMu_j, propSigma_j, propPi_j];
            ll(ind) = prop_innerSmcObject_j.outerLogLikelihood;
            gatingProbabilities(ind) = propGateProb_j;
            innerSMCs{ind} = prop_innerSmcObject_j;
            
            logTarget(ind) = prop_logZ_tM_j + logPropGatingPrior_j;
            
            accepted(ind) = 1;
        end
        
        update();
    end
    
    smc2Object.psiParticles = psiParticles;
    smc2Object.outerLL = ll;
    smc2Object.gatingProbabilities = gatingProbabilities;
    smc2Object.innerSmcObjects = innerSMCs;
    smc2Object.logTarget_K = logTarget;
    
    acceptedAmount = sum(accepted);
end
