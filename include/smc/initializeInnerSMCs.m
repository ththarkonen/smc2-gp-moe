function [innerSmcs, outerLogLikelihoods, gatingProbabilities, logPosterior] = initializeInnerSMCs( data, psi, settings)
    
    J = settings.outerParticleAmount;
    N = settings.nInputDimensions;
    K = settings.numberOfExperts;
    
    innerSmcs = cell( J, 1);
    
    outerLogLikelihoods = ones( J, 1);
    gatingProbabilities = zeros( J, 1);
    logPosterior = -Inf( J, 1);
    
    NK = N * K;
    
    update = progressInitialization( J, " ");
    parfor jj = 1:J

        muInds = 1:NK;
        sigmaInds = ( NK + 1 ):2*NK;
        piInds = ( 2*NK + 1 ):(2*NK + K);

        mu_jj = psi( jj, muInds);
        sigma_jj = psi( jj, sigmaInds);
        pi_jj = psi( jj, piInds);

        [ partition_jj, gp_jj, ~, error] = partition( data, mu_jj, sigma_jj, pi_jj);

        if( error )
            ll_j = -Inf;
            w_j = 0.0;
            gp_jj = 1;

            innerSmc_jj.C = [];
            innerSmc_jj.K = K;
            innerSmc_jj.gatingLogPrior = -Inf;

            innerSmcs{jj} = innerSmc_jj;
        else

            tempInnerSmcObj = struct();
            tempInnerSmcObj.K = K;
            tempInnerSmcObj.C = partition_jj.C;
            tempInnerSmcObj.R_K = partition_jj.R_K;
            tempInnerSmcObj.kappaVector = [0.0];
            
            logGatingProbability_jj = log( gp_jj );
            logPrior_jj = psiLogPrior( mu_jj, sigma_jj, pi_jj, settings);
            logPrior_jj = logPrior_jj + logGatingProbability_jj;
            
            innerSmc_jj = innerSMC( tempInnerSmcObj, settings);
            
            innerSmc_jj.gatingLogPrior = logPrior_jj;
            innerSmcs{jj} = innerSmc_jj;

            ll_j = innerSmc_jj.outerLogLikelihood;
            
            logPosterior(jj) = log( innerSmc_jj.Z_tM ) + logPrior_jj;
        end

        outerLogLikelihoods(jj) = ll_j;
        gatingProbabilities(jj) = gp_jj;

        update();
    end
end