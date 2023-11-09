function [ theta, p0, w_K, ll_K, logPrior_K, proposalMats, mcmcMP_K] = innerProposalAndResample( theta, p0, w_K, ll_K, logPrior_K, ESS, ar, mcmcMP_K, settings)
    
    N = settings.nInputDimensions;
    K = settings.numberOfExperts;
    
    minESS = settings.innerMinimumEffectiveSampleSize;
    targetAR = settings.innerTargetAcceptanceRate;

    [ nParticles, nParameters] = size( theta );
    nParameterPerGP = nParameters / K;
    
    w_K( w_K == 0 ) = eps;
    proposalMats = zeros( nParameterPerGP, nParameterPerGP, K);

    if( isfinite( sum(ar) ) )    
        exponent = -5 * ( targetAR - ar );
        mcmcMP_K = 2.^exponent .* mcmcMP_K;
    end
    
    for kk = 1:K

        inds_kk = computeThetaInds( kk, K, N);
        
        weights_kk = w_K(:,kk);
        
        p_kk = theta( :, inds_kk);
        p0_kk = p0( :, inds_kk);
        ll_kk  = ll_K(:,kk);
        logPrior_kk = logPrior_K(:,kk);
        
        mcmcMP_kk = mcmcMP_K(kk);
        
        if( ESS < minESS )
            resampleInds = resampleResidual( weights_kk );

            p_kk = p_kk( resampleInds, :);
            ll_kk = ll_kk( resampleInds );
            logPrior_kk = logPrior_kk( resampleInds );

            weights_kk = ones( nParticles, 1) / nParticles;
        end
    
        try
            K_mat = weightedcov( p_kk, weights_kk);
        catch
            K_mat = eye( nParameterPerGP, nParameterPerGP);
        end

        invD = inv( sqrt( diag(diag(K_mat)) ) );
        corMat = invD * K_mat * invD;

        pMedian = median(p_kk);
        absoluteDeviations = abs( p_kk - repmat( pMedian, nParticles, 1) );
        medianAbsDevs = median( absoluteDeviations );

        K_mat = 1.4826^2 * ( medianAbsDevs' * medianAbsDevs ) .* corMat;

        covMH = mcmcMP_kk * 2.38^2 / nParameterPerGP * K_mat;

        try
            proposalMats(:,:,kk) = chol(covMH);
        catch
            I = eye( nParameterPerGP, nParameterPerGP);
            proposalMats(:,:,kk) = sqrt( mcmcMP_kk ) * I;
        end
        
        w_K(:,kk) = weights_kk;
        
        theta( :, inds_kk) = p_kk;
        p0( :, inds_kk) = p0_kk;
        ll_K(:,kk) = ll_kk;
        logPrior_K(:,kk) = logPrior_kk;
    end
end