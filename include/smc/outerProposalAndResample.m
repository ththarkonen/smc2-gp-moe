function [ psi, innerSMCs, ll, weights, proposalMat, mcmcMP] = outerProposalAndResample( psi, innerSMCs, ll, weights, ESS, ar, mcmcMP, settings)
    
    minESS = settings.outerMinimumEffectiveSampleSize;
    targetAR = settings.outerTargetAcceptanceRate;
    
    N = settings.nInputDimensions;
    K = settings.numberOfExperts;

    [nParticles, nParameters] = size( psi );
    
    NK = N * K;
    piInds = ( 2*NK + 1 ):( 2*NK + K );
    
    weights( weights == 0 ) = eps;
    psi( psi == 0 ) = eps;
    
    logPsi = psi;
    logPsi( :, piInds) = log( logPsi( :, piInds) );
    
    Kmat = weightedcov( logPsi, weights);
    Kmat = Kmat + eps * eye( nParameters, nParameters);
    
    invD = inv( sqrt( diag(diag(Kmat)) ) );
    corMat = invD * Kmat * invD;
    
    if( ESS < minESS )
        
        resampleInds = resampleResidual(weights);
        
        psi = psi( resampleInds, :);
        innerSMCs = innerSMCs( resampleInds );
        
        ll = ll( resampleInds );
        weights = ones( nParticles, 1) / nParticles;
        
        uniqueParticles = unique( psi, 'rows' );
        nUniques = size( uniqueParticles, 1);
        
        disp(['Resampling outer particles...'])
        disp([ num2str(nUniques), ' unique outer particles.'])
    end
    
    logpMedian = median( logPsi );
    
    absoluteDeviations = abs( logPsi - repmat( logpMedian, nParticles, 1) );
    medianAbsDevs = median( absoluteDeviations );
    
    Kmat = 1.4826^2 * ( medianAbsDevs' * medianAbsDevs ) .* corMat;
    
    if( isfinite(ar) )    
        exponent = -5 * ( targetAR - ar );
        mcmcMP = 2^exponent * mcmcMP;
    end
    
    covMH = mcmcMP * 2.38^2 / nParameters * Kmat;
    
    try
        proposalMat = chol( covMH );
    catch
        proposalMat = sqrt(mcmcMP) * eye( nParameters );
        disp('Fallback outer proposal');
    end
end