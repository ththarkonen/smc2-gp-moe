function [ theta, weights, ll_K] = resampleInner( theta, weights, ll_K, settings)

    K = settings.numberOfExperts;
    N = settings.nInputDimensions;

    [nParticles, ~] = size( theta );
    weights( weights == 0 ) = eps;

    for kk = 1:K
        
        weights_kk = weights( :, kk);
        ll_kk  = ll_K( :, kk);
        
        inds_kk = computeThetaInds( kk, K, N);
        resampleInds_kk = resampleResidual( weights_kk );
        
        theta_kk = theta( resampleInds_kk, inds_kk);
        ll_kk = ll_kk( resampleInds_kk, :);

        weights_kk = ones( nParticles, 1) / nParticles;
        
        theta( :, inds_kk) = theta_kk;
        weights( :, kk) = weights_kk;
        ll_K( :, kk) = ll_kk;
    end
end