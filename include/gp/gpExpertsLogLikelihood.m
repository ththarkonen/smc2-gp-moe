function [ll, ll_K] = gpExpertsLogLikelihood( C, R_K, theta, nInputDimensions)

    K = size( C, 2);
    
    ll = 0;
    ll_K = zeros( 1, K);
    
    for kk = 1:K
        
        inds_kk = computeThetaInds( kk, K, nInputDimensions);
        
        XY_kk = C{kk};
        
        if( isempty(XY_kk) )
            continue;
        end
        
        Y_k = XY_kk( :, end);
        
        R_kk = R_K{kk};
        theta_kk = theta(inds_kk);
        
        ll_kk = computeGpLogLikelihood( Y_k, R_kk, theta_kk);

        ll_K(kk) = ll_kk;
        ll = ll + ll_kk;
    end
end