function [fStar, fStarSigma] = createResultPredictions( x, f, xStar, theta)

    nData = length(f);
    nPredictions = length(xStar);
    
    nDimensions = size( x, 2);

    R_xx = ones( nData, nData, nDimensions);
    R_xSx = ones( nPredictions, nData, nDimensions);
    R_xSxS = ones( nPredictions, nPredictions, nDimensions);
    
    for ii = 1:nData
        for jj = 1:nData
            for kk = 1:nDimensions

                x_ii_kk = x( ii, kk);
                x_jj_kk = x( jj, kk);
    
                R_xx( ii, jj, kk) = abs( x_ii_kk - x_jj_kk );
            end
        end
    end
    
    for ii = 1:nPredictions
        for jj = 1:nData
            for kk = 1:nDimensions

                x_ii_kk = xStar( ii, kk);
                x_jj_kk = x( jj, kk);
    
                R_xSx( ii, jj, kk) = abs( x_ii_kk - x_jj_kk );
            end
        end
    end

    for ii = 1:nPredictions
        for jj = 1:nPredictions
            for kk = 1:nDimensions

                x_ii_kk = xStar( ii, kk);
                x_jj_kk = xStar( jj, kk);
    
                R_xSxS( ii, jj, kk) = abs( x_ii_kk - x_jj_kk );
            end
        end
    end
    
    mu = theta(1);
    s2 = theta(2).^2;
    ls = theta(3:end-1);
    s2F = theta(end).^2;
    
    K_xSxS = s2F * sqrExpCovMatrix( R_xSxS, ls);
    K_xSx = s2F * sqrExpCovMatrix( R_xSx, ls);
    K_xx = s2F * sqrExpCovMatrix( R_xx, ls) + s2 * eye(nData);

    [fStar, fStarSigma] = gpPrediction( f - mu, K_xSxS, K_xSx, K_xx );
    
    fStar = fStar + mu;
end