function [fStar, fStarSigma] = createPredictions( x, f, theta)

    nPredictions = length(x);
    R_xS = ones( nPredictions, nPredictions);

    for ii = 1:nPredictions
        for jj = 1:nPredictions

            xII = x(ii);
            xJJ = x(jj);

            R_xS(ii,jj) = abs( xII - xJJ );
        end
    end
    
    mu = theta(1);
    s2 = theta(2).^2;
    ls = theta(3);
    s2F = theta(4).^2;

    K_xSxS = s2F * sqrExpCovMatrix( R_xS, ls);
    K_xSx = K_xSxS;
    K_xx = K_xSxS + s2 * eye(nPredictions);

    [fStar, fStarSigma, fCov] = gpPrediction( f - mu, K_xSxS, K_xSx, K_xx );
    
    fStar = fStar + mu;
    fStarSigma = sqrt( diag( fCov ) );
end