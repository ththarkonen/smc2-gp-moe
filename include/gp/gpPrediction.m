function [fStar, fSigma, fCov] = gpPrediction( f, K_xSxS, K_xSx, K_xx )

    fStar = K_xSx * inv(K_xx) * f;
    fCov  = K_xSxS - K_xSx * inv(K_xx) * K_xSx';
    
    fSigma = sqrt( diag(fCov) );
end