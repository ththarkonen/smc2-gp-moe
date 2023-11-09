function [logP] = thetaLogPrior( theta, kk, K, settings)
    
    nInputDimensions = settings.nInputDimensions;
    lengthScaleInds = 3:(3 + nInputDimensions - 1);
    
    maxY = settings.prior.theta.maxY;
    noiseSigmaSigma = settings.prior.theta.noiseSigmaSigma;
    lengthScaleSigma = settings.prior.theta.lengthScaleSigma;
    signalSigmaSigma = settings.prior.theta.signalSigmaSigma;

    gpMean = theta(1);
    noiseSigma = theta(2);
    lengthScales = theta( lengthScaleInds );
    signalSigma = theta(end);
    
    logP = ( 0 <= gpMean ) * ( gpMean <= maxY );
    logP = logP * ( noiseSigma > 0 );
    logP = logP * all( lengthScales > 0);
    logP = logP * ( signalSigma > 0 );
    logP = log( logP );

    if( isinf( logP ) )
        return;
    end
    
    logP = logP - 0.5 * ( noiseSigma / noiseSigmaSigma ).^2;
    logP = logP - 0.5 * log( noiseSigmaSigma );
    logP = logP - 0.5 * log( 2*pi );
    
    logP = logP - 0.5 * ( signalSigma / signalSigmaSigma ).^2;
    logP = logP - 0.5 * log( signalSigmaSigma );
    logP = logP - 0.5 * log( 2*pi );
    
    for ii = 1:nInputDimensions
    
        lengthScale_ii = lengthScales(ii);
    
        logP = logP - 0.5 * ( lengthScale_ii / lengthScaleSigma ).^2;
        logP = logP - 0.5 * log( lengthScaleSigma );
        logP = logP - 0.5 * log( 2*pi );
    end
end