function [ data, settings] = preprocessNasa( data )

    % Include dependencies to path
    includeFolders = genpath('include');
    addpath( includeFolders );

    nInputDimensions = size( data.x, 2);
    
    for ii = 1:nInputDimensions
        
        x_ii = data.x(:,ii);
        
        minX = min( x_ii );
        maxX = max( x_ii );
    
        data.x(:,ii) = ( x_ii - minX ) / ( maxX - minX );
    end
    
    data.R = computeDistanceMatrix( data.x );
    
    meanY = mean( data.y );

    data.y = data.y - meanY;
    data.y = data.y / std( data.y );
    data.y = data.y - min( data.y );

    data.x = data.x(1:10:end,:);
    data.y = data.y(1:10:end);

    maxY = max( data.y );
        
    settings = struct();
    settings.nInputDimensions = nInputDimensions;
    settings.numberOfExperts = 27;
        
    settings.prior.psi.muSigma = 0.05;
    settings.prior.psi.sigmaSigma = 0.01;
    settings.prior.psi.concentrationParameter = 0.5 * settings.numberOfExperts;

    settings.prior.theta.maxY = maxY;
    settings.prior.theta.noiseSigmaSigma = 0.25 * maxY;
    settings.prior.theta.lengthScaleSigma = 0.125;
    settings.prior.theta.signalSigmaSigma = 0.25 * maxY;
    
    settings.outerParticleAmount = 100;
    settings.innerParticleAmount = 20;
    
    settings.outerMinimumEffectiveSampleSize = 1.1 * settings.outerParticleAmount;
    settings.innerMinimumEffectiveSampleSize = 1.1 * settings.innerParticleAmount;
    
    settings.outerMinMcmcSteps = 3;
    settings.outerMaxMcmcSteps = 10;

    settings.innerMinMcmcSteps = 3;
    settings.innerMaxMcmcSteps = 30;

    settings.outerTargetAcceptanceRate = 0.15;
    settings.innerTargetAcceptanceRate = 0.30;

    settings.learningRate = 0.80;
    settings.innerDistanceRatio = 0.1;
    settings.outerDistanceRatio = 0.1;

    settings.initialOuterProposalMultiplier = 0.01;
    settings.initialInnerProposalMultiplier = 0.01;

    settings.showDebugFigure = false;
    settings.saveVideo = false;
end