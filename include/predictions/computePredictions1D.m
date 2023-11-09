function [ smc2Object ] = computePredictions1D( smc2Object, settings)

    data = smc2Object.data;
    psi = smc2Object.psiParticles;
    
    gatingProbabilities = smc2Object.gatingProbabilities;
    innerSMCs = smc2Object.innerSmcObjects;
    
    N = smc2Object.settings.nInputDimensions;
    K = smc2Object.K;

    xMin = min( data.x );
    xMax = max( data.x );
    
    yMin = min( data.y );
    yMax = max( data.y );
    
    yInterval = yMax - yMin;
    
    nX = settings.nPredictionLocations;
    nY = settings.nDensityLocations;
    yPad = settings.yPadding * yInterval;
    
    xStar = linspace( xMin, xMax, nX)';
    yStar = linspace( yMin - yPad, yMax + yPad, nY)';

    dataStar = {};
    dataStar.x = xStar;
    dataStar.y = yStar;
    dataStar.R = computeDistanceMatrix( xStar );
    
    yDataGrid = zeros( nY, nX);
    yGrid = zeros( nY, nX);

    J = smc2Object.settings.outerParticleAmount;
    M = smc2Object.settings.innerParticleAmount;

    for jj = 1:J
        
        psi_jj = psi( jj, :);
        innerSmcObject_jj = innerSMCs{jj};
        thetaParticles_jj = innerSmcObject_jj.thetaParticles;
        
        muInds = 1:K;
        sigmaInds = (K+1):2*K;
        piInds = (2*K + 1):3*K;
        
        mu_jj = psi_jj( muInds );
        sigma_jj = psi_jj( sigmaInds );
        pi_jj = psi_jj( piInds );
        
        C_jj = innerSmcObject_jj.C;
        gatingProbability_jj = gatingProbabilities(jj);

        [~, ~, classProbability, ~] = partition( dataStar, mu_jj, sigma_jj, pi_jj);
        
        for kk = 1:K
           
            C_jj_kk = C_jj{kk};
        
            if( isempty( C_jj_kk ) )
                continue;
            end
            
            x_kk = C_jj_kk(:,1);
            y_kk = C_jj_kk(:,2);
            
            classProbability_kk = classProbability(:,kk);
        
            thetaInds_kk = computeThetaInds( kk, K, N);
            
            for m = 1:M
            
                theta_m = thetaParticles_jj( m, thetaInds_kk);
                measurementSigma = theta_m(2);
                measurementVariance = measurementSigma.^2;
                
                [yStar_m, yStarSigma_m] = createResultPredictions( x_kk, y_kk, xStar, theta_m);
            
                parfor ii = 1:nX

                    mean_ii = yStar_m(ii);
                    sigma_ii = yStarSigma_m(ii);
                    dataPredictiveSigma_ii = sqrt( measurementVariance + sigma_ii.^2 );

                    yPredictive_ii = normpdf( yStar, mean_ii, sigma_ii);
                    yDataPredictive_ii = normpdf( yStar, mean_ii, dataPredictiveSigma_ii);
                    
                    yGrid(:,ii) = yGrid(:,ii) + gatingProbability_jj * classProbability_kk(ii) .* yPredictive_ii;
                    yDataGrid(:,ii) = yDataGrid(:,ii) + gatingProbability_jj * classProbability_kk(ii) .* yDataPredictive_ii;
                end
            end
        end
        
        progress = jj / J;
        disp( progress )
    end
    
    yGridIntegral = sum( yGrid );
    yDataGridIntegral = sum( yDataGrid );
    
    for ii = 1:nX
        yGrid(:,ii) = yGrid(:,ii) / yGridIntegral(ii);
    end
    
    for ii = 1:nX
        yDataGrid(:,ii) = yDataGrid(:,ii) / yDataGridIntegral(ii);
    end

    [ X, Y] = meshgrid( xStar, yStar);
    
    PSM = computePosteriorSimilarityMatrix( smc2Object );
    posteriorK = computePosteriorNumberOfExperts( smc2Object );
    
    posteriorMedian = computePosteriorMedian( yStar, yGrid );
    posteriorLowerBound = computePosteriorQuantile( yStar, yDataGrid, 0.05);
    posteriorUpperBound = computePosteriorQuantile( yStar, yDataGrid, 0.95);

    posteriorDataMedian = computePosteriorMedian( yStar, yDataGrid );
    posteriorDataLowerBound = computePosteriorQuantile( yStar, yDataGrid, 0.05);
    posteriorDataUpperBound = computePosteriorQuantile( yStar, yDataGrid, 0.95);

    modelQuantiles = struct();
    modelQuantiles.median = posteriorMedian;
    modelQuantiles.lowerBound = posteriorLowerBound;
    modelQuantiles.upperBound = posteriorUpperBound;

    dataQuantiles = struct();
    dataQuantiles.median = posteriorDataMedian;
    dataQuantiles.lowerBound = posteriorDataLowerBound;
    dataQuantiles.upperBound = posteriorDataUpperBound;

    smc2Object.densities = struct();
    smc2Object.densities.model = yGrid;
    smc2Object.densities.data = yDataGrid;

    smc2Object.densities.grids = struct();
    smc2Object.densities.grids.X = X;
    smc2Object.densities.grids.Y = Y;

    smc2Object.posteriorSimilarityMatrix = PSM;
    smc2Object.posteriorNumberOfExperts = posteriorK;

    smc2Object.quantiles = struct();
    smc2Object.quantiles.data = dataQuantiles;
    smc2Object.quantiles.model = modelQuantiles;
end