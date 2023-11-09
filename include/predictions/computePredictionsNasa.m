function [ smc2Object ] = computePredictionsNasa( smc2Object, settings)

    psi = smc2Object.psiParticles;
    
    gatingProbabilities = smc2Object.gatingProbabilities;
    innerSMCs = smc2Object.innerSmcObjects;
    
    N = smc2Object.settings.nInputDimensions;
    K = smc2Object.K;
    NK = N * K;
    
    xData = smc2Object.data.xSubset;
    yData = smc2Object.data.ySubset;
    
    nX = settings.nPredictionLocations;
    nY = settings.nDensityLocations;

    xMin = min( xData(:,1) );
    xMax = max( xData(:,1) );

    x2 = xData(1,2);
    x3 = xData(1,3);

    x2Temp = repmat( x2, [ nX, 1]);
    x3Temp = repmat( x3, [ nX, 1]);

    yMin = min( yData );
    yMax = max( yData );
    
    yInterval = yMax - yMin;
    yPad = settings.yPadding * yInterval;
    
    xStar = linspace( xMin, xMax, nX)';
    yStar = linspace( yMin - yPad, yMax + yPad, nY)';

    xStarMat = [xStar, x2Temp, x3Temp];

    data = {};
    data.x = xStarMat;
    data.y = yStar;
    data.R = computeDistanceMatrix( xStarMat );
    
    yDataGrid = zeros( nY, nX);
    yGrid = zeros( nY, nX);

    J = smc2Object.settings.outerParticleAmount;
    M = smc2Object.settings.innerParticleAmount;

    for jj = 1:J
        
        psi_jj = psi( jj, :);
        innerSMC_jj = innerSMCs{jj};
        theta_jj = innerSMC_jj.thetaParticles;
        
        muInds = 1:( NK );
        sigmaInds = ( NK + 1 ):2*NK;
        piInds = ( 2*NK + 1 ):(2*NK + K);
        
        mu_jj = psi_jj( muInds );
        sigma_jj = psi_jj( sigmaInds );
        pi_jj = psi_jj( piInds );
        
        C_jj = innerSMC_jj.C;
        gatingProbability_jj = gatingProbabilities(jj);

        [~, ~, classProbability, ~] = partition( data, mu_jj, sigma_jj, pi_jj);
        
        for kk = 1:K
           
            C_jj_kk = C_jj{kk};
        
            if( isempty( C_jj_kk ) )
                continue;
            end
            
            x_kk = C_jj_kk( :, 1:3);
            y_kk = C_jj_kk( :, 4);
            
            classProbability_kk = classProbability(:,kk);
        
            thetaInds_kk = computeThetaInds( kk, K, N);
            
            for m = 1:M
            
                theta_m = theta_jj( m, thetaInds_kk);
                measurementSigma = theta_m(2);
                measurementVariance = measurementSigma.^2;
                
                [yStar_m, yStarSigma_m] = createResultPredictions( x_kk, y_kk, xStarMat, theta_m);
            
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
        
        progress = jj / J
        disp('asd')
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