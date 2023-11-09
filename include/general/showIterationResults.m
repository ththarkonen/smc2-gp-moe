function [frame] = showIterationResults( smc2Object, prevPsiParticles)

    blue = [52, 121, 211] / 255;
    fuchsia = [211, 52, 121] / 255;
    green = [121, 211, 52] / 255;
    black = [170 170 170] / 255;
    
    colors = [ blue; fuchsia; green; black];
    nColors = size( colors, 1);

    data = smc2Object.data;
    
    K = smc2Object.K;
    psi = smc2Object.psiParticles;
    
    innerSmcObjects_K = smc2Object.innerSmcObjects;
    settings = smc2Object.settings;
    
    mu0 = linspace( 0, 1, K + 2);
    mu0 = mu0( 2:(end-1) );
    dMu = mu0(2) - mu0(1);
    
    gatingNetworkMuSigma = settings.prior.psi.muSigma * dMu;
    gatingNetworkSigmaSigma = settings.prior.psi.sigmaSigma * dMu;
    
    [ J, nParameters] = size(psi);
    [ M, ~] = size( innerSmcObjects_K{1}.thetaParticles );
    
    nGpParameters = 4*K;
       
    gpParameters = zeros( J * M, nGpParameters);

    expertsUsedVector = zeros( K, 1);
    
    for jj = 1:J
        
        startInd = 1 + ( jj - 1 ) * M;
        endInd = jj * M;
        inds = startInd:endInd;

        tempObj = innerSmcObjects_K{jj};
        
        nExperts = 0;
        
        for kk = 1:K 
        
            C_kk = tempObj.C{kk};
            
            if( ~isempty( C_kk ) )
                nExperts = nExperts + 1;
            end
        end
        
        expertsUsedVector(nExperts) = expertsUsedVector(nExperts) + 1;
        
        gpParameters( inds, :) = tempObj.thetaParticles;
    end
    
    expertsUsedVector = expertsUsedVector / sum(expertsUsedVector);
    
    muInds = 1:K;
    sigmaInds = ( K + 1 ):2*K;
    piInds = ( 2*K + 1 ):3*K;
    
    meanPsi = median(psi);
    
    meanMu = meanPsi(muInds);
    meanSigma = meanPsi(sigmaInds);
    meanPi = meanPsi(piInds);
    
    meanTheta = median( gpParameters );
    
    meanGpMu = meanTheta(1:K);
    meanGpSigma = meanTheta( (K+1):(2*K) );
    meanGpLS = meanTheta( ( 2*K + 1 ):(3*K) );
    meanGpSigmaF = meanTheta( ( 3*K + 1 ):(4*K) );

    [partedData, ~, ~, errorFlag] = partition( data, meanMu, meanSigma, meanPi);
                                                                 
    psiTitles = {};
    gpTitles = {};
    
    counter = 0;
    for kk = 1:K
        
        counter = counter + 1;
        
        meanInd = counter;
        sigmaInd = K + counter;
        piInd = 2*K + counter;
        
        gpMeanInd = counter;
        gpSigmaInd = K + counter;
        gpLengthInd = 2*K + counter;
        gpSignalSigmaInd = 3*K + counter;
        
        meanStr_k = ['$\mu_{',num2str(kk),'}$'];
        sigmaStr_k = ['$\sigma_{',num2str(kk),'}$'];
        piStr_k = ['$v_{',num2str(kk),'}$'];
        
        gpMeanStr_k = ['$m_{',num2str(kk),'}$'];
        gpSigmaStr_k = ['$\sigma_{',num2str(kk),'}$'];
        gpLengthStr_k = ['$l_{',num2str(kk),'}$'];
        gpSignalSigmaStr_k = ['$\sigma_{f,',num2str(kk),'}$'];
        
        psiTitles{meanInd} = meanStr_k;
        psiTitles{sigmaInd} = sigmaStr_k;
        psiTitles{piInd} = piStr_k;
        
        gpTitles{gpMeanInd} = gpMeanStr_k;
        gpTitles{gpSigmaInd} = gpSigmaStr_k;
        gpTitles{gpLengthInd} = gpLengthStr_k;
        gpTitles{gpSignalSigmaInd} = gpSignalSigmaStr_k;
    end
    
    hFigure = figure(1);
    clf
    
    hFigure.Name = '';
    hFigure.Units = 'normalized';
    hFigure.OuterPosition = [0 0 1 1];
    
    rootPanel = uix.GridFlex( 'Parent', hFigure, 'Spacing', 5, 'Padding', 5 );
    
    dataPanel = uipanel( 'Parent', rootPanel, 'Title', '');
    thetaPanel = uipanel( 'Parent', rootPanel, 'Title', '');
    psiPanel = uipanel( 'Parent', rootPanel, 'Title', '');
    extraPanel = uipanel( 'Parent', rootPanel, 'Title', '');
    
    set( rootPanel, 'Widths', [-1 -1], 'Heights', [-1 -1]);
   
    counter = 1;
    counterParam = 1;
    
    for ii = 1:nParameters
        
        subplot( 3, K, ii, 'Parent', psiPanel);
        hold on
        
        nBins = 20;
        
        h = histogram( psi(:,ii), 'NumBins', nBins);
        
        colorInd = mod( counter, nColors);
        
        if( colorInd == 0)
            colorInd = nColors;
        end
        
        h.FaceColor = colors( colorInd, :);
        h.FaceAlpha = 1.0;
        
        if( counterParam == 1 )
            
            mu_kk = mu0(counter);
            s = gatingNetworkMuSigma;
            
            xTemp = linspace( mu_kk - 3*s, mu_kk + 3*s, 100);
            muTemp = normpdf( xTemp, mu_kk, s);
            
            yyaxis right
            hPrior = plot( xTemp, muTemp);
            yyaxis left
            
            hPrior.Color = 'black';
        end
        
        if( counterParam == 2 )
            
            s = gatingNetworkSigmaSigma;
            
            xTemp = linspace( 0, 3*s, 100);
            sigmaTemp = normpdf( xTemp, 0, s);
            
            yyaxis right
            hPrior = plot( xTemp, sigmaTemp);
            yyaxis left
            
            hPrior.Color = 'black';
        end
        
        h = title( psiTitles{ii}, 'interpreter', 'latex');
        h.FontSize = 15;
        
        axis tight
        
        counter = counter + 1;
        
        if( counter > K )
            counter = 1;
            counterParam = counterParam + 1;
        end
    end
    
    sgtitle('$\Psi$', 'interpreter', 'latex');
    
    counter = 1;
    for ii = 1:nGpParameters
        
        subplot( 4, K, ii, 'Parent', thetaPanel);
        hold on
        
        h = histogram( gpParameters(:,ii) );
        
        colorInd = mod( counter, nColors);
        
        if( colorInd == 0)
            colorInd = nColors;
        end
        
        h.FaceColor = colors( colorInd, :);
        h.FaceAlpha = 1.0;
        
        h = title( gpTitles{ii}, 'interpreter', 'latex');
        h.FontSize = 15;
        
        counter = counter + 1;
        
        if( counter > K )
            counter = 1;
        end
    end
    
    sgtitle('$\Theta$', 'interpreter', 'latex');
    
    subplot( 1, 1, 1, 'Parent', dataPanel);
    cla
    hold on
    
    for kk = 1:K
        
       	xy = partedData.C{kk};
        
        if( isempty(xy) )
            continue;
        end

        x_k = xy(:,1);
        y_k = xy(:,2);
        
        theta_k = [ meanGpMu(kk), meanGpSigma(kk), meanGpLS(kk), meanGpSigmaF(kk)];
        
        [yStar, yStarSigma] = createPredictions( x_k, y_k, theta_k); 
        
        areaX = [x_k; flipud(x_k)];
        areaY = [yStar + 2 * yStarSigma; flipud(yStar - 2 * yStarSigma)];
        
        colorInd = mod( kk, nColors);
        
        if( colorInd == 0)
            colorInd = nColors;
        end
        
        tempColor = colors( colorInd, :);
        
        h = plot( x_k, y_k, '.');
        h.Color = tempColor;
        h.MarkerSize = 25;
        
        hF = fill( areaX, areaY, tempColor);
        hF.FaceAlpha = 0.15;
        hF.LineStyle = 'none';
        
        h = plot( x_k, yStar);
        h.Color = tempColor;
        h.LineWidth = 2;
    end
    
    psi0 = prevPsiParticles;
    X = randn( J, 1);
    
    if( ~isempty(psi0) )
        
        counter = 1;
        
        for ii = 1:nParameters
        
            colorInd = mod( counter, nColors);

            if( colorInd == 0)
                colorInd = nColors;
            end

            tempColor = colors( colorInd, :);

            subplot( 4, K, ii, 'Parent', extraPanel);
            h = plot( psi(:,ii),  psi0(:,ii), '.');
            h.Color = tempColor;
            axis equal

            h = title( psiTitles{ii}, 'interpreter', 'latex');
            h.FontSize = 15;

            counter = counter + 1;

            if( counter > K )
                counter = 1;
            end
                
        end
    end
    
    subplot( 4, K, 4*K)
    bar( expertsUsedVector )
    h = title( '$K$', 'interpreter', 'latex');
    h.FontSize = 15;
    
    drawnow();
    frame = getframe(gcf); 
    drawnow();
end