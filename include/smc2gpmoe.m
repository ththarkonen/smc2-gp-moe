function [resultsObj] = smc2gpmoe( data, settings)
    
    J = settings.outerParticleAmount;
    N = settings.nInputDimensions;
    
    minMcmcSteps = settings.outerMinMcmcSteps;
    maxMcmcSteps = settings.outerMaxMcmcSteps;
    
    showDebugFigure = settings.showDebugFigure;
    saveVideo = settings.saveVideo;

    timeStamp = datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-SS' );
    timeStamp = string( timeStamp );
    timeStart = tic();
        
    disp('Starting SMC2...');

    if( showDebugFigure )

        fileName = ".\videos\" + timeStamp;
        writerObj = VideoWriter( fileName, 'Uncompressed AVI');
        writerObj.FrameRate = 4;

        % open the video writer
        if( saveVideo )
            open(writerObj);
        end
    end

    psiParticles = initializePsi( settings ); 
    [innerSMCs, outerLogLikelihoods, gatingProbabilities, logTarget] = initializeInnerSMCs( data, psiParticles, settings);

    eta = settings.learningRate;
    targetAR = settings.outerTargetAcceptanceRate;

    essHist = J;
    mcmcAcceptanceRateHist = targetAR;
    essAcceptanceRateHist = J;

    smcWeights = ones( J, 1) / J;

    smc2Object = struct();
    smc2Object.data = data;

    smc2Object.psiParticles = psiParticles;
    smc2Object.proposalMat = [];
    smc2Object.outerLL = outerLogLikelihoods;
    smc2Object.gatingProbabilities = gatingProbabilities;
    smc2Object.logTarget = logTarget;

    smc2Object.ESS = {};
    smc2Object.K = settings.numberOfExperts;
    smc2Object.acceptanceRates = {};

    smc2Object.kappaVector = [0];
    smc2Object.innerSmcObjects = innerSMCs;
    smc2Object.settings = settings;

    if( N == 1 && showDebugFigure )

        frame = showIterationResults( smc2Object, psiParticles);

        if( saveVideo )
            writeVideo( writerObj, frame);
        end
    end

    counter = 1;
    likelihoodEvaluations = 0;
    mcmcMP = settings.initialOuterProposalMultiplier;

    while(true)

        psiParticles = smc2Object.psiParticles;
        outerLogLikelihoods = smc2Object.outerLL;
        kappaVector = smc2Object.kappaVector;
        innerSMCs = smc2Object.innerSmcObjects;

        kappa = kappaVector(counter);
        essAR = essAcceptanceRateHist(counter);
        ar = mcmcAcceptanceRateHist(counter);

        [smcWeights, kappa, ESS] = weightTempering( smcWeights, outerLogLikelihoods, eta, essAR, kappa, psiParticles);

        disp(['Next kappa: ', num2str(kappa)]);
        disp(' ');

        kappaVector(end+1) = kappa;
        essHist(end+1) = ESS;

        smc2Object.kappaVector = kappaVector;
        
        tempESS = -1;
        [psiParticles, innerSMCs, outerLogLikelihoods, smcWeights, mcmcProposalMat, ~] = outerProposalAndResample( psiParticles, innerSMCs, outerLogLikelihoods, smcWeights, tempESS, ar, mcmcMP, settings);

        update = progressSingleSMC( J, " ");
        tic
        for jj = 1:J

            innerSmcObject_jj = innerSMCs{jj};
            innerSmcObject_jj.kappaVector = kappaVector;

            innerSmcObject_jj = singleStepInnerSMC( innerSmcObject_jj, settings);

            outerLogLikelihoods(jj) = innerSmcObject_jj.outerLogLikelihood;
            innerSMCs{jj} = innerSmcObject_jj;

            update();
        end
        toc
        smc2Object.psiParticles = psiParticles;
        smc2Object.proposalMat = mcmcProposalMat;
        smc2Object.innerSmcObjects = innerSMCs;
        smc2Object.outerLL = outerLogLikelihoods;

        amountAccepted = 0;
        mcmcDistance = 0;

        flagStopMCMC = false;
        mcmcCounter = 0;
        prevPsiParticles = psiParticles;
        prevLogTarget = smc2Object.logTarget;

        while( ~flagStopMCMC || (mcmcCounter < minMcmcSteps) )

            mcmcCounter = mcmcCounter + 1;

            [smc2Object, tempAccepted] = psiMcmcUpdate( smc2Object );
            amountAccepted = amountAccepted + tempAccepted;

            [flagStopMCMC, mcmcDistance] = stopOuterMCMC( smc2Object, prevLogTarget, mcmcDistance);

            acceptanceRate = tempAccepted / (J);
            ar = acceptanceRate;

            psiParticles = smc2Object.psiParticles;
            outerLogLikelihoods = smc2Object.outerLL;
            innerSMCs = smc2Object.innerSmcObjects;

            if( mcmcCounter > maxMcmcSteps )
                flagStopMCMC = true;
            end

            [~, ~, ~, ~, mcmcProposalMat, mcmcMP] = outerProposalAndResample( psiParticles, innerSMCs, outerLogLikelihoods, smcWeights, ESS, ar, mcmcMP, settings);

            smc2Object.psiParticles = psiParticles;
            smc2Object.outerLL = outerLogLikelihoods;
            smc2Object.innerSmcObjects = innerSMCs;
            smc2Object.proposalMat = mcmcProposalMat;

            disp(['Iteration acceptance rate: ', num2str(acceptanceRate), ' with ',num2str(tempAccepted),' accepted of ', num2str(J), ' ', num2str(kappa)]);

            if( N == 1 && showDebugFigure )

                frame = showIterationResults( smc2Object, prevPsiParticles);

                if( saveVideo )
                    writeVideo( writerObj, frame);
                end
            end
            
        end

        acceptanceRate = amountAccepted / (J * mcmcCounter);
        ESS = computeESS( psiParticles, smcWeights);

        disp(' ');
        disp(['Iteration: ', num2str(length(kappaVector))]);
        disp(['Acceptance rate: ', num2str(acceptanceRate)]);
        disp(['Outer MCMC steps: ', num2str(mcmcCounter)]);
        disp(['Kappa: ', num2str(kappa)]);

        counter = counter + 1;

        mcmcAcceptanceRateHist(counter) = acceptanceRate;
        essAcceptanceRateHist(counter) = ESS;

        smc2Object.ESS = essAcceptanceRateHist;
        smc2Object.acceptanceRates = mcmcAcceptanceRateHist;

        if( kappa < 1 && acceptanceRate < (1/J) )
            disp('Low MH acceptance rate');
        end

        if( kappa >= 1 )
            break;
        end
    end

    if range(smcWeights) ~= 0
        
        [psiParticles, innerSMCs, outerLogLikelihoods, ~, mcmcProposalMat, ~] = outerProposalAndResample( psiParticles, innerSMCs, outerLogLikelihoods, smcWeights, tempESS, ar, mcmcMP, settings);
    
        smc2Object.psiParticles = psiParticles;
        smc2Object.outerLL = outerLogLikelihoods;
        smc2Object.innerSmcObjects = innerSMCs;
        smc2Object.proposalMat = mcmcProposalMat;
    end

    % close the writer object
    if( saveVideo )
        close( writerObj );
    end
    
    wallTime = toc(timeStart);
    for ii = 1:J

        innerSMC_jj = smc2Object.innerSmcObjects{ii};
        likelihoodEvaluations = likelihoodEvaluations + innerSMC_jj.likelihoodEvaluations;
    end


    resultsObj = smc2Object;
    resultsObj.wallTimeSeconds = wallTime;
end






