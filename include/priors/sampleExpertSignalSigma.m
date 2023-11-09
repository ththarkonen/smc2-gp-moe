function [expertSignalSigmas] = sampleExpertSignalSigma( J, kk, K, settings)
    
    expertSignalSigmas = settings.prior.theta.signalSigmaSigma * randn( J, 1);
    expertSignalSigmas = abs( expertSignalSigmas );
end