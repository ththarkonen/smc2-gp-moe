function [expertNoiseSigmas] = sampleExpertNoiseSigma( J, kk, K, settings)
    
    expertNoiseSigmas = settings.prior.theta.noiseSigmaSigma * randn( J, 1);
    expertNoiseSigmas = abs( expertNoiseSigmas );
end