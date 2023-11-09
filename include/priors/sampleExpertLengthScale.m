function [expertLengthScales] = sampleExpertLengthScale( J, kk, K, settings)
    
    expertLengthScales = settings.prior.theta.lengthScaleSigma * randn( J, 1);
    expertLengthScales = abs( expertLengthScales );
end