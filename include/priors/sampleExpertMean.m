function [expertMeans] = sampleExpertMean( J, kk, K, settings)

    expertMeans = settings.prior.theta.maxY * rand( J, 1);
end