function [inds_kk] = computeThetaInds( kk, K, nInputDimensions)

    mu_k_Ind = kk;
    noiseSigma_k_Ind = kk + K;
    lengthScales_k_Inds = (1 + (kk - 1) * nInputDimensions + 2 * K):(1 + (kk - 1) * nInputDimensions + 2 * K + nInputDimensions - 1);
    signalSigma_k_Ind = (kk + 2 * K + K * nInputDimensions);

    inds_kk = [ mu_k_Ind, noiseSigma_k_Ind, lengthScales_k_Inds, signalSigma_k_Ind];
end