function [pi_kk] = sampleGatingPi( J, kk, K, settings)

    alpha = settings.prior.psi.concentrationParameter;
    pi_kk = drchrnd( alpha / K * ones( 1, K), J);
end