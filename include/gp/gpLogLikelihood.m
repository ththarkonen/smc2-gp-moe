function [ll] = gpLogLikelihood( f, K)

    n = length(f);
    ll = -0.5 * ( f'*( K \ f ) + logdet(K) + n * log(2*pi) );
end