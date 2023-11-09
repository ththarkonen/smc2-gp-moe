function [ psi, ll, innerSmcObjects, weights] = resampleOuter( psi, ll, innerSmcObjects, weights)

    [nParticles, ~] = size(psi);
    weights( weights == 0 ) = eps;

    resampleInds = resampleResidual( weights );
        
    psi = psi( resampleInds, :);
    innerSmcObjects = innerSmcObjects( resampleInds );

    ll = ll( resampleInds );
    weights = ones( nParticles, 1) / nParticles;
end