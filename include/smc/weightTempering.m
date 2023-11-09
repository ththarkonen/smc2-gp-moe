function [temperedWeights, kappa, ESS] = weightTempering( w, ll, learningRate, essAcceptanceRate, kappaPrev, particles )

    minKappa = kappaPrev;
    maxKappa = 1.0;
    tempKappa = 1.0;  
    
    w(isnan(w)) = eps;
    w(w == 0) = eps;
    
    tempWeights = exp( ( tempKappa - kappaPrev ) * ( ll - max(ll) ) );
    tempWeights = tempWeights / sum(tempWeights);

    tempESS = computeESS( particles, tempWeights);
    
    alphaESS = learningRate * essAcceptanceRate;
    absDiffESS = abs( tempESS - alphaESS );
    
    if( tempESS < alphaESS )
        while( absDiffESS > 1 )
            
            if( tempESS < alphaESS )
                maxKappa = tempKappa;
            else
                minKappa = tempKappa;
            end
            
            tempKappa = 0.5 * ( maxKappa + minKappa );
            
            tempWeights = exp( ( tempKappa - kappaPrev ) * ( ll - max(ll) ) );
            tempWeights = tempWeights / sum(tempWeights);

            tempESS = computeESS( particles, tempWeights);
            absDiffESS = abs( tempESS - alphaESS );
        end
    end
    
    temperedWeights = tempWeights;
    kappa = tempKappa;
    ESS = tempESS;
end