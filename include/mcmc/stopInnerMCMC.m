function [stopFlag, distanceMean] = stopInnerMCMC( particles, prevParticles, prevDistanceMean, settings)

    dr = settings.innerDistanceRatio;
    
    residuals2 = ( exp(particles) - exp(prevParticles) ).^2;
    distances2 = sum( residuals2, 2);
    
    distances = sqrt( distances2 );
    distanceMean = mean(distances);
    
    distanceChange = abs(distanceMean - prevDistanceMean);
    
    if( distanceChange < dr * prevDistanceMean )
        stopFlag = true;
    else
        stopFlag = false;
    end
    
end