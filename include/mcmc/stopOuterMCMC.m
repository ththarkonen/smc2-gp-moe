function [stopFlag, distanceMean] = stopOuterMCMC( smc2Object, prevParticles, prevDistanceMean)

    particles = smc2Object.logTarget_K;
    
    dr = smc2Object.settings.outerDistanceRatio;
    
    residuals2 = ( particles - prevParticles ).^2;
    distances2 = sum( residuals2, 2);
    
    distances = sqrt( distances2 );
    distanceMean = mean(distances);
    
    distanceChange = abs(distanceMean - prevDistanceMean);
    
    relativeDistanceChange = distanceChange / prevDistanceMean;
    
    if( relativeDistanceChange < dr  )
        stopFlag = true;
    else
        stopFlag = false;
    end
    
end
