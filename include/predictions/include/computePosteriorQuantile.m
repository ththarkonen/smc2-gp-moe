function [posteriorMedian] = computePosteriorQuantile( Y, yDensity, p)

    ySum = cumsum( yDensity, 1);
    
    [nY, nX] = size( ySum );
    
    posteriorMedian = zeros( 1, nX);
    
    for ii = 1:nX
        
        tempIntgY = ySum(:,ii);
        
        d = abs( tempIntgY - p );
        
        [~, ind] = min( d );
        
        yValue = tempIntgY(ind);
        
        if( yValue > p )    
            secondInd = ind - 1;
        else
            secondInd = ind + 1;
        end
        
        yValueSecond = tempIntgY(secondInd);
        
        dy = abs( yValue - yValueSecond );
        y0 = min( [ yValue, yValueSecond] );
        ind0 = min( [ ind, secondInd] );
        
        dInd = ( p - y0 ) / dy;
        
        medianInd = ind0 + dInd;
        
        posteriorMedian(ii) = interp1( Y, medianInd);
    end
end