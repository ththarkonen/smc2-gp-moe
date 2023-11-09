function [R] = computeDistanceMatrix( x )

    [ n, nInputDimensions] = size( x );

    R = zeros( n, n, nInputDimensions);

    for ii = 1:n
        for jj = 1:n
            for d = 1:nInputDimensions
            
                xII = x( ii, d);
                xJJ = x( jj, d);
                
                dist = abs( xII - xJJ );
                
                R( ii, jj, d) = dist;
                R( jj, ii, d) = dist;
            end
        end
    end
end