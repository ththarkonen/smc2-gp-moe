function [ highDensityMatrix ] = computeHighDensityArea( result )
    
    [ nY, nX] = size( result.densities.grids.X );
    highDensityMatrix = NaN( nY, nX);

    for ii = 1:nX

        C_ii = result.densities.data(:,ii);

        [cOpt, ~] = fminsearch( @(c) ( sum(C_ii( C_ii <= c )) - 0.10 )^2, 0.001);
        highDensityMatrix( C_ii > cOpt, ii) = 1;
    end
end