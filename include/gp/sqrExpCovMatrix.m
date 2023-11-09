function [K] = sqrExpCovMatrix( R, lengthScales)

    [n, m, nDims] = size( R );
    K = zeros( n, m);
    
    for d = 1:nDims
        
        R_d = R(:,:,d);
        ls_d = lengthScales(d);
    
        K = K + exp( -0.5 * ( R_d / ls_d ).^2 );
    end
end