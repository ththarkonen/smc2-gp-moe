function [ll] = computeGpLogLikelihood( f, R, theta)
    
    nDims = size( R, 3);

    lengthScaleInds = 3:( 3 + nDims - 1 );
    
    m = theta(1);
    s2 = theta(2).^2;
    ls = theta( lengthScaleInds );
    s2F = theta(end).^2;
    
    f = f - m;
    nData = length(f);
    
    s2I = s2 * eye( nData );
    K = s2F * sqrExpCovMatrix( R, ls);

    try
        warning ( 'off', 'all');
        ll = gpLogLikelihood( f, K + s2I);
        warning ( 'on', 'all');
    catch
        ll = -Inf;
    end
end