function [K] = sqrExpCovMatrixFull( R, lengthScale)

    K =  exp( -R.^2 / lengthScale );
end