function [partedData, gatingProbability, classProbability, errorFlag] = partition( data, mu, sigma, pi_kk)

    x = data.x;
    y = data.y;
    R = data.R;

    n = size( x, 1);
    nDimensions = size( x, 2);
    
    K = length(mu) / nDimensions;
    errorFlag = false;

    C_kk = cell(1,K);
    R_kk = cell(1,K);

    pi_kk = pi_kk(:)';
    pi_kk = pi_kk / sum( pi_kk );
    
    pi_kk = repmat( pi_kk, n, 1);
    
    % Covariance matrix

    logClassProbability = zeros( n, K);
    
    mu = mu';
    
    startInds = 1:nDimensions:( K * nDimensions );
    
    for kk = 1:K
        
        startInd_kk = startInds(kk);
        
        inds = startInd_kk:( startInd_kk + nDimensions - 1 );
        
        mu_k = mu(inds);
        sigma_k = sigma(inds);
        
        % Covariance matrix
        invC = diag( sigma_k.^(-2) );
        detC = prod( sigma_k.^2 );
        
        % Done using log for numerical stability
        for ii = 1:n

            x_ii = x( ii, :)';
            logClassProbability( ii, kk) = -0.5 * ( ( x_ii - mu_k )' * invC * ( x_ii - mu_k ) );
        end
    
        logClassProbability(:,kk) = logClassProbability(:,kk) - 0.5 * nDimensions * log( ( 2*pi ) ) - 0.5 * log( detC );
    end
    
    maxLCP = max( logClassProbability, [], 2);
    maxLCP = repmat( maxLCP, 1, K);

    classProbability = exp( logClassProbability - maxLCP );
    classProbability = pi_kk .* classProbability;
    
    rowSum = sum(classProbability, 2);
    rowSums = repmat( rowSum, 1, K);
    classProbability = classProbability ./ rowSums;
    
    upperLimits = cumsum(classProbability, 2);
    lowerLimits = [zeros( n, 1), upperLimits(:,1:K-1)];
    
    randomClusterer = rand( n, 1);
    randomClusterer = repmat( randomClusterer, 1, K);
    
    clusterInds = (randomClusterer > lowerLimits) & (randomClusterer <= upperLimits);
    
    gatingPointWiseProbability = classProbability(clusterInds);
    gatingProbability = prod(gatingPointWiseProbability);
    
    for kk = 1:K
        
        kInds = clusterInds(:,kk);
        C_kk{kk} = [x(kInds,:), y(kInds)];
        R_kk{kk} = R( kInds, kInds);
    end

    partedData.C = C_kk;
    partedData.R_K = R_kk;
end