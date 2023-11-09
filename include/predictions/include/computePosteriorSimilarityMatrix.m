function [posteriorSimilarityMatrix] = computePosteriorSimilarityMatrix( smc2Object )

    xData = smc2Object.data.x;
    yData = smc2Object.data.y;
    
    K = smc2Object.K;
    
    data = [xData, yData];
    nData = size( data, 1);
    
    posteriorSimilarityMatrix = zeros( nData, nData);
    
    innerSmcObjects = smc2Object.innerSmcObjects;

    J = length( innerSmcObjects );
    
    for jj = 1:J
       
        innerSmcObject_jj = innerSmcObjects{jj};
        C_jj = innerSmcObject_jj.C;
        
        for kk = 1:K

            C_jj_kk = C_jj{kk};

            if( isempty( C_jj_kk ) )
                continue;
            end
            
            nTemp = size( C_jj_kk, 1);
            
            for n = 1:nTemp
                for m = 1:nTemp
                    
                    xyTempN = C_jj_kk( n, :);
                    xyTempM = C_jj_kk( m, :);

                    indN = find( ismember( data, xyTempN, 'rows') );
                    indM = find( ismember( data, xyTempM, 'rows') );
                    
                    posteriorSimilarityMatrix( indN, indM) = posteriorSimilarityMatrix( indN, indM) + 1;
                end
            end
        end
        
        progress = jj / J
    end
    
    posteriorSimilarityMatrix = posteriorSimilarityMatrix / J;
end