function [expertsUsedVector] = computePosteriorNumberOfExperts( smc2Object )

    K = smc2Object.K;
    innerSmcObjects_K = smc2Object.innerSmcObjects;
    
    [ J, ~] = size(innerSmcObjects_K)
    [ M, ~] = size( innerSmcObjects_K{1}.thetaParticles )

    expertsUsedVector = zeros( K, 1);
    
    for jj = 1:J
        
        startInd = 1 + ( jj - 1 ) * M;
        endInd = jj * M;
        inds = startInd:endInd;

        tempObj = innerSmcObjects_K{jj};
        
        nExperts = 0;
        
        for kk = 1:K 
        
            C_kk = tempObj.C{kk};
            
            if( ~isempty( C_kk ) )
                nExperts = nExperts + 1;
            end
        end
        
        expertsUsedVector(nExperts) = expertsUsedVector(nExperts) + 1;
        
        gpParameters( inds, :) = tempObj.thetaParticles;
    end
    
    expertsUsedVector = expertsUsedVector / sum(expertsUsedVector);
%     figure();
%     
%     bar( expertsUsedVector )
%     h = title( '$K$', 'interpreter', 'latex');
%     h.FontSize = 15;
end