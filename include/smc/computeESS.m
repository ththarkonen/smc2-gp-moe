function [ESS] = computeESS( p, w)

%     uniqueParticles = unique( p, 'rows' );
%     nUniques = length(uniqueParticles);
% 
%     tempWeightsUniques = zeros( nUniques, 1);
% 
%     for ii = 1:nUniques
% 
%         tempParticle = uniqueParticles(ii);
%         indsDuplicates = logical( p(:,1) == tempParticle(:,1) );
% 
%         tempWeightsDuplicates = w( indsDuplicates );
%         tempWeightsUniques(ii) = sum( tempWeightsDuplicates );
%     end
% 
%     ESS = 1 / sum( tempWeightsUniques.^2 );
    ESS = 1 / sum( w.^2 );
end