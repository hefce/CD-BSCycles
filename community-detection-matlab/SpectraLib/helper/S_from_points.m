function S=S_from_points(points,sigma,smoothing_constant,self_similarity_zero)
% function S=S_from_points(points,sigma,smoothing_constant,self_similarity_zero)

  % calculate the similarity in case sigma > 0 
  if sigma > 0 
  	similarity=AffinitySimilarityMatrix(points,sigma) ; 
  end

  % smooth the similarity in case the smoothing_constant > 0 
  if smoothing_constant > 0 
	similarity=smooth_similarity(similarity,smoothing_constant); 
  end 
  
  if nargin > 3
	if self_similarity_zero > 0 
	  similarity=makeDiagonalZero(similarity); 
	end
  end
  S=similarity; 