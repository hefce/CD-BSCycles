function S=AffinitySimilarityMatrix(points,sigma) 

% function S=DistanceSimilarityMatrix(points,sigma) 
% Returns the similarity between two points based on the distance between
% the two points. NOTE: The vectors are assumed to be COLUMN VECTORS 

  S=zeros(size(points,2)); 
  exp_mult_factor=-(0.5)./(sigma*sigma); 
  S=exp(exp_mult_factor*squareform(pdist(points')));   %points' are ROW VECTORS
  
%   num_vectors=size(points,2);
%   for i=1:num_vectors
% 	for j=i:num_vectors
%         diff=points(:,i)-points(:,j);        
% 	    S(i,j)=exp( exp_mult_factor*sum(diff.^2) );
%         S(j,i)=S(i,j); 
% 	end
%   end
  
