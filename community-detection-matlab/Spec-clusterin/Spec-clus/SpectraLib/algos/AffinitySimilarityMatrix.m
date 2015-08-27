function S=AffinitySimilarityMatrix(points,sigma) 

% function S=DistanceSimilarityMatrix(points,sigma) 
% Returns the similarity between two points based on the distance between
% the two points. NOTE: The vectors are assumed to be COLUMN VECTORS 
  
  exp_mult_factor=-0.5/(sigma*sigma); 
  num_vectors=size(points,2); 

  S=squareform(pdist(points')); % need to transpose it pdist expects row vectors
  S=exp(exp_mult_factor*S);   
  
%   S=zeros(num_vectors,num_vectors); 
%   for i=1:num_vectors
% 	for j=i:num_vectors
%         diff=points(:,i)-points(:,j);        
% 	    S(i,j)=exp( exp_mult_factor*sum(diff.^2) );
%         S(j,i)=S(i,j); 
% 	end
%   end
  
