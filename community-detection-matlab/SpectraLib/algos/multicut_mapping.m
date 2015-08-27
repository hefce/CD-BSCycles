function [vv]=multicut_mapping(S,k)
% function vv=multicut_mapping(S,k)
% Maps the similarity to the n points of dimension k-1 using the multicut
% lemma. 
  
%  global CLUSTER_MULTICUT_NORMALIZE_EV

  n=size(S,2);    

  % NOTE THIS IS VERY IMPORTANT. I AM NOT SURE WHETHER MATHAMATICALLY
  % THIS IS THE CORRECT THING TO DO. THE RESULTS TEND TO GET MUCH BETTER
  % THIS WAY. 

  % S=makeDiagonalZero(S); 

  %% CHANGED.....
  % Compute the Stochastic Matrix  P
  % P=S_to_P(S);
  %  D=diag(sum(S)); 
  %  P=D\S; 
  
  D=diag(sum(S,2));
  % Compute the top k eigenvectors
  %  [ vv indices lambda_save ] = eigsort_gen( S , D, k ); 
  % lambda=lambda_save(indices(1:k)); 

  % use eigs instead of eig() 
  [ vv lambda] = myeigs_gen( S , D, k );

  %lambda=diag(lambda);

  %lambda = lambda( 2:k );              % keep only top k  
  % fprintf('mutlicut first EV mean = %g, std= %g\n', mean(vv(:,1)), std(vv(:,1))); 
%  if use_all
	vv = vv( :, 1:k ) ; %% Now vv is nxk
%  else 
%	vv = vv( :, 2:k ) ; %% Now vv is nx(k-1)
%  end
  
  vv=vv'; % makes it kxn i.e. n column vectors of dimension k
  
