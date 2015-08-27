function assignment=cluster_spectral_general(S,k,mapping_method,assign_method)
% function assignment=cluster_spectral_general(S,k,mapping_method,assign_method)
% 
% Generic method for clustering. note that "S" can either be the
% similarity matrix or the points depending on the input argument. 
  xx=''; 
  switch (mapping_method) 
%    case {'ang_classic'}
% 	xx=ang_mapping(S,k,0); 
   case {'ang_gen'}
	xx=ang_mapping(S,k); 
%    case {'mcut_one_less'}
% 	xx=multicut_mapping(S,k,0);  
   case {'mcut_all'}
	xx=multicut_mapping(S,k);  
   case {'none'}
	xx=''; 
   otherwise
	warning(' invalid mapping method in cluster_spectral_general'); 
  end
  if length(xx) > 0 
	global SPECTRAL_SIGMA; 
	S=AffinitySimilarityMatrix(xx,SPECTRAL_SIGMA);
  end
  
  % Now S is either the original S or the new mapped S. In either casewe
  % can apply the new algorithm to get the assignment 
  
  global KMEANS_THRESHOLD KMEANS_MAX_ITER;
  threshold = KMEANS_THRESHOLD;
  max_iter  = KMEANS_MAX_ITER; 
  switch (assign_method )
   
	% four methods based on points in the space (use xx) 
   case {'anchor'}
	global ANCHOR_NMIN;
	assignment=cluster_point_anchor(xx,k,ANCHOR_NMIN);
	
   case {'kmeans'}
	assignment=cluster_point_kmeans(xx,k,5,20); 
   
   case {'kmeans_ortho'}
	assignment=cluster_point_kmeans(xx,k,1,0); 

   case {'kmeans_random'}
	assignment=cluster_point_kmeans(xx,k,0,20); 

	% hiearchical clustering based on  S 
   case {'ward'}
	assignment=cluster_ward_linkage(S,k); 
	
   case {'single'}
	assignment=cluster_single_linkage(S,k); 
	

	% two level spectral methods 	
   case {'ang_ward'}
	assignment=cluster_spectral_general(S,k,'ang_gen','ward');

   case {'ang_kmeans'}
	assignment=cluster_spectral_general(S,k,'ang_gen','kmeans');

   case {'mcut_ward'}
	assignment=cluster_spectral_general(S,k,'mcut_all','ward');

   case {'mcut_kmeans'}
	assignment=cluster_spectral_general(S,k,'mcut_all','kmeans');
	
   otherwise
	warning(' *** invalid assignment *** in cluster_spectral_general'); 
	
  end
  
	
  
